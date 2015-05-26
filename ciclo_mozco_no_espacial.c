/*
Copyright 2015 Jorge Velazquez
*/
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#include "libPP_5.0.h"
#include "EntSalArb_MP.h"
#include "GNA.h"
#include "MC_sweep_mozquito.h"
#include "conn_mysql.h"

main(){	
	
///////////////////////////Inicializa parametros de la simulacion
int NDX=50;
int NDY=NDX;
int T_max = 3250; //1400; 
int NoEnsambles=8;

int INI_FEMALE=0;
int INI_MALE=0;
float INI_DENSITY_PUPAE=0.3;

mozquito_parameters global_parameters;
global_parameters.PupaDeadRate=0.12;			
global_parameters.PupaOffspringRate=0.05;
global_parameters.FemaleOffspringFraction=0.5;
global_parameters.FemaleDeadRate=0.14*1.05042;
global_parameters.MaleDeadRate=0.14*1.05042;
global_parameters.FemaleOffspringRate=1.0;
global_parameters.Metabolic_Time = obtain_metabolic_time(&global_parameters);


omp_set_num_threads(4);

//date time of simulation
time_t raw_now = time(NULL);
struct tm * now;
now = localtime ( &raw_now );
char sim_time[50];
sprintf(sim_time,"'%d-%d-%d %d:%d:%d'",now->tm_year + 1900, now->tm_mon + 1, now->tm_mday,now->tm_hour,now->tm_min,now->tm_sec);
	
////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:

Float2D_MP MP_RhoVsT_1;		
InicializaFloat2D_MP(&MP_RhoVsT_1, T_max, 3, 0);
float FisicalTime[T_max+1];
FisicalTime[0]=0.0;
float R_dyn[T_max+1];
float R_Delta[T_max+1];
float R_Delta_2[T_max + 1];
float Temp_dyn[T_max+1];
			
			///////////////////////////////////// INICIA PARALLEL

			#pragma omp parallel			///////Estado INICIAL:
			{
				init_JKISS(); //Inicializa la semilla de cada proceso.
				
				int num_threads = omp_get_num_threads();
				int id = omp_get_thread_num();
				int MaxPar = NoEnsambles/num_threads;
				#pragma omp master
				{
						MaxPar+= NoEnsambles - MaxPar * num_threads;	 
				}
				estado e[MaxPar];
				mozquitos_state mozquitos[MaxPar];
				mozquito_parameters param;
				param = global_parameters;
				
				//MaxPar=CargaEstado_MP(contenedorLec,"T_001",e,NDX,NDY,id,MaxPar);
				
				int Par;
				for(Par=0;Par<MaxPar;Par++)
				{
					AlojaMemoria(&e[Par], NDX, NDY);  
					ResetEstado(&e[Par]);
					mozquitos[Par].female=INI_FEMALE;
					mozquitos[Par].male=INI_MALE;
				}

				
					for(Par=0;Par<MaxPar;Par++)
					{
					//InsertaIndividuosAleatorio(&e[Par],100,NoEspecie);
					GeneraEstadoAleatorio(&e[Par], INI_DENSITY_PUPAE, 1);
					}
				
			///////////////////////////////////////Termina Estado INICIAL
		Float2D_MP MP_RhoVsT;	
		InicializaFloat2D_MP(&MP_RhoVsT, T_max, 3, MaxPar);
		ResetFloat2D_MP(&MP_RhoVsT);
		MP_RhoVsT.NoEnsambles=MaxPar;
		float Area;
		float temperature, humidity;
		float R;
		temperature = calendar_temperature(50.0);
			
			////////////////////////////////Barrido Monte CARLO:
				int i;
				for(i=0;i<T_max;i++)
				{
					//temperature = calendar_temperature(FisicalTime[e[0].T]);
					humidity = calendar_humidity(FisicalTime[e[0].T]);
					set_param_temperature_dependent(&param,temperature);
					set_diapause_humidity_dependent(&param,humidity);
					for(Par=0;Par<MaxPar;Par++)
					{
						Area=(float)(e[Par].NDX*e[Par].NDY);
						MP_RhoVsT.array[e[Par].T][1]+=((float)e[Par].ON/Area);
						MP_RhoVsT.array[e[Par].T][2]+=((float)mozquitos[Par].female/Area);
						MP_RhoVsT.array[e[Par].T][3]+=((float)mozquitos[Par].male/Area);	
						MC_sweep_mozquito(&e[Par], &mozquitos[Par], &param);
					}
					
					
					#pragma omp single
					{
						FisicalTime[e[0].T]=FisicalTime[e[0].T-1] + 1.0/param.Metabolic_Time;
						if(i>2)
						{
						R=param.FemaleOffspringFraction*param.PupaOffspringRate*(param.FemaleOffspringRate + (calendar_humidity(FisicalTime[i]) - calendar_humidity(FisicalTime[i-2]))/(FisicalTime[i] - FisicalTime[i-2]));
						R=R/(param.FemaleDeadRate*(param.PupaDeadRate + param.PupaOffspringRate));			
						}else{R=0.0;}
						if(R>1.0)
						{
							R_dyn[i]=param.FemaleOffspringFraction*param.PupaOffspringRate*(R-1.0)/(R*param.FemaleDeadRate);
						}else{
							R_dyn[i]=0.0;
						}
						#pragma omp flush
						if(i>52)
						{
							//printf("Ri=%f , Ri-1=%f, DR =%f , DT = %f , D=%f , 1/eps= %f \n",R_dyn[e[0].T],R_dyn[e[0].T - 1],(R_dyn[i]-R_dyn[i-1]),(FisicalTime[i]-FisicalTime[i-1]),(R_dyn[i]-R_dyn[i-1])/(FisicalTime[i]-FisicalTime[i-1]),(1.0/param.FemaleDeadRate));	
							R_Delta[i]=((R_dyn[i]-R_dyn[i - 52])/(FisicalTime[i]-FisicalTime[i - 52]))/param.FemaleDeadRate;
						}else{ R_Delta[i]=0.0;}
						#pragma omp flush
						if(i>(52 + 100)) // 52 + 254
						{
							R_Delta_2[i]=((R_Delta[i]-R_Delta[i - 100])/(FisicalTime[i]-FisicalTime[i - 100]))/param.FemaleDeadRate;
						}else{ R_Delta_2[i]=0.0; }
						Temp_dyn[i]=temperature;
						//printf("i=%d, r=%f, r1=%f, r2=%f, rt=%f\n",i, R_dyn[i],R_Delta[i],R_Delta_2[i], R_dyn[i]-R_Delta[i]+R_Delta_2[i]);
					}
						if((i-(i/500)*500)==499)    //Inicializa cada 500 pasos
						{
							init_JKISS();
						}
		
				}
				
				for(Par=0;Par<MaxPar;Par++)
				{
					Area=(float)(e[Par].NDX*e[Par].NDY);
					MP_RhoVsT.array[e[Par].T][1]+=((float)e[Par].ON/Area);
					MP_RhoVsT.array[e[Par].T][2]+=((float)mozquitos[Par].female/Area);
					MP_RhoVsT.array[e[Par].T][3]+=((float)mozquitos[Par].male/Area);
				}
				
			////////////////////////////////Termina Monte CARLO
			#pragma omp barrier
			#pragma omp single
			{
				ResetFloat2D_MP(&MP_RhoVsT_1);
			}
			
				SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
				
			//	#pragma omp master
				//{
				//	PD_GuardaEstadoEn_MP(contenedorCompleto, e, id, 1);
				//}
									
						//Libera Memoria
						for(Par=0;Par<MaxPar;Par++)
						{
							LiberaMemoria(&e[Par]);
						}
						LiberaMemoriaFloat2D_MP(&MP_RhoVsT);	

			}	////////////////////////////////////////////////////////////////////TERMINA PARALLEL
int i;
		for(i=1;i<T_max;i++)
		{
			R_dyn[i]=R_dyn[i] - R_Delta[i] + R_Delta_2[i];
		}

	//// Guarda parametros en MySql	y crea CONTENEDOR
	char contenedor[300];
	sprintf(contenedor,"DependenciaEnTemperatura");
	CreaContenedor(contenedor);	
		
	char values[300];
	sprintf(values,"NULL,%s,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%d,%d,%d,%d,'%s',0",sim_time,
	global_parameters.PupaDeadRate,			
	global_parameters.PupaOffspringRate,
	global_parameters.FemaleOffspringFraction,
	global_parameters.FemaleDeadRate,
	global_parameters.MaleDeadRate,
	global_parameters.FemaleOffspringRate,
	global_parameters.Metabolic_Time,
	INI_FEMALE,
	INI_MALE,
	INI_DENSITY_PUPAE,
	NDX,
	NDY,
	T_max,
	NoEnsambles,
	contenedor
	);
	
	printf("MySQL client version: %s\n", mysql_get_client_info());
	MYSQL *con = connect_db("localhost","dengue_fcfm", "EPAFJV", "dengue_fcfm");
	int inserted_id = insert_into_db(con, "sim_non_spatial",values);
	mysql_close(con);
	char contenedorCompleto[200];
	sprintf(contenedorCompleto,"%s/%d",contenedor,inserted_id);
	CreaContenedor(contenedorCompleto);
	store_density_evolution(contenedorCompleto,&MP_RhoVsT_1, 0, FisicalTime,R_dyn,Temp_dyn);
	FILE *aA;
	char archivo[200];
	sprintf(archivo,"Graficas/simulation.tex");	
	aA=fopen(archivo, "w");
	fprintf(aA,"\\newcommand{\\data}{../%s/density_evolution}\n\\newcommand{\\plotTitle}{id=%d  }",contenedorCompleto,inserted_id);
	fclose(aA);
	sprintf(archivo,"Graficas/include_make");
	aA=fopen(archivo, "w");
	fprintf(aA,"id = %d\n",inserted_id);
	fclose(aA);
	system("cd Graficas; make figure");
	//////
	
	LiberaMemoriaFloat2D_MP(&MP_RhoVsT_1);

						
return;
}
