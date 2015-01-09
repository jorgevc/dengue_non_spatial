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
int NDX=300;
int NDY=NDX;
int T_max = 30;
int NoEnsambles=4;

int INI_FEMALE=0;
int INI_MALE=0;
float INI_DENSITY_PUPAE=0.3;

mozquito_parameters global_parameters;
global_parameters.PupaDeadRate=0.0;			
global_parameters.PupaOffspringRate=1.0;
global_parameters.FemaleOffspringFraction=0.5;
global_parameters.FemaleDeadRate=0.0;
global_parameters.MaleDeadRate=0.0;
global_parameters.FemaleOffspringRate=0.0;
global_parameters.Metabolic_Time = obtain_metabolic_time(&global_parameters);


omp_set_num_threads(4);
////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:

Float2D_MP MP_RhoVsT_1;		
InicializaFloat2D_MP(&MP_RhoVsT_1, T_max, 3, 0);		
			
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
		
			////////////////////////////////Barrido Monte CARLO:
				int i;
				for(i=0;i<T_max;i++)
				{
					for(Par=0;Par<MaxPar;Par++)
					{
						Area=(float)(e[Par].NDX*e[Par].NDY);
						MP_RhoVsT.array[e[Par].T][1]+=((float)e[Par].ON/Area);
						MP_RhoVsT.array[e[Par].T][2]+=((float)mozquitos[Par].female/Area);
						MP_RhoVsT.array[e[Par].T][3]+=((float)mozquitos[Par].male/Area);
						MC_sweep_mozquito(&e[Par], &mozquitos[Par], &param);
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

	
	//// Guarda parametros en MySql	y crea CONTENEDOR
	char contenedor[300];
	sprintf(contenedor,"TEST");
	CreaContenedor(contenedor);	
	
	time_t raw_now = time(NULL);
	struct tm * now;
	now = localtime ( &raw_now );
	char sim_time[50];
	sprintf(sim_time,"'%d-%d-%d %d:%d:%d'",now->tm_year + 1900, now->tm_mon + 1, now->tm_mday,now->tm_hour,now->tm_min,now->tm_sec);
		
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
	store_density_evolution(contenedorCompleto,&MP_RhoVsT_1, 0);	
	//////
	
	LiberaMemoriaFloat2D_MP(&MP_RhoVsT_1);

						
return;
}
