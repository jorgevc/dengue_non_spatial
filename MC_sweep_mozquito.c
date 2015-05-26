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

#include "MC_sweep_mozquito.h"
#include "GNA.h"
#include <stdio.h>
#include <time.h>

#define DIAPAUSIC 2
#define NODIAPAUSIC 1

void MC_sweep_mozquito(estado *es, mozquitos_state *mozquitos, mozquito_parameters *param)
{
//	start2 = clock();  //clock comentar!!
int Indice;
float DT=0.0;
int Total_Individuals=es->ON + mozquitos->female + mozquitos->male;
	
	while(DT<1.0){
		if(Total_Individuals != 0){
			DT+=1.0/Total_Individuals; 
			Indice = I_JKISS(1,Total_Individuals);
			Update_mozquito(es, mozquitos, param, Indice);
			Total_Individuals = es->ON + mozquitos->female + mozquitos->male;	
		}else{
			DT=2.0;
		}
	}
	
(es->T)++;		

//     end2 = clock();			//clock comentar!!
  //   tiempo2 = ((double) (end2 - start2))/ CLOCKS_PER_SEC;;  //clock comentar!!	
  //   printf("Tiempo en elejir vecinos en una corrida:%f, tiempo corrida:%f\n",cpu_time_used,tiempo2);
   //  cpu_time_used=0.0;
		
return;
}

void Update_mozquito(estado *es,mozquitos_state *mozquitos, mozquito_parameters *parameters, int N)
{
double Rand; 
float pPupaDead, pPupaOffspring, pMaleDead, pFemaleDead, pMozquitoOffspring, pPupaAntiDiapasue,pPupaDiapause;
sitio vecino;
int NDX = es->NDX;
int NDY = es->NDY;
int i,j;

	Rand = F_JKISS();
	
	if(N <= es->ON) //pupae
	{
		i=es->SO[N].i;
		j=es->SO[N].j;
		pPupaDead=parameters->PupaDeadRate/parameters->Metabolic_Time; //Asignar Max_Metabolic, si no hay division entre cero.
		pPupaOffspring=parameters->PupaOffspringRate/parameters->Metabolic_Time;
		pPupaDiapause=parameters->DiapauseRate/parameters->Metabolic_Time;
		
		if(es->TIPO[i][j] == DIAPAUSIC)
		{
			pPupaAntiDiapasue=parameters->AntiDiapauseRate/parameters->Metabolic_Time;
			if(Rand <= pPupaAntiDiapasue )
			{
				es->TIPO[i][j] = NODIAPAUSIC ;
			}
		}else{
			if(Rand<=(pPupaOffspring + pPupaDead )) //dead or offspring
			{	
				es->s[i][j]=0;
				es->SO[N]=es->SO[(es->ON)];
				es->INDICE[es->SO[es->ON].i][es->SO[es->ON].j]=N;
				(es->ON)--;	
				
				if(Rand <= pPupaOffspring )
				{
					if(Rand <= (pPupaOffspring*(parameters->FemaleOffspringFraction)))
					{
						(mozquitos->female)++;
					}else{
						(mozquitos->male)++;
					}
				}	
			}else{	//diapause or nothing
				if(Rand <= pPupaDiapause)	//diapause
				{
						es->TIPO[i][j]= DIAPAUSIC ;
				}
			}
		}
	}else{		//mozquitos
		if(N <= es->ON + mozquitos->female) //female mozquitos
		{
			pFemaleDead=parameters->FemaleDeadRate/parameters->Metabolic_Time;
			pMozquitoOffspring=parameters->FemaleOffspringRate/parameters->Metabolic_Time;
			
			if(Rand <= pFemaleDead)
			{
				(mozquitos->female)--;
			}else{
				if(Rand <= (pFemaleDead + pMozquitoOffspring))
				{
					vecino.i=I_JKISS(1,es->NDX);
					vecino.j=I_JKISS(1,es->NDY);
					
					if(es->s[vecino.i][vecino.j]==0)
					{
						es->s[vecino.i][vecino.j]=1;
						(es->ON)++;
						es->SO[(es->ON)]=vecino;
						es->INDICE[vecino.i][vecino.j]=(es->ON);
					}	
				}
			}	
		}else{		//male mozquitos
			pMaleDead=parameters->MaleDeadRate/parameters->Metabolic_Time;
			if(Rand <= pMaleDead)
			{
				(mozquitos->male)--;
			}
		}
	}
return;
}

float obtain_metabolic_time(mozquito_parameters *param)
{
	float metaTime,pupaTime,femaleTime,maleTime,diapausicTime;
	
	pupaTime = param->PupaDeadRate + param->PupaOffspringRate;
	femaleTime = param->FemaleDeadRate + param->FemaleOffspringRate;
	maleTime = param->MaleDeadRate;
	
	if(pupaTime > femaleTime)
	{
		metaTime=pupaTime;
	}else{
		metaTime=femaleTime;
	}
	if(metaTime < maleTime)
	{
		metaTime=maleTime;
	}
	
return metaTime;
}

void store_density_evolution(char *contenedor, Float2D_MP *RhoVsT, short int Fecha, float *time_map,float *R,float *Temp)
{
FILE *datos;
char archivo[200];

if(RhoVsT!=NULL)
{
	if(RhoVsT->NoEnsambles <1)
	{
		printf("Numero de ensambles 0 para escribir!!");
		return;
	}
	
	int T_max=RhoVsT->i_max;
	int NoEspecies=RhoVsT->j_max;
	float NoEnsambles=(float)RhoVsT->NoEnsambles;
	
	if(Fecha==1)
	{
		time_t now = time(NULL);
		struct tm * timeinfo;
		timeinfo = localtime ( &now );
		
		sprintf(archivo,"%s/density_evolution_%d-%d-%d_%d:%d:%d",contenedor,timeinfo->tm_mday,(timeinfo->tm_mon + 1),(timeinfo->tm_year + 1900),timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
	}else{
		sprintf(archivo,"%s/density_evolution",contenedor);
	}
	
	datos=fopen(archivo, "w");
	fprintf(datos,"t pupae(1) female(2) male(3) total fisical_time R Temp\n");
	
	
	int T,e;
	float Total;
	for(T=0;T<=T_max;T++)
	{
		fprintf(datos,"%d",T);
		Total=0.0;
		for(e=1;e<=NoEspecies;e++)
		{
			fprintf(datos," %f",RhoVsT->array[T][e]/NoEnsambles);
			Total+=RhoVsT->array[T][e]/NoEnsambles;
		}
		fprintf(datos, " %f",Total);
		if(time_map != NULL)
		{
			fprintf(datos, " %f",time_map[T]);
		}
		if(R != NULL)
		{
			fprintf(datos, " %f",R[T]);
		}
		if(Temp != NULL)
		{
			fprintf(datos, " %f",Temp[T]);
		}
		fprintf(datos, "\n");
		//printf("guadando:%d %f %d NoEnsambles=%d \n",T,RhoVsT->array[T][e]/NoEnsambles,e, RhoVsT->NoEnsambles);	
	}
	fclose(datos);
}
return;
}


float calendar_temperature(float FisicalTime)
{
float Temperature;
int nextMonth;
float MonthlyTemp[12];
float UnitsPerMonth = 30.0;
int month = (int)(FisicalTime/UnitsPerMonth);
float day = FisicalTime - month*30.0;

month=(month-(month/12)*12);

MonthlyTemp[0]=13.9;
MonthlyTemp[1]=15.1;
MonthlyTemp[2]=17.1;
MonthlyTemp[3]=19.0;
MonthlyTemp[4]=19.8;
MonthlyTemp[5]=19.4;
MonthlyTemp[6]=18.4;
MonthlyTemp[7]=18.4;
MonthlyTemp[8]=18.2;
MonthlyTemp[9]=17.3;
MonthlyTemp[10]=15.8;
MonthlyTemp[11]=14.5;

nextMonth= (month + 1) - ((month + 1)/12)*12;
Temperature = MonthlyTemp[month] + day*(MonthlyTemp[nextMonth] - MonthlyTemp[month])/30.0;

return Temperature;
}

float calendar_humidity(float FisicalTime)
{
float H;
float MonthlyHumidity[12];
int nextMonth;
float UnitsPerMonth = 30.0;
int month = (int)(FisicalTime/UnitsPerMonth);
float day = FisicalTime - month*30.0;

month=(month-(month/12)*12);

MonthlyHumidity[0]=1.0;
MonthlyHumidity[1]=1.0;
MonthlyHumidity[2]=1.0;
MonthlyHumidity[3]=1.0;
MonthlyHumidity[4]=1.0;
MonthlyHumidity[5]=1.0;
MonthlyHumidity[6]=0.0;
MonthlyHumidity[7]=0.0;
MonthlyHumidity[8]=0.0;
MonthlyHumidity[9]=0.0;
MonthlyHumidity[10]=0.0;
MonthlyHumidity[11]=0.0;

nextMonth= (month + 1) - ((month + 1)/12)*12;

H = MonthlyHumidity[month] + day*(MonthlyHumidity[nextMonth] - MonthlyHumidity[month])/30.0;

return H;
}

void set_param_temperature_dependent(mozquito_parameters *param,float temperature)
{
	param->FemaleDeadRate = feamale_mortality_rate(temperature);
	param->FemaleOffspringRate = feamale_oviposition_rate(temperature);
	param->PupaDeadRate = aquatic_mortality_rate(temperature);
	param->PupaOffspringRate = aquatic_transition_rate(temperature);
	param->Metabolic_Time = obtain_metabolic_time(param);
	return;
}

void set_diapause_humidity_dependent(mozquito_parameters *param,float humidity)
{
	float pupaPartialRate,PupaDeadRate,PupaOffspringRate, excess;
	PupaDeadRate = param->PupaDeadRate;
	PupaOffspringRate = param->PupaOffspringRate;
	
	if(humidity > 1.0)
	{
		humidity=1.0;
	}
	if(humidity < 0.0)
	{
		humidity=0.0;
	}
	param->AntiDiapauseRate = param->Metabolic_Time*humidity;
	param->DiapauseRate  = param->Metabolic_Time*(1.0 - humidity);
	pupaPartialRate = PupaDeadRate + PupaOffspringRate ;
	excess = pupaPartialRate + param->DiapauseRate - param->Metabolic_Time ;	
	if(excess > 0.0)
	{
		param->PupaDeadRate = PupaDeadRate - excess*(PupaDeadRate/pupaPartialRate);
		param->PupaOffspringRate = PupaOffspringRate - excess*(PupaOffspringRate/pupaPartialRate);
	}
return;
}

float feamale_mortality_rate(float T)
{
	float a,a1,a2,a3,a4,R;
	a=0.8692;
	a1=-0.159;
	a2=0.01116;
	a3=-0.0003408;
	a4=0.000003809;
	R=a + a1 * T + a2 * T*T + a3*T*T*T + a4*T*T*T*T;
	return R;
}

float feamale_oviposition_rate(float Temp)
{
	float a,a1,a2,a3,a4,R;
	a=-5.4;
	a1=1.8;
	a2=-0.2124;
	a3=0.01015;
	a4=-0.0001515;
	R=a + a1 * Temp + a2 * Temp*Temp + a3*Temp*Temp*Temp + a4*Temp*Temp*Temp*Temp;
	return R;
}

float aquatic_mortality_rate(float Temp)
{
	float a,a1,a2,a3,a4,R;
	a=2.130;
	a1=-0.3797;
	a2=0.02457;
	a3=-0.0006778;
	a4=0.000006794;
	R=a + a1 * Temp + a2 * Temp*Temp + a3*Temp*Temp*Temp + a4*Temp*Temp*Temp*Temp;
	return R;
}

float aquatic_transition_rate(float Temp)
{
	float a,a1,a2,a3,a4,a5,a6,a7,R;
	a=0.131;
	a1=-0.05723;
	a2=0.01164;
	a3=-0.001341;
	a4=0.00008723;
	a5=-0.000003017;
	a6=0.00000005153;
	a7=-0.000000000342;
	R=a + a1 * Temp + a2 * Temp*Temp + a3*Temp*Temp*Temp + a4*Temp*Temp*Temp*Temp;
	return R;
}
