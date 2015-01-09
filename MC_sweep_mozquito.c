/*
Copyright 2012 Jorge Velazquez
*/
#include "MC_sweep_mozquito.h"
#include "GNA.h"
#include <stdio.h>
#include <time.h>

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
float pPupaDead, pPupaOffspring, pMaleDead, pFemaleDead, pMozquitoOffspring;
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

		if(Rand<=(pPupaOffspring + pPupaDead )) //dead or offspring or nothing
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
	float metaTime,pupaTime,femaleTime,maleTime;
	
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

void store_density_evolution(char *contenedor, Float2D_MP *RhoVsT, short int Fecha)
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
	fprintf(datos,"t pupae(1) female(2) male(3) total \n");
	
	
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
		fprintf(datos, " %f\n",Total);		
		//printf("guadando:%d %f %d NoEnsambles=%d \n",T,RhoVsT->array[T][e]/NoEnsambles,e, RhoVsT->NoEnsambles);	
	}
	fclose(datos);
}
return;
}
