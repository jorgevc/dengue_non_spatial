#include "libPP_5.0.h"

typedef struct MOZQUITO_PARAMETERS_STRUCT mozquito_parameters;
typedef struct MOZQUITOS_STATE_STUCT mozquitos_state;

struct MOZQUITO_PARAMETERS_STRUCT {
float PupaDeadRate;			
float PupaOffspringRate;
float FemaleOffspringFraction;
float FemaleDeadRate;
float MaleDeadRate;
float FemaleOffspringRate;
float Metabolic_Time;
float DiapauseRate;
float AntiDiapauseRate;		
};

struct MOZQUITOS_STATE_STUCT {
int female;
int male;
};

void MC_sweep_mozquito(estado *es, mozquitos_state *mozquitos, mozquito_parameters *param);

void Update_mozquito(estado *es,mozquitos_state *mozquitos, mozquito_parameters *parameters, int N);

float obtain_metabolic_time(mozquito_parameters *param);

void store_density_evolution(char *contenedor, Float2D_MP *RhoVsT, short int Fecha, float *time_map, float *R, float *Temp);

float calendar_temperature(float FisicalTime);

float calendar_humidity(float FisicalTime);

void set_param_temperature_dependent(mozquito_parameters *param,float temperature);

void set_diapause_humidity_dependent(mozquito_parameters *param,float humidity);

float feamale_mortality_rate(float T);

float feamale_oviposition_rate(float Temp);

float aquatic_mortality_rate(float Temp);

float aquatic_transition_rate(float Temp);
