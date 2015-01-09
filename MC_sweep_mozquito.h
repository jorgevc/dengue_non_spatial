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
};

struct MOZQUITOS_STATE_STUCT {
int female;
int male;
};

void MC_sweep_mozquito(estado *es, mozquitos_state *mozquitos, mozquito_parameters *param);

void Update_mozquito(estado *es,mozquitos_state *mozquitos, mozquito_parameters *parameters, int N);

float obtain_metabolic_time(mozquito_parameters *param);

void store_density_evolution(char *contenedor, Float2D_MP *RhoVsT, short int Fecha);
