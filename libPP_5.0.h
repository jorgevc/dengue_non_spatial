/** Contains a coordinate pair. 
 * Generaly it is used to label a site on a lattice.
*/

#ifndef LIBPP_50
#define LIBPP_50

typedef struct {
int i;			/**< i coordinate, could be thinked as the X coordinate*/ 
int j;			/**< j coordinate, could be thinked as the Y coordinate*/ 
} sitio;

/** Contains the parameters describing the individuals of a species. */
typedef struct {
float Birth;			/**< Birth rate of the species*/
float Coagulation;		/**< Inter-Competition rate of the species. 
							Is the rate at wich the precence of an individual of this specie led a dead 
							of an idividual of OTHER species*/
float Dead;				/**< Intrinsic Dead rate.*/
float CoagulationIntra;	/**< Intra-Competition rate.
							Is the rate at wich the precence of an individual of this specie led a dead 
							of an idividual of the SAME species*/
int RadioBirth;			/**< Radio within a descendant could born in units of lattice sites. */
int RadioCoa;			/**< Radio within Inter-Competion acts in units of lattice sites. */
int RadioCoaIntra;		/**< Radio within Intra-Competion acts in units of lattice sites. */
} especie;

/** Contains the state of the system(lattice). */
typedef struct {
int NDX;			/**< The size of the lattice in the i coordinate (X coordinate). */
int NDY;			/**< The size of the lattice in the j coordinate (Y coordinate). */
int T;				/**< The Time-steeps that the sistem has been evolve since the begining of the simulation. 
						It is incremented by one, each time the instance is passed to BarrMCcRyCamp(estado *es) */
int ON;				/**< The total number of occupied sites in the lattice. */
int **s; 			/**< A 2 dimensional array of int representing the actual lattice(system). 
						A value of 0 in s[i][j] represent an empy site in (i,j). A value > 0 represent an occupy site.*/
int **INDICE;		/**< A 2 dimensional array which value at INDICE[i][j] is the index of an OCCUPIED site at (i,j) on a list 
						of the occupied sites. If (i,j) is not occupyed the return value is undetermined. */
sitio *SO;			/**< A 1 dimensional array of scructs sitio that is a list of the occupied sites. 
						To find the index of an specific occupied site at (i,j) in this list, one have to look for it in INDICE[i][j]. */
int **TIPO;			/**< A 2 dimensional array which value at TIPO[i][j] is used to label the species of the individual in (i,j). 
						If the site in (i,j) is not occupyed the return value is undetermined. */
} estado;

/** General porpuse 2 dimensinal array of float suitable to be used for ensemble averages. */
typedef struct {
float **array;		/**< The actual 2 dimenasional array of float: <Float2D_MP>.array[i][j]. */
int i_max;			/**< The maximun value in the i dimension */
int j_max;			/**< The maximun value in the j dimension */
int T;				/**< Used to store the time-step of the property that is stored in the instance(if aplicable). */
int NoEnsambles;	/**< The number of ensambles that are stored in the instance. This number is generally used to obtain the mean in a
						multi-thread program, or simplie the ensemble average. */
} Float2D_MP;

/** General porpuse 1 dimensinal array of float suitable to be used for ensemble averages. */
typedef struct {
float *array;		/**< The actual 1 dimenasional array of float: <Float1D_MP>.array[i]. */
int i_max;			/**< The maximun value of i. It is equal to (array_size + 1) */
int T;				/**< Used to store the time-step of the property that is stored in the instance(if aplicable). */
int NoEnsambles;	/**< The number of ensambles that are stored in the instance. This number is generally used to obtain the mean in a
						multi-thread program, or simplie the ensemble average. */
} Float1D_MP;

/** General porpuse 1 dimensinal array of float suitable to be used for ensemble averages. */
typedef struct {
int **array;		/**< The actual 2 dimenasional array of float: <Int2D_MP>.array[i][j]. */
int i_max;			/**< The maximun value in the i dimension */
int j_max;			/**< The maximun value in the j dimension */
int NoEnsambles;	/**< The number of ensambles that are stored in the instance. This number is generally used to obtain the mean in a
						multi-thread program, or simplie the ensemble average. */
} Int2D_MP;

/** Struct that represent a distribution on a 1 dimensional real space. i.e. Represents a distribution f(x) where x could be real. */
typedef struct {
int *array;			/**< The values of the distribution. */
int T;				/**< Used to store the time-step of the property that is stored in the instance(if aplicable). */
int i_max;			/**< The maximun value of i. It is equal to (array_size + 1) */
float TamParticion; /**< The maximum resolution that the distribution can handle on domain.  */
int NoEnsambles;	/**< The number of ensambles that are stored in the instance. This number is generally used to obtain the mean in a
						multi-thread program, or simplie the ensemble average. */
float xIni;			/**< The initial value of X where the distribution is "defined". */
float xFin;			/**< The final value of X where the distribution is "defined". */
} Dist_MP;

typedef struct {   /** TIPO = 0 : Todos los tipos, s = 0 : Todos los tamanos, NEG = 0 : No negacion, on : resultado del Numero de elementos del grupo que se asigno la ultima vez que se proceso con CFFT_Univ_MP */
int TIPO;
int s;
int NEG;
int on;
} Grupo;

typedef struct {
int MeanSquare;
int NoEnsambles;
int NoMuestras;
int Muestra;	/** Muestra = 0 : Toma todas las muestras. */
} CorrDescriptor;

extern Grupo GRUPO_INI;

/** Set the value of the Birth rate of a species. 
 * @param[in] L the value of the Birth rate.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the crated species to save memory. 
 * 			After seting all the values for a specie it should be call the function @see EscalaTiempoMetabolico(int tipo).
 * @see especie
 * */
void SetBirth(float L, int tipo);

/** Set the value of the Inter-Coagulation rate of a species. 
 * @param[in] e the value of the Inter-Coagulation rate.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero). 
 * 			After seting all the values for a species it should be call the function @see EscalaTiempoMetabolico(int tipo).
 * @see especie
 * */
void SetCoagulation(float e, int tipo);

/** Set the value of the Intra-Coagulation rate of a species. 
 * @param[in] e the value of the Intra-Coagulation rate.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero). 
 * 			After seting all the values for a species it should be call the function @see EscalaTiempoMetabolico(int tipo).
 * @see especie
 * */
void SetCoagulationIntra(float e,int tipo);

/** Set the value of the Intrinsic Dead rate of a species. 
 * @param[in] d the value of the Intrinsic Dead rate.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero). 
 * 			After seting all the values for a species it should be call the function @see EscalaTiempoMetabolico(int tipo).
 * @see especie
 * */
void SetDead(float d, int tipo);

/** Set the value of the Birth range of a species in number of lattice sites. 
 * @param[in] rb the radio in number of lattice sites of the range of the new offsprings.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero).
 * @see especie 
 * */
void SetRadioBirth(int rb, int tipo);

/** Set the value of the Inter-Coagulation range of a species in number of lattice sites. 
 * @param[in] rc the radio in number of lattice sites of the range of Inter-Coagulation..
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero).
 * @see especie 
 * */
void SetRadioCoa(int rc, int tipo);

/** This shuld be called when a rate of the specie parameters is updated, or after setting up all the rates for the first time.  
 * It sets an internal scale factor of the library that is used to map the rates given in units of phisical time, to rates in units of computational Time-steps.
 * @param tipo the number that identifies the specie. 
 * */
void EscalaTiempoMetabolico(int tipo);

/** allocate the necesary memory for a lattice (system). 
 * @param[out] *es pointer to the lattice (system)
 * @param NDX the size of one side of the lattice: i=X coordinate 
 * @param NDY the size of the other side of the lattice: j=Y coordinate
 * @see estado 
 * */
void AlojaMemoria(estado *es, int NDX, int NDY);

/** Initialize a lattice (system) to be empty. 
 * @param *es pointer to the lattice (system).
 * */
void ResetEstado(estado *es);

/** Writtes to a <estado> a uniform random state of occupied sites of a species. 
 * @param *es lattice (system) where the random state will be "written".
 * @param[in] frac the fraction of the lattice that will be occupied for the generated random state.
 * @param[in] tipo the specie of the generated random state.
 * */
void GeneraEstadoAleatorio(estado *es, float frac, int tipo);

/**
 * Update an individual of a system <*es>. 
 * @param *es the system (lattice)
 * @param N the index of the individual that is going to be updated @see estado::INDICE
 * @param campo In the niche model it is used to simulate the enviromental niche.
 */
void ActualizaRyC(estado *es, int N, int campo);

/**
 * Makes the system <*es> evolve one Time-step. 
 * One Time-step is counted when on average all the individuals in the system has been updated once.
 * @param *es pointer to the struct 'estado' that is going to be evolved. 
 */
void BarrMCcRyCamp(estado *es);

/**
 * Chouse a neighbour(occupied site) of (i,j) within a radio with equal probabiliy of each neighbour to be picked. 
 * @param[in] i,j the site from wich the neighbour is going to be picked up.
 * @param[in] radio the radio within the neighbour is going to be picked up.
 * @param[out] *vecino pointer to a struct vecino where the coordinates of the neighbour are going to be stored.
 */
void EligeUniforme(int i,int j,int radio, sitio *vecino);

/**
 * Inserts individuals of a species in uniform distributed random sites.
 * @param *es System where the individuals are going to be inserted.
 * @param N Number of individuals to be inserted.
 * @param tipo The kind of species to be inserted @see especie , @see SetSpecie2
 */
void InsertaIndividuosAleatorio(estado *es, int N, int tipo);

/**
 * Deprecated(Not efficient): Gives the mean number of neighbours of the system at a given radio.
 * @param *es the system to analize
 * @param radio the radio at wich the mean neighbours are going to be counted.
 */
float OnPromRadio(estado *es, int radio);  //Deprecated

/**
 * Deprecated(Not efficient), Gives the correlation at a radio
 * @param *es System to be analized
 * @param radio distance at wich obtain the correlation.
 */
float FuncionCorrelacion(estado *es,int radio);			//No eficiente

/**
 * Deprecated(Not efficient), Gives the correlation at a radio
 * @param *es System to be analized
 * @param radio distance at wich obtain the correlation.
 */
float FuncionCorrelacion2(estado *es,int radio);		//No eficiente


/**
 * Deprecated(Not efficient), Gives the correlation at a radio but just taking into account 4 points. 
 * @param *es System to be analized
 * @param radio distance at wich obtain the correlation.
 */
float Correlacion(estado *es,int radio);  //No eficiente

/**
 * Deprecated(Not efficient), Gives the correlation at a radio betwen two species
 * @param *es System to be analized
 * @param radio distance at wich obtain the correlation.
 * @param TipoOrigen The origin species from wich the correlation is going to be obtained.
 * @param TipoDistante Te target species to wich the correlation is going to be obtained.
 */
float FuncionCorrelacionSpecies(estado *es,int radio,int TipoOrigen, int TipoDistante); //No eficiente

/**
 * Deprecated(Not efficient), Gives the correlation at a radio betwen two species with just 4 points. 
 * @param *es System to be analized
 * @param radio distance at wich obtain the correlation.
 * @param TipoOrigen The origin species from wich the correlation is going to be obtained.
 * @param TipoDistante Te target species to wich the correlation is going to be obtained.
 */
float CorrelacionEspecies(estado *es,int radio,int TipoOrigen, int TipoDistante);

/**
 * Returns the number of individuals of a single species in a system. 
 * @param *es The system to be analized.
 * @param tipo the kind of species which indiviudals are going to be counted.
 */
int CuentaEspecie(estado *es, int tipo);

/**
 * Inserts and Individual of a specie to a specified location. 
 * @param *es The system in wich the individual is going to be inserted.
 * @param i,j the site where the individual is going to be inserted.
 * @param tipo The kind of specie of the individual to be inserted. @see especie , @see SetSpecie2. 
 */
void InsertaIndividuoEn(estado *es,int i,int j,int tipo);

/**
 * Allocate memory to store the parameters of a specie. Lowest numbers as possible should be used to save memory.
 * @param tipo The species to be allocated. 
 * @see especie
 */
void AlojaMemoriaEspecie(int tipo);

/**
 * Writes the density of each specie into a <Float2D_MP> and the distribution of densities among species in a <Dist_MP>.
 * @param *es System to be analized.
 * @param *RhoVsT The density of each specie is stored in <RhoVsT>[T][i] where T is the Time-steep taken from <es> and i is the species. 
 * This function does not change any other values of <Float2D_MP>, so it can be used each time steep to store the evolution of the densities.
 * This parameter can be NULL, if just the distribution is needed.
 * @param *RhoDist The distribution of densities is stored in <RhoDist> and assigned the Time-steep of <es>. This parameter can be NULL, 
 * if just the densities are needed.
 * @see Float2D_MP
 * @see Dist_MP
 */
void ActualizaRhoVsT_MP(estado *es,Float2D_MP *RhoVsT,Dist_MP *RhoDist);

/**
 * Allocates the necesary memory for a <Float2D_MP>.
 * @param *ARRAY The descriptor <Float2D_MP> of the array that is going to be allocated.
 * @see Float2D_MP
 */
void IniciaMemoriaFloat2D_MP(Float2D_MP *ARRAY);

/**
 * Allocates the necesary memory for a <Int2D_MP>.
 * @param *ARRAY The descriptor <Int2D_MP> of the array that is going to be allocated.
 * @see Int2D_MP
 */
void IniciaMemoriaInt2D_MP(Int2D_MP *ARRAY);

/**
 * Allocates the necesary memory for a <Dist_MP>.
 * @param *Dist The descriptor <Dist_MP> of the distribution that is going to be allocated.
 * @see Dist_MP
 */
void IniciaMemoriaDist_MP(Dist_MP *Dist);

/**
 * Adds in pairs each element of two arrays <Float2D_MP>
 * @param *Origen One array to be added. This array is not changed.
 * @param *Destino Other array to be added. The result is stored in this same array. The final <Destino>.NoEnsambles is the addition of the original two arrays.
 * The final <Destino>.T takes the value of <Origen>.T 
 * @see Float2D_MP
 */
void SumaFloat2D_MP(Float2D_MP *Origen, Float2D_MP *Destino);

/**
 * Adds two <Dist_MP> to make ensemble avearges.
 * @param *Origen One distribution to be added. This distribution is not changed.
 * @param *Destino Other distribution to be added. The result is stored in this distribution. <Destino>.NoEnsambles thakes the value of the addition of the
 * original <Destino>.NoEnsambles plus <Origen>.NoEnsambles
 * @see Dist_MP
 */
void SumaDist_MP(Dist_MP *Origen, Dist_MP *Destino);

/**
 * Initialize and allocate the memory of a <Float2D_MP> object. 
 * Note: This shold be called without previously have been called IniciaMemoriaFloat2D_MP on the same object !
 * The initial values of the array are set to 0.
 * @param *Objeto Pointer to the <Float2D_MP> to be inicialized.
 * @param i_max,j_max Maximun values of i,j in the array <Objeto>.array[i][j]. This is needed to allocate the neccesary Memory.
 * @param NoEnsambles The value of <Objeto>.NoEnsambles that is going to be set.
 * @see Float2D_MP
 */
void InicializaFloat2D_MP(Float2D_MP *Objeto, int i_max, int j_max, int NoEnsambles);

/**
 * Initialize and allocate the memory of a <Dist_MP> object.
 * Note: This shold be called without previously have been called IniciaMemoriaDist_MP on the same object !
 * The initial values of the distribution are set to 0.
 * @param *Objeto Pointer to the <Dist_MP> to be inicialized.
 * @param TamParticion Maximun resolution in the domain that the distribution is going to be able to store. Remember: More resolution translate in more memory needs.
 * @param xIni The closest value to 0 of the domain that the distribution is going to be able to store. Remember: Bigger domain is equal to more memory needs.
 * @param xFin The farest value to 0 of the domain that the distribution is able to store.  Remember: Bigger domain is equal to more memory needs.
 * Example: if you one to store a distribution f(x) where 3.5<x<10.7 with a resolution of 0.1, then TamParticion=0.1, xIni=3.5 and xFin=10.7
 */
void InicializaDist_MP(Dist_MP *Objeto, float TamParticion, float xIni, float xFin);

/**
 * Deprecated: Same as SetSpecie2 but without setting Intra-Competition parameters.
 * @see SetSpecie2
 */
void SetSpecie(int NoEspecie, float Birth, float Coagulation, float Dead, float RadioBirth, float RadioCoa);

/**
 * Resets to 0 the values of a <Float2D_MP> object.
 * The values of the array in are set to 0 , T is set to 0 and NoEnsambles is set to 0.
 * @param *ARRAY pointer to the <Float2D_MP> to be reset.
 * @see Float2D_MP
 */
void ResetFloat2D_MP(Float2D_MP *ARRAY);

/**
 * Resets to 0 the values of a <Dist_MP> object.
 * The values of the array in are set to 0 , T is set to 0 and NoEnsambles is set to 0. TamParticion, xIni, xFin are not changed.
 * @param *Dist pointer to the <Dist_MP> to be reset.
 * @see Dist_MP
 */
void ResetDist_MP(Dist_MP *Dist);

/**
 * Set the values and allocate the necesary memory to define a new specie in to the library.
 * @param NoEspecie The label of the specie to be set.
 * @param Birth The Birth rate.
 * @param Coagulation The Inter-Competion rate. Is the rate at wich the precence of an individual of this specie led a dead of an idividual of OTHER species.
 * @param CoagulationIntra The Intra-Competition rate. Is the rate at wich the precence of an individual of this specie led a dead of an idividual of the SAME species.
 * @param Dead The Intrinsic Dead rate.
 * @param RadioBirth Radio within a descendant could born in units of lattice sites.
 * @param RadioCoa Radio within Inter-Competion acts in units of lattice sites.
 * @param RadioCoaIntra Radio within Intra-Competion acts in units of lattice sites.
 */
void SetSpecie2(int NoEspecie, float Birth, float Coagulation, float CoagulationIntra, float Dead, float RadioBirth, float RadioCoa, float RadioCoaIntra);

/**
 * Deprecated: Employs Correlacion that is Deprecated. 
 */
void ActualizaCorrelacion_MP(estado *es, Float1D_MP *corr);

/**
 * Deprecated: Employs CorrelacionEspecies that is Deprecated. 
 */
void ActualizaCorrelacionTipo_MP(estado *es, Float1D_MP *corr, int TipoOrigen, int TipoObjetivo);

/**
 * Add the arrays of two <Float1D_MP> objects.
 * @param *Origen pointer to one <Float1D_MP> to be added. This object is not changed.
 * @param *Destino pointer to other <Float1D_MP> to be added. The result is stored in <Destino>. The resulting <Destino>.NoEnsambles is the sum of the original two
 * objects. The final <Destino>.T is set to be the same as <Origen>.T 
 */
void SumaFloat1D_MP(Float1D_MP *Origen,Float1D_MP *Destino);

/**
 * Initializes and allocate the memory for a <Float1D_MP> object.
 * The values of the array are set to 0 with calloc (check before critical use), the member T is set to 0 and also the Member NoEnsambles is set to 0.
 * @param *Objeto A pointer to the struct <Float1_MP> to be initialized.
 * @param i_max The maximum value of the array in <Float1_MP>. This is needed for memory allocation.
 */
void InicializaFloat1D_MP(Float1D_MP *Objeto, int i_max);

/**
 * Deprecated: Obains the 2 dimensional correlation and writtes it to a hardcoded file.
 * Was used for testing porpuses only.
 */
void CorrXY(estado *es);

/**
 * Deprecated: to few options plus not treat safe. Obtains the 2 dimensional correlation of all individuals of a system. 
 * @param *es state to be analized.
 * @param *correlacion Where the correlation is stored.
 */
void CFFT(estado *es, Float2D_MP *correlacion);

/**
 * Obtains the 2 dimensional correlation between all the individuals of a system.
 * @param *es array of systems from wich the correlation is going to be obtained. The final correlation is the mean of the correlation of each system.
 * @param NoEnsambles Number of systems to analizes. This must not be greater than the number of systems in the the array <*es> .
 * @param *correlation Here is where the correlation is going to be add. The final result is the addition of what it is inside <correlacion> originaly and the result
 * of the correlation. This is used to obtain the mean of the correlation of many treats. The resulting <correlacion>.NoEnsambles is its original value plus
 * the correlations that the function were able to obtain. CFFT_Tipos_MP is a more powerful function.
 * @see CFFT_Tipos_MP 
 */
void CFFT_MP(estado *es, int NoEnsambles, Float2D_MP *correlacion);

/**
 * Obtains the radial correlation from a 2 dimensional correlation
 * @param *corr2D The 2 dimensional correlation obtained from @see CFFT_Tipos_MP or @see CFFT_MP
 * @param *corrRadial The obtained radial correlation is stored in this object. It inherits the values of the members T and NoEnsambles .
 * This correlation is the halve of the shortest side of <corr2D>, because it is assumed periodic boundaries.
 */
void CompactaCorrelacion(Float2D_MP *corr2D, Float1D_MP *corrRadial);

/**
 * Obtains the 2 dimensional correlation between the individuals of two species of the system.
 */
void CFFT_Tipos_MP(estado *es, int NoEnsambles, Float2D_MP *correlacion,int TipoOrigen, int TipoDestino);

/**
 * Reset the values of a <Float1D_MP> object to 0.
 * The values of the array are set to 0, the value of T is set to 0, the value of NoEnsambles is set to 0.
 * @param *ARRAY the object to be reset.
 */
void ResetFloat1D_MP(Float1D_MP *ARRAY);

/**
 * Obtain the mean of the four cuadrants of a 2 dimensional correlation.
 * This is useful for periodic boundary conditions where the 2 dimensional correlation is periodic in each cuadrant.
 * @param *corr2D the 2 dimensional correlation to be folded. The result is stored in this same object in the cuadrant with lowest indices(left-down cuadrant).
 */
void DoblaCorrelacion(Float2D_MP *corr2D);

/**
 * Makes the integral of an array, with a hard coded integrand. 
 * @param *Function array to be integrated
 * @param inicial lower bound of the integral.
 * @param final upper bound of the integral
 */
float Integra(Float1D_MP *Funcion, int inicial, int final);

/**
 * Free the memory allocated for a system <estado>.
 * @param *es System to be freed.
 */
void LiberaMemoria(estado *es);

/**
 * Free the memory allocated for a <Float2D_MP> object.
 * @param *ARRAY the object to be freed.
 */
void LiberaMemoriaFloat2D_MP(Float2D_MP *ARRAY);

/**
 * Free the memory allocated for a <Float1D_MP> object.
 * @param *Ojbeto the object to be freed.
 */
void LiberaMemoriaFloat1D_MP(Float1D_MP *Objeto);
void CFFT_Univ_MP(estado *es, CorrDescriptor *Especifica, Float2D_MP *correlacion, Grupo *TipoOrigen, Grupo *TipoDestino);

#endif /* !LIBPP_50 */
