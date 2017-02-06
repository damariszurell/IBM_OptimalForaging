/* * Individual-based foraging model of breeding white storks
 * -----------------------------
 * Dr. Damaris Zurell
 * Plant Ecology and Nature Conservation
 * University of Potsdam, Germany
 *
 *
 * The model simulates the daily foraging trips of breeding white storks.
 * It simulates the spatial structure and breeding success of white stork populations 
 * in heterogeneous landscapes by explicitly simulating foraging behaviour and 
 * home range formation. Because resource depletion is modelled explicitly, the model can 
 * predict the maximum carrying capacity of white stork breeding populations in different 
 * landscapes, and density-dependent breeding success by inducing fixed stork density 
 * levels below carrying capacity. 
 *
 * Important files:
 *  * breeding.cpp              : main program
 *  * simulation.h              : class definition
 *  * simulation.cpp            : class code
 *  * Cell.h					: class definition
 *  * Cell.cpp					: class code
 *  * Stork.h					: class definition
 *  * Stork.cpp					: class code 
 *  * Colony.h					: class definition
 *  * Colony.cpp				: class code
 *	* GPS.h						: class definition
 *  * GPS.cpp					: class code
 *  * functions.h				: class definition
 *  * functions.cpp				: class code
 *  * Parameter.h				: class definition
 *  * Parameter.cpp				: class code
 *  * INPUT->lc_33_stat         : example landscape
 *
 * Copyright (C) 2014  Dr. Damaris Zurell
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * See the GNU General Public License under <http://www.gnu.org/licenses/>.
 */


#pragma once	//specifies that the file will be included only once by compiler (can reduce build times)

#define max_x 33	
#define max_y 33

//----------------------------------------------------------------------
class Parameter
{
public:
	Parameter(void);
	~Parameter(void);

	// parameters
	static float cellwidth;	//width of grid cell
	static int max_days;	//number simulation days
	static int spinup;	//number of spinup days

	static int economy; // 1=area/time-minimising; 2 = resource/energy-maximising
	static int forageChoice; //1= mix random&optimal patch selection (with prrob_randomPatch); 2=lottery patch selection
	static int DistmaxSetting; //1=fixed distmax for all patch selection strategies; 2= emergent distmax for optimal&lottery patch selection
	static bool GPStag; //are birds GPS tagged?
	static bool MemoryMap; //do storks remember/prefer previously visited patches?
	static int memoryweight; //memory weighting factor 

	static int distmax;	//maximum foraging distance [in cells]
	static int distnearby;	//maximum distance of nearby patches [in cells]
	static double prob_nearby; //probability of nearby random patch selection
	static double prob_randomPatch;	//probability to randomly select foraging patch
	static int max_birds;	//maximum number of birds

	static double pMale;

	static double emax;	//max. energy intake
	static int em;	//existence metabolism of adults
	static int breeding_cost; 
	static int ef;	//flight costs
	static int hours_prey;	//how many hours per day is prey active?

	static bool feed_fullbrood;	//always aim to feed six young or not?
	static bool marry_quick;	
	static bool pseudocolonizer;

	static double tf;	//flight time per hour per kilometer
	static int foragetime_max;	//available foraging time per day
	static int foragetime_brooding; 
	static double tau_max;	//maximum duration of foraging trip
	static double tau_max_good;	//maximum duration of foraging trip in good patches
	static double tau_max_bad;		//maximum duration of foraging trip in bad patches
	static double delta_tau;	//step size to increase duration foraging trip when energy requirements not fulfilled
	static double tol_starve;	//tolerate food shortage before starvation [%]
	static double tol_juvstarve;	//juvenile tolerate food shortage before starvation [%]

	static int init_nchicks; // initial number of juveniles
	static double delta_mass20;	//daily increase in chick body mass until day 20
	static float delta_mass30;
	static float delta_mass40;
	static int age_brooded;	//age until which chicks need brooding

	//functions to manipulate parameters
	static void Set_tau_max(double q, int n_previous);
	static void Set_max_birds(int n);
	static void Set_fullbrood(bool n);
	static void Set_quickwedding(bool n);
	static void Set_pseudocolonizer(bool n);
	static void Set_probRandomPatch(double p);
	static void Set_breedingCost(int x);
	static void Set_emax(double e);
	static void Set_spinup(int d);
};

