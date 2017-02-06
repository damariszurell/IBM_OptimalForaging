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


#include <math.h>

#include "Parameter.h"
#include "functions.h"

//-----------------------------------------------

//set parameters
	float Parameter::cellwidth=0.5;	//cellwidth [km]
	int Parameter::max_days=60;
	int Parameter::spinup=5;

	int Parameter::economy=1;	//1=area/time-minimising; 2=resource/energy-maximising
	int Parameter::forageChoice=1; //1= mix random&optimal patch selection (with prob_randomPatch); 2=lottery patch selection
	int Parameter::DistmaxSetting=2; //1=fixed distmax for all patch selection strategies; 2= emergent distmax for optimal&lottery patch selection
	bool Parameter::GPStag=false; // Are birds GPS tagged?
	bool Parameter::MemoryMap=false; //remember previously visites patches?
	int Parameter::memoryweight=1;	//memory weighting factor

	int Parameter::distmax=10;
	int Parameter::distnearby=5;
	double Parameter::prob_nearby=0.72; 
	double Parameter::prob_randomPatch=0.25;
	int Parameter::max_birds=100000;	/*set default very large! this parameter
										is only for manipulation purposes 
										e.g. single home range simulation*/
	double Parameter::pMale=0.5;

	double Parameter::emax=1330;	//max. intake rate [kilojoules per hour] 1330
	int Parameter::em=1613;	//existence metabolism of adults[kilojoules per day]
	int Parameter::breeding_cost=1613; //extra energy costs during nest buidling, breeding, incubation [kJ per day]
	int Parameter::ef=40;	//flight costs [kilojoules per km]
	int Parameter::hours_prey=18;	//active hours prey
	
	bool Parameter::feed_fullbrood=false;
	bool Parameter::marry_quick=true;
	bool Parameter::pseudocolonizer=true;

	double Parameter::tf=0.111;	//flight time per hour per km
	int Parameter::foragetime_max=18;	//available foraging time per day
	int Parameter::foragetime_brooding=9;	//available foraging time per day
	double Parameter::tau_max_good=2.0;	//max. duration foraging trip [h]
	double Parameter::tau_max_bad=tau_max_good/2;
	double Parameter::tau_max=tau_max_good;
	double Parameter::delta_tau=0.01;	//step size to increase duration of foraging trip [h]
	double Parameter::tol_starve=0.2;	//tolerate food shortage
	double Parameter::tol_juvstarve=0.2;	//juvenile tolerate food shortage

	int Parameter::init_nchicks=4;
	double Parameter::delta_mass20=75.26316;
	float Parameter::delta_mass30=110.0;
	float Parameter::delta_mass40=20.0;
	int Parameter::age_brooded=20;

//-----------------------------------------------
//Constructor
Parameter::Parameter(void)
{
}

//-----------------------------------------------
//Destructor
Parameter::~Parameter(void)
{
}

//---------------------------------------------------------------
//set maximum foraging time depending on patch quality and previous visits to low quality patches
void Parameter::Set_tau_max(double q, int n_previous)
{
	float n;
	double ran, prob;
	
	if (n_previous<1)	
	{
		if (q<0.65) tau_max=tau_max_bad;
		if (q>=0.65) tau_max=tau_max_good;
	}
	else
	{
		if (n_previous<=6) n= (float) n_previous; else n=6.0;
		prob=exp(n)/exp(6.0);
		ran= rand01();
		
		if (ran<=prob) tau_max=tau_max_good; else tau_max=tau_max_bad;
	}
}

//------------------------------------------------------------------------
// set maximum number of birds
void Parameter::Set_max_birds(int n)
{
	max_birds=n;
}

//------------------------------------------------------------------------
// aiming to feed full brood or not?
void Parameter::Set_fullbrood(bool n)
{
	feed_fullbrood=n;
}

//------------------------------------------------------------------------
// is quick wedding allowed?
void Parameter::Set_quickwedding(bool n)
{
	marry_quick=n;
}

//------------------------------------------------------------------------
// does storks pseudo-colonize?
void Parameter::Set_pseudocolonizer(bool n)
{
	pseudocolonizer=n;
}

//------------------------------------------------------------------------
// set probability/proportion of random patch selection
void Parameter::Set_probRandomPatch(double p)
{
	prob_randomPatch=p;
}

//------------------------------------------------------------------------
// anticipated chick weight for home range search
void Parameter::Set_breedingCost(int x)
{
	breeding_cost=x;
}

//--------------------------------------------------------------------------
// set emax
void Parameter::Set_emax(double e)
{
	emax=e;
}

//--------------------------------------------------------------------------
// set number of spinup days
void Parameter::Set_spinup(int d)
{
	spinup=d;
}