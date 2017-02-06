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
#include <algorithm>
using namespace std;

#include "Stork.h"
#include "Parameter.h"
#include "functions.h"

int Stork::popsize=0;	//initialise number of individuals

//------------------------------------------------------------
//Constructor
Stork::Stork(int sexed)
{
	popsize++;
	identity=popsize;

	if (sexed==1) sex=Male;
	if (sexed==2) sex=Female;
	status=searching;
	x=-99;
	y=-99;
	foraging_time=foragetime_max;
	egained=0;
	n_lowquality=0;
	n_chicks=0;
	age_chicks=-1;
	mass_chicks=0.0;
	enest=0.0;
	paired=false;
	alive=true;
	MateID=0;

	if (MemoryMap) SetArraySize2d(memorymap,max_x,max_y,1.0);

	InitERequired();
	InitEConsumed();
}

//------------------------------------------------------------
//Destructor
Stork::~Stork(void)
{
	//popsize--;
}

//------------------------------------------------------------
//return number of initiated individuals
int Stork::GetPopsize(void)
{
	return popsize;
}

//------------------------------------------------------------
// Set breeding status
void Stork::SetStatus(int s)
{
	if (s==0) status=searching;
	if (s==1) status=pre_breeding;
	if (s==2) status=breeding;
	if (s==3) status=non_breeding;
}

//------------------------------------------------------------
//chicks age one day, mass + energy requirements + foraging time updated
void Stork::ChicksAge(void)
{
	if (n_chicks>0)
	{
		age_chicks++;

		if (age_chicks<20) mass_chicks+=delta_mass20;
		if ((age_chicks>=20)&&(age_chicks<30)) mass_chicks+=delta_mass30;
		if ((age_chicks>=30)&&(age_chicks<40)) mass_chicks+=delta_mass40;

		enest=5.6*pow(mass_chicks,0.81);

		if (age_chicks<age_brooded) foraging_time=foragetime_brooding; else foraging_time=foragetime_max;
	}
}

//--------------------------------------------------------------
//starve chick
bool Stork::StarveChick(void)
{
	if (n_chicks>0) 
	{
		n_chicks--;
		if (n_chicks==0) age_chicks=-1;
		return(true);
	} 
	else return(false);
}

//--------------------------------------------------------------
//set energy requirements of the day - if paired for chicks and both adults 
void Stork::InitERequired(void)
{
	if (status==searching)
	{
		if (marry_quick==true) erequired=2*(breeding_cost+em); else erequired=breeding_cost+em;
	}
	else if (status==pre_breeding)	//note: only pairs enter pre-breeding phase!
	{
		erequired=2*(breeding_cost+em); 
	}
	else if (status==breeding)
	{
		if (feed_fullbrood==true) erequired=init_nchicks*enest+2*em;
		else if (feed_fullbrood==false) erequired=n_chicks*enest+2*em;
	}
	else if (status==non_breeding)
	{
		erequired=em;
	}
}

//--------------------------------------------------------------
//set truely consumed energy of day - if paired for chicks and both adults
void Stork::InitEConsumed(void)
{
	if (status==searching)
	{
		if (marry_quick==true) econsumed=2*(breeding_cost+em); else econsumed=breeding_cost+em;
	}
	else if (status==pre_breeding)	//note: only pairs enter pre-breeding phase!
	{
		econsumed=2*(breeding_cost+em); 
	}
	else if (status==breeding)
	{
		econsumed=n_chicks*enest+2*em;
	}
	else if (status==non_breeding)
	{
		econsumed=em;
	}
}

//---------------------------------------------------------------
// increase required energy in case flight to any patch was of no avail
void Stork::IncrERequired(double e)
{
	erequired+=e;
	econsumed+=e;
}

//------------------------------------------------------------------------
//decrease foraging time, but not less than zero
void Stork::DecreaseForagingTime(double incr)
{
	foraging_time=std::max(0.0,(foraging_time-incr));
}