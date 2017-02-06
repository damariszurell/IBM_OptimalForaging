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


#include <algorithm>
using namespace std;

#include "Cell.h"
#include "Parameter.h"


//----------------------------------------------------------------------
//Constructor
Cell::Cell(int x_init, int y_init)
{
	x=x_init;
	y=y_init;
	maxResource=0;
	pseudoResource=0;
	resource=0;
	occupants=0;
	pairs=0;
}

//----------------------------------------------------------------------
//Destructor
Cell::~Cell(void)
{
}

//---------------------------------------------------------------------
//
void Cell::DecrResource(double e)
{
	pseudoResource=std::max(0.0,pseudoResource-e);
}