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


#pragma once //specifies that the file will be included only once by compiler (can reduce build times)

#include <math.h>
#include <vector>

#include "Parameter.h"
#include "functions.h"

//----------------------------------------------------------------------------------------------------
class GPS : public Parameter
{
	
public:
	GPS(void);	//Constructor
	~GPS(void);	//Destructor
	
	//functions to access private variables
	std::vector<int> GetGPSday(void) const {return day;}
	std::vector<int> GetGPSx(void) const {return GPS_x;}
	std::vector<int> GetGPSy(void) const {return GPS_y;}
	std::vector<double> GetGPStime(void) const {return GPS_time;}
	std::vector<double> Getflighttime(void) const {return flight_time;}
	std::vector<double> GetAvailableEnergy(void) const {return AvailableEnergy;}
	std::vector<double> GetConsumedEnergy(void) const {return ConsumedEnergy;}
	int GetGPSsize(void) const{return day.size();}
	void AddGPSday(int x) {day.push_back(x);}
	void AddGPSx(int x) {GPS_x.push_back(x);}
	void AddGPSy(int y) {GPS_y.push_back(y);}
	void AddGPStime(double time) {GPS_time.push_back(time);}
	void Addflighttime(double time) {flight_time.push_back(time);}
	void AddAvailableEnergy(double en) {AvailableEnergy.push_back(en);}
	void AddConsumedEnergy(double en) {ConsumedEnergy.push_back(en);}

	//------------------------------------------------------------
private:	
	vecint day;
	vecint GPS_x;
	vecint GPS_y;
	vecdouble GPS_time;
	vecdouble flight_time;
	vecdouble AvailableEnergy;
	vecdouble ConsumedEnergy;
};

