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


#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#pragma once //specifies that the file will be included only once by compiler (can reduce build times)

#include<vector>

//-------------------------------------------------------------------
double rand01(void);
int randMinMax(int min,int max);

//-------------------------------------------------------

	typedef std::vector<std::vector<double> > array2Ddouble;
	typedef std::vector<std::vector<int> > array2Dint;

	typedef std::vector<double> vecdouble;
	typedef std::vector<int> vecint;

	void SetArraySize2d(std::vector<std::vector<int>> &, int, int);
	void SetArraySize2d(std::vector<std::vector<int>> &, int, int, int);
	void SetArraySize2d(std::vector<std::vector<double>> &, int, int);
	void SetArraySize2d(std::vector<std::vector<double>> &, int, int, double);
	void SetVectorSize(std::vector<int> & , int);
	void SetVectorSize(std::vector<double> & , int);
	
	double mean(array2Dint values,int rows,int index_col);
	double sd(array2Dint values,int rows,int index_col);
	double mean(array2Ddouble values,int rows,int index_col);
	double sd(array2Ddouble values,int rows,int index_col);
	int max_array(array2Dint values,int rows,int index_col);
	int min_array(array2Dint values,int rows,int index_col);
	double max_array(array2Ddouble values,int rows,int index_col);
	double min_array(array2Ddouble values,int rows,int index_col);
	int max_vec(vecint vec);
	double max_vec(vecdouble vec);

//-------------------------------------------------------------------
#endif