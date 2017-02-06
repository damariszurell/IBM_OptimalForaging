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
#include <cstdlib>
#include<time.h>
#include<vector>
#include<algorithm>

#include "functions.h"


//-------------------------------------------------------------------
//generates pseudorandom number between 0 and 1
double rand01(void)
{
	return(rand()/32767.0);
}

//-------------------------------------------------------------------
//generates pseudorandom number between min and max
int randMinMax(int min, int max)
{
	int modulo=max-min+1;

	if (modulo>0) return(rand()%modulo+min); else return(min);
}

//----------------------------------------------------------------
void SetArraySize2d(std::vector<std::vector<int > > &ar, int xsize, int ysize)
{
    ar.resize(xsize);
    for (int i = 0; i < xsize; ++i)
    {
        ar[i].resize(ysize);
    }

	for (int i=0; i< xsize; i++)
		for (int j=0; j< ysize; j++)
		{
			ar[i][j]=0;
		}
}

//----------------------------------------------------------------
void SetArraySize2d(std::vector<std::vector<int > > &ar, int xsize, int ysize, int init)
{
    ar.resize(xsize);
    for (int i = 0; i < xsize; ++i)
    {
        ar[i].resize(ysize);
    }

	for (int i=0; i< xsize; i++)
		for (int j=0; j< ysize; j++)
		{
			ar[i][j]=init;
		}
}

//----------------------------------------------------------------
void SetArraySize2d(std::vector<std::vector<double > > &ar, int xsize, int ysize)
{
    ar.resize(xsize);
    for (int i = 0; i < xsize; ++i)
    {
        ar[i].resize(ysize);
    }

	for (int i=0; i< xsize; i++)
		for (int j=0; j< ysize; j++)
		{
			ar[i][j]=0.0;
		}
}

//----------------------------------------------------------------
void SetArraySize2d(std::vector<std::vector<double > > &ar, int xsize, int ysize, double init)
{
    ar.resize(xsize);
    for (int i = 0; i < xsize; ++i)
    {
        ar[i].resize(ysize);
    }

	for (int i=0; i< xsize; i++)
		for (int j=0; j< ysize; j++)
		{
			ar[i][j]=init;
		}
}

//----------------------------------------------------------------
void SetVectorSize(std::vector<int > &vec, int size)
{
    vec.resize(size);

	for (int i=0; i<size; i++)
	{
		vec[i]=0;
	}
}

//----------------------------------------------------------------
void SetVectorSize(std::vector<double > &vec, int size)
{
    vec.resize(size);

	for (int i=0; i<size; i++)
	{
		vec[i]=0.0;
	}
}

//-------------------------------------------------------------------
//compute maximum of array
int max_vec(vecint vec)
{
	int maxval=vec[0];

	for (int i=1;i<vec.size();i++) if (vec[i]>maxval) maxval=vec[i];

	return(maxval);
}

//-------------------------------------------------------------------
//compute maximum of array
int max_array(array2Dint values,int rows,int index_col)
{
	int maxval=values[0][index_col];

	for (int i=1;i<rows;i++)
	{
		if (values[i][index_col]>maxval) maxval=values[i][index_col];
	}

	return(maxval);
}

//-------------------------------------------------------------------
//compute maximum of array
double max_array(array2Ddouble values,int rows,int index_col)
{
	double maxval=values[0][index_col];

	for (int i=1;i<rows;i++)
	{
		if (values[i][index_col]>maxval) maxval=values[i][index_col];
	}

	return(maxval);
}

//-------------------------------------------------------------------
//compute minimum of array
int min_array(array2Dint values,int rows,int index_col)
{
	int minval=values[0][index_col];

	for (int i=1;i<rows;i++)
	{
		if (values[i][index_col]<minval) minval=values[i][index_col];
	}

	return(minval);
}

//-------------------------------------------------------------------
//compute minimum of array
double min_array(array2Ddouble values,int rows,int index_col)
{
	double minval=values[0][index_col];

	for (int i=1;i<rows;i++)
	{
		if (values[i][index_col]<minval) minval=values[i][index_col];
	}

	return(minval);
}

//-------------------------------------------------------------------
//calculate arithmetic mean
double mean(array2Dint values, int rows,int index_col)
{
	int sum=0;
	
	for (int i=0;i<rows;i++)
	{
		sum+=values[i][index_col];
	}
	
	return (double) sum/rows;
}

//--------------------------------------------------------------------
//calculate standard deviation
double sd(array2Dint values,int rows, int index_col)
{
	double av=mean(values,rows,index_col);
	double temp=0;

	for (int i=0;i<rows;i++)
	{
		temp += (av-values[i][index_col])*(av-values[i][index_col]);
	}

	return temp/rows;
}

//-------------------------------------------------------------------
//calculate arithmetic mean
double mean(array2Ddouble values, int rows,int index_col)
{
	double sum=0;
	
	for (int i=0;i<rows;i++)
	{
		sum+=values[i][index_col];
	}
	
	return (double) sum/rows;
}

//--------------------------------------------------------------------
//calculate standard deviation
double sd(array2Ddouble values,int rows, int index_col)
{
	double av=mean(values,rows,index_col);
	double temp=0;

	for (int i=0;i<rows;i++)
	{
		temp += (av-values[i][index_col])*(av-values[i][index_col]);
	}

	return temp/rows;
}



