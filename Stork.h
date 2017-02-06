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

#include "Parameter.h"
#include "GPS.h"
#include "functions.h"


//----------------------------------------------------------------------------------------------------
class Stork : public Parameter
{
	
public:
	Stork(int sexed);	//Constructor
	~Stork(void);	//Destructor
	
	//functions to access private variables
	static int GetPopsize(void);

	int GetId(void) const {return identity;}
	int GetSex(void) const {return sex;}
	int GetX(void) const {return x;}
	int GetY(void) const {return y;}
	int GetChicks(void) const {return n_chicks;}
	int GetAgeChicks(void) const {return age_chicks;}
	bool IsPaired(void) const {return paired;}
	bool IsAlive(void) const {return alive;}
	int GetMateID(void) const {return MateID;}
	double GetEGained(void) const {return egained;}
	double GetERequired(void) const {return erequired;}
	double GetEConsumed(void) const {return econsumed;}
	double GetForagingtime(void) const {return foraging_time;}
	double GetForagingtimeperday(void) const {if ((status==2)&&(age_chicks<age_brooded)) return foragetime_brooding; else return foragetime_max;}
	int GetNLowquality(void) const {return n_lowquality;}
	int GetStatus(void) const{return status;}

	void SetXY(int xnew,int ynew) {x=xnew;y=ynew;}	//set location
	void SetStatus(int s);
	void InitERequired(void);
	void IncrERequired(double e); 
	void InitEConsumed(void);
	void InitEGained(void) {egained=0;}
	void IncrEGained(double e)  {egained+=e;}
	void SyncEs(double er, double ec, double eg) {erequired=er; econsumed=ec;egained=eg;}
	void InitChicks(void) {n_chicks=init_nchicks; age_chicks=0; mass_chicks=70;enest=5.6*pow(70,0.81);foraging_time=9;status=breeding;}
	void ChicksAge(void);
	bool StarveChick(void);
	void InitForagingTime(void) {foraging_time=18;}
	void DecreaseForagingTime(double incr); 
	void InitNLowquality(void) {n_lowquality=0;}
	void IncrNLowquality(void) {n_lowquality++;}
	void Pair(int id) {paired=true; MateID=id;}
	void Unpair(void) {paired=false;MateID=0;}
	void Die(void) {alive=false;}

	//------------------------------------------------------------
	static int popsize;	//total number of individuals in simulation

	GPS GPSdata;
	array2Dint memorymap;

private:	
	unsigned int identity;
	enum SEX {Male=1, Female=2} sex; 
	int x;
	int y;
	int n_chicks;
	enum STATUS {searching=0,pre_breeding=1,breeding=2,non_breeding=3} status;
	int age_chicks;
	double mass_chicks;
	double enest;
	double erequired;
	double econsumed;
	double egained;
	double foraging_time;
	int n_lowquality;
	bool paired;
	bool alive;
	int MateID;
};

