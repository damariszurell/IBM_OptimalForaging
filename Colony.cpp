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


#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;
#include <math.h>
#include <list>
#include <vector>

#include "Colony.h"
#include "Parameter.h"
#include "Cell.h"
#include "Stork.h"
#include "functions.h"

#define PI 3.14159265358979323846264338328

//------------------------------------------------------------------------
//Constructor
Colony::Colony(void)
{
	Cell *pCell=0;	//sets pointer to class Cell

	//initiate grid
	for (int x=0;x<max_x;x++)
		for (int y=0;y<max_y;y++)
		{
			pCell = new Cell(x,y);
			Grid[x][y] = pCell;
		}
	pCell=0;

	//ReadLandscape("e:/Projects/DIP_WhiteStork/HomeRangeSim/Code/CiconiaHomeRange/INPUT/lc_33_stat.txt");
	ReadLandscape("INPUT/lc_33_stat_rescaled.txt");
}

//-------------------------------------------------------------------------
//Destructor
Colony::~Colony(void)
{
	Birdlist.clear();
}

//--------------------------------------------------------------------------
//read landscape file and multiply by emax (maximum food availability)
void Colony::ReadLandscape(std::string Landsc_file)
{
	std::ifstream file (Landsc_file.c_str());
	float hsi;

	if (file.is_open())
	{
		//loop over all grid cells
		for (int x=0;x<max_x;x++)
			for (int y=0;y<max_y;y++)
			{
				//read in relative habitat suitability and multiply by maximum resource availability emax
				file >> hsi;
				Grid[x][y]->SetMaxResource(hsi/100*emax);
				Grid[x][y]->InitResource();
			}

		file.close();
	}
	else cout <<"Unable to open " << Landsc_file <<endl;
}

//-------------------------------------------------------------------------------
//add storks until model world is saturated
void Colony::InitBirds(void)
{
	int sex=2, xopt=-99,yopt=-99,count=0;

	bool homerange_found=false, grid_saturated=false, colonizer=false,quickwedding=false,fullfed=false;
	double enorm,eopt,eopt_pertime,er_opt, ft_opt;

	//initialise static member variable popsize (only important for multiple simulations to reset IDs starting from 1)
	Stork dummy(1);
	dummy.popsize=0;

	// start home range search
	while ((grid_saturated==false)&&(count<max_birds))
	{
	eopt=0.0;eopt_pertime=0.0;er_opt=100000.0;ft_opt=0.0;
	colonizer=false;quickwedding=false;fullfed=false;

	if (rand01()<pMale) sex=1; else sex=2;
	//sex=sex%2+1;	//for equal number of males and females entering simulation in turn
	Stork ind(sex);

	Birdlist.push_back(ind);
	list<Stork>::iterator firstbird=Birdlist.begin();
	list<Stork>::iterator lastbird=Birdlist.end();	//refers to element one past last element
	lastbird--;	//decrement iterator to refer to last element of list

	//-------------------------

	count= (int) Birdlist.size();
	enorm= (*lastbird).GetERequired();

	//------------------------------------------
	//test all suitable grid cells as potential central home range cell
	//start with cells that are already occupied (pseudo-colonizers!)
	if (count>1)
	{
	colonizer=true;
	
	//quick marriage?
	if (marry_quick==true)
	{
		list<Stork>::iterator potpartner;
		for (potpartner = Birdlist.begin(); potpartner != lastbird; ++potpartner)
		{
			if ((*lastbird).IsPaired()==false)
			{
				if (((*potpartner).IsPaired()==false) && ((*potpartner).GetSex()!=(*lastbird).GetSex()))
					{
					(*lastbird).Pair((*potpartner).GetId());
					(*potpartner).Pair((*lastbird).GetId());
					
					//set central cell of optimal home range
					int xmarried=(*potpartner).GetX();
					int ymarried=(*potpartner).GetY();
					(*lastbird).SetXY(xmarried,ymarried);
					Grid[xmarried][ymarried]->IncrOccupants();
				}
			} else break;
		}

		if ((*lastbird).IsPaired()==true) {quickwedding=true;homerange_found=true;}
	}

	if (quickwedding==false)
	{
	for (int x=0;x<max_x;x++)
		for (int y=0;y<max_y;y++)
		{
			bool search=false;
			bool better=false;

			if (pseudocolonizer==true) 
			{if ((Grid[x][y]->GetOccupants()>0)&&(IsPatchAvailable(lastbird,x,y)==true)) search=true;}
			else if (pseudocolonizer==false)
			{if (IsPatchAvailable(lastbird,x,y)==true) search=true;}

			if (search==true)
			{
				while (((*lastbird).GetForagingtime()>0)&&((*lastbird).GetEGained()<(*lastbird).GetERequired())&&(IsPatchAvailable(lastbird,x,y)==true))
				{
					if (forageChoice==1) {if (rand01()<prob_randomPatch) RandomPatch(lastbird,x,y); else OptimalPatch(lastbird,x,y);}
					else if (forageChoice==2) LotteryPatch(lastbird,x,y);
				}

				if ((*lastbird).GetForagingtime()<foragetime_max)
				{
				if ((economy==1)&&(((*lastbird).GetEGained()-(*lastbird).GetERequired())>=0) && 
					((*lastbird).GetEGained()/(foragetime_max-(*lastbird).GetForagingtime())>eopt_pertime)) better=true; //area-minimising
				else if ((economy==2)&&(((*lastbird).GetEGained()-(*lastbird).GetERequired())>eopt) && 
					(((*lastbird).GetERequired()-enorm)<=(er_opt-enorm))) better=true; //resource-maximising; the er_opt/enorm accounts for flight costs = increasing energy required
				else better=false;

				if (better==true)	
				{ 
					eopt_pertime=(*lastbird).GetEGained()/(foragetime_max-(*lastbird).GetForagingtime());
					eopt=(*lastbird).GetEGained()-(*lastbird).GetERequired(); 
					xopt=x; yopt=y;
					er_opt=(*lastbird).GetERequired();
					ft_opt=(*lastbird).GetForagingtime();
					fullfed=true;

					for (int xi=0;xi<max_x;xi++)
						for (int yi=0;yi<max_y;yi++)
						{
							Grid[xi][yi]->RememberResource();
						}
				}}

				//restore resources
				for (int xi=0;xi<max_x;xi++)
					for (int yi=0;yi<max_y;yi++)
					{
						Grid[xi][yi]->RestoreResource();
					}
				
				//reset stork's energy levels etc.
				(*lastbird).InitNLowquality();
				(*lastbird).InitForagingTime();
				(*lastbird).InitERequired();
				(*lastbird).InitEConsumed();
				(*lastbird).InitEGained();
			}
		}	
	}
	}	//end loop over occupied cells
	
	//--------------------------------		
	//resources sufficient to cover energy requirements?
	if ((fullfed==false)&&(quickwedding==false))
	{
	//optimal home range found? if not, now search all unoccupied cells (if pseudocolonizer==true ad for first bird)
	homerange_found=false;
	colonizer=false;

	if ((count==1)||(pseudocolonizer==true))
	{
	for (int x=0;x<max_x;x++)
		for (int y=0;y<max_y;y++)
		{
			bool better=false;

			if ((Grid[x][y]->GetOccupants()==0)&&(IsPatchAvailable(lastbird,x,y)==true))
			{

				while (((*lastbird).GetForagingtime()>0)&&((*lastbird).GetEGained()<(*lastbird).GetERequired())&&(IsPatchAvailable(lastbird,x,y)==true))
				{
					if (forageChoice==1) {if (rand01()<prob_randomPatch) RandomPatch(lastbird,x,y); else OptimalPatch(lastbird,x,y);}
					else if (forageChoice==2) LotteryPatch(lastbird,x,y);
				}

				
				if ((*lastbird).GetForagingtime()<foragetime_max)
				{
				if ((economy==1)&&(((*lastbird).GetEGained()-(*lastbird).GetERequired())>=0) && 
					((*lastbird).GetEGained()/(foragetime_max-(*lastbird).GetForagingtime())>eopt_pertime)) better=true; //area-minimising
				else if ((economy==2)&&(((*lastbird).GetEGained()-(*lastbird).GetERequired())>eopt) && 
					(((*lastbird).GetERequired()-enorm)<=(er_opt-enorm))) better=true; //resource-maximising;
				else better=false;

				if (better==true)	
				{ 
					eopt_pertime=(*lastbird).GetEGained()/(foragetime_max-(*lastbird).GetForagingtime());
					eopt=(*lastbird).GetEGained()-(*lastbird).GetERequired(); 
					xopt=x; yopt=y;
					er_opt=(*lastbird).GetERequired();
					ft_opt=(*lastbird).GetForagingtime();
					fullfed=true;

					for (int xi=0;xi<max_x;xi++)
						for (int yi=0;yi<max_y;yi++)
						{
							Grid[xi][yi]->RememberResource();
						}
				}}

				//restore resources
				for (int xi=0;xi<max_x;xi++)
					for (int yi=0;yi<max_y;yi++)
					{
						Grid[xi][yi]->RestoreResource();
					}
				
				//reset stork's energy levels etc.
				(*lastbird).InitNLowquality();
				(*lastbird).InitForagingTime();
				(*lastbird).InitERequired();
				(*lastbird).InitEConsumed();
				(*lastbird).InitEGained();
				//}
			}
		}
	}
	}	//end loop over unoccupied cells


	//---------------------------------------
	//resources sufficient to cover energy requirements?
	if ((fullfed==true)&&(quickwedding==false))
	{
		homerange_found=true;
		//set central cell of optimal home range
		(*lastbird).SetXY(xopt,yopt);
		Grid[xopt][yopt]->IncrOccupants();

		//now truely decrease resources
		for (int xi=0;xi<max_x;xi++)
				for (int yi=0;yi<max_y;yi++)
				{
					Grid[xi][yi]->DecrTrueResource();
					Grid[xi][yi]->RestoreResource();
				}

		//----------------------------------
		//if settled in occupied cell search for partner 
		if (colonizer==true)
		{
		list<Stork>::iterator potpartner;
		for (potpartner = Birdlist.begin(); potpartner != lastbird; ++potpartner)
		{
			if ((*lastbird).IsPaired()==false)
			{
				if (((*potpartner).GetX() ==(*lastbird).GetX()) && ((*potpartner).GetY()==(*lastbird).GetY()) &&
					((*potpartner).IsPaired()==false) && ((*potpartner).GetSex()!=(*lastbird).GetSex()))
					{
					(*lastbird).Pair((*potpartner).GetId());
					(*potpartner).Pair((*lastbird).GetId());
				}
			} else break;
		}
		}
	} 

		//home range found? if not model world saturated
	if (homerange_found==false) 
	{
		Birdlist.erase(lastbird);
		grid_saturated=true;
	}
	//----------------------------------------
	//count number of birds
	count= (int) Birdlist.size();

	} //model world saturated

	//-------------------------------------------
	//loop over all cells to reset resources to max. 
	for (int x=0;x<max_x;x++)
		for (int y=0;y<max_y;y++)
		{
			Grid[x][y]->InitResource();
		}

	//--------------------------------------------
	//loop over all storks to set status and energy requirements etc. 
	list<Stork>::iterator ind;

	for(ind=Birdlist.begin(); ind != Birdlist.end(); ++ind)
	{
		if ((*ind).IsPaired()==true) (*ind).SetStatus(1); else (*ind).SetStatus(3);
		(*ind).InitForagingTime();
		(*ind).InitNLowquality();
		(*ind).InitERequired();
		(*ind).InitEConsumed();
		(*ind).InitEGained();
	}
}

//-------------------------------------------------------------------------------
// initiate breeding
void Colony::InitBreeding(void)
{
	//loop over all storks to find paired ones and start brood
	list<Stork>::iterator ind;

	for(ind=Birdlist.begin(); ind != Birdlist.end(); ++ind)
	{
		if ((*ind).IsPaired()==true) (*ind).InitChicks();
		(*ind).InitERequired();
		(*ind).InitEConsumed();
	}
}

//-------------------------------------------------------------------------------
// at end of day check for starved individuals and chicks, initialise variables
void Colony::GoodNight(void)
{
	bool starve;

	list<Stork>::iterator ind;
	vector<list<Stork>::iterator> deleteList; 
	
	//loop over all birds to check for starving individuals and starved chicks
	for(ind=Birdlist.begin(); ind != Birdlist.end(); ++ind)
	{
		if ((*ind).IsAlive()==true)
		{
		//check for starved individuals and chicks
		if ((*ind).GetEGained()<(*ind).GetEConsumed())
		{
			//paired?
			if ((*ind).IsPaired()==true)
			{
				//find partner
				list<Stork>::iterator partner;
				for (partner = Birdlist.begin(); partner != Birdlist.end(); ++partner)
				{
					if ((*partner).GetId()==(*ind).GetMateID()) break;
				}

				//----------------
				//any chicks?
				if ((*ind).GetChicks()>0)
				{
					starve=true;
					//starve chicks until energy balance ok or until no chicks left
					while (((*ind).GetEGained()<((*ind).GetEConsumed()*(1-tol_juvstarve)))&&(starve==true))
					{
						starve=(*ind).StarveChick();
						(*partner).StarveChick();

						//reset energy requirements
						(*ind).InitEConsumed();
						(*partner).InitEConsumed();
					}
				}

				//-----------------
				//still not enough energy gained? then starve partner (note! adult tolerance against food shortage)
				if ((*ind).GetEGained()<((*ind).GetEConsumed()*(1-tol_starve)))
				{
					int x=(*partner).GetX();
					int y=(*partner).GetY();
					deleteList.push_back(partner);
					(*partner).Die();
					Grid[x][y]->DecrOccupants();
					(*ind).Unpair();
					(*ind).InitEConsumed();
				}

			}	//end statement over pairs and chicks
			
			//-----------------------------------
			// not enough energy for single ind.? (note! adult tolerance against food shortage)
			if ((*ind).GetEGained()<((*ind).GetEConsumed()*(1-tol_starve))) 
			{
				int x=(*ind).GetX();
				int y=(*ind).GetY();
				deleteList.push_back(ind);
				(*ind).Die();
				Grid[x][y]->DecrOccupants();
			}
		}}
	}	//end starvation loop

	//------------------------------------------
	//loop over deleteList to remove dead birds
	for (vector<list<Stork>::iterator>::size_type ix=0;ix!=deleteList.size();ix++)
	{
		Birdlist.erase(deleteList[ix]);
	}

	//------------------------------------------
	//loop over all remaining birds to reset variables
	for(ind=Birdlist.begin(); ind != Birdlist.end(); ++ind)
	{
		(*ind).InitForagingTime();
		(*ind).ChicksAge();	//if n_chicks>0 increment chicks' age, body mass + energy requirements, set foraging time
		(*ind).InitNLowquality();
		(*ind).InitERequired();
		(*ind).InitEConsumed();
		(*ind).InitEGained();
	}

	//-------------------------------------------
	//loop over all cells to reset resources to max. 
	for (int x=0;x<max_x;x++)
		for (int y=0;y<max_y;y++)
		{
			Grid[x][y]->InitResource();
		}

	//-----------------------------
		deleteList.clear();
}

//--------------------------------------------------------------------------------
// foraging over available foraging time
void Colony::Forage(int day)
{
	list<Stork>::iterator ind;
	vector<list<Stork>::iterator> shuffledBirdlist;

	//fill shuffledBirdlist with iterators of Birdlist
	for(ind=Birdlist.begin(); ind != Birdlist.end(); ++ind)
	{
		shuffledBirdlist.push_back(ind);
	}

	// shuffle entries in shuffledBirdlist
	for (int i=0;i<((int) Birdlist.size()-1);i++)
	{
		int r = i + (rand() % ((int) Birdlist.size()-i));
		list<Stork>::iterator temp = shuffledBirdlist[i];
		shuffledBirdlist[i]=shuffledBirdlist[r];
		shuffledBirdlist[r]=temp;
	}

	//loop over all storks in shuffledBirdlist
	for (vector<list<Stork>::iterator>::size_type ix=0;ix!=shuffledBirdlist.size();ix++)
	{
	ind=shuffledBirdlist[ix];

	if (IsPatchAvailable(ind,(*ind).GetX(),(*ind).GetY())==true)
	{
		//if paired, make sure that partners approx. share foraging/feeding:
		int erequired_denom=1;
		if ((*ind).IsPaired()==true) 
		{
			list<Stork>::iterator partner;
			for (vector<list<Stork>::iterator>::size_type iy=0;iy!=shuffledBirdlist.size();iy++)
			{
				partner=shuffledBirdlist[iy];

				if ((*partner).GetId()==(*ind).GetMateID()) 
				{
					if (iy>ix) {erequired_denom=2;}	
					break;
				}
			}	 
		}

		//foraging
		while (((*ind).GetForagingtime()>0)&&((*ind).GetEGained()<((double) (*ind).GetERequired()/erequired_denom))&&(IsPatchAvailable(ind,(*ind).GetX(),(*ind).GetY())==true))
		{
			if (forageChoice==1) {if (rand01()<prob_randomPatch) RandomPatch(ind,(*ind).GetX(),(*ind).GetY()); 
				else OptimalPatch(ind,(*ind).GetX(),(*ind).GetY());}
			else if (forageChoice==2) LotteryPatch(ind,(*ind).GetX(),(*ind).GetY());

			if (GPStag) (*ind).GPSdata.AddGPSday(day);
		}

	//-------------------------------
	//now truely decrease resources
	for (int xi=0;xi<max_x;xi++)
		for (int yi=0;yi<max_y;yi++)
		{
			Grid[xi][yi]->RememberResource();
			Grid[xi][yi]->DecrTrueResource();
			Grid[xi][yi]->RestoreResource();
		}

	}}

	shuffledBirdlist.clear();
}

//----------------------------------------------------------------------------------
//random patch selection
void Colony::RandomPatch(list<Stork>::iterator ind, int x, int y)
{
	int  dist,xnew=x,ynew=y;
	double angle, ran, flighttime,flightx, tau=0.0, tau_f, tau_pot, enet=0.0,resource_lost;
	bool foraged=true;

	Parameter *param=0;
	
	//-----------------------------------------------
	do {
	//choose random direction and distance (with prob_nearby that nearby patch is chosen)
	ran=rand01();
	if (ran<=prob_nearby) dist=randMinMax(0,distnearby); else dist=randMinMax(distnearby+1,distmax); 

	if (dist<=1) flighttime=0;	//when foraging in central home range cell or 1 next to it then no time costs for flight
	else {flightx = ((dist/2)-0.5); flighttime = tf * flightx;}	//divide dist by two to obtain km instead of number of cells
	
	if (dist>0)
	{
	angle=2*PI*randMinMax(1,360)/360;

	//get new x and y coordinates
	if (sin(angle)>=0) xnew+=(int) ceil(sin(angle)*dist); else xnew+=(int) floor(sin(angle)*dist);
	if (cos(angle)>=0) ynew+=(int) ceil(cos(angle)*dist); else ynew+=(int) floor(cos(angle)*dist);

	//cyclic boundary conditions (through remainder %max_x resp. %max_y)
	if (xnew<0) xnew=(xnew+max_x)%max_x;
	else if (xnew>=max_x) xnew=xnew%max_x;
	if (ynew<0) ynew=(ynew+max_y)%max_y;
	else if (ynew>=max_y) ynew=ynew%max_y;
	}
	} while (Grid[xnew][ynew]->GetResource()==0); // stork is not completely stupid but only choose random patch from patch with resources > 0

	tau_pot = (*ind).GetForagingtime()-flighttime; 
	if (tau_pot>0)	//enough foraging time left for flight?
	{
	//Foraging time tau to cover flight costs:
	if (Grid[xnew][ynew]->GetResource()>0) tau_f = dist * ef / Grid[xnew][ynew]->GetResource();

	if ((Grid[xnew][ynew]->GetResource()>0)&&(tau_f<=tau_pot))	//enough foraging time left?
	{
	//set tau_max
	param = new Parameter;
	param->Set_tau_max(Grid[xnew][ynew]->GetMaxResource()/emax,(*ind).GetNLowquality());

	//net energy gain is zero after accounting for flight costs - need to further increase tau for enet>0
	enet =0.0;
	tau=tau_f;
	//increase tau until tau_max (resp. tau_pot) or until energy requirements of adults and nestlings covered
	while ((tau<=(tau_max-0.01))&&(tau<=(tau_pot-0.01))&&(enet < ((*ind).GetERequired() - (*ind).GetEGained())) )
	{
		tau+=0.01;
		enet = Grid[xnew][ynew]->GetResource() * tau - (dist * ef);
	}

	(*ind).DecreaseForagingTime(flighttime);	//substract flighttime from available foraging time of day
	(*ind).DecreaseForagingTime(tau);
	(*ind).IncrEGained(enet);

	//remember if low quality patch was visited
	if (Grid[xnew][ynew]->GetMaxResource()/emax<0.65) (*ind).IncrNLowquality(); else (*ind).InitNLowquality();

	delete param;
	param=0;

	//foraging trip successful
	} 	//foraging trip failed - not enough time left for foraging or no resources found
	else 
	{
		(*ind).DecreaseForagingTime(flighttime); 
		if (Grid[xnew][ynew]->GetResource()<=0.0) {(*ind).IncrERequired(dist*ef);(*ind).IncrNLowquality();}
		else if (tau_f>tau_pot) {tau=tau_pot; (*ind).IncrERequired((tau_f-tau_pot)*Grid[xnew][ynew]->GetResource());}
	}
	
	resource_lost=Grid[xnew][ynew]->GetResource()*tau/hours_prey;	//resource taken by bird averaged over daytime hours
	Grid[xnew][ynew]->DecrResource((double) resource_lost);

	if ((GPStag)&&((*ind).GetStatus()>0)) {(*ind).GPSdata.Addflighttime(flighttime); 
		(*ind).GPSdata.AddAvailableEnergy(Grid[xnew][ynew]->GetResource()); (*ind).GPSdata.AddConsumedEnergy(enet);}

	if ((MemoryMap)&&((*ind).GetStatus()>0)) (*ind).memorymap[xnew][ynew]+=1*memoryweight;

	} 	//foraging trip failed - not enough flighttime left
	else 
	{
		double time, disttravelled;

		time=(*ind).GetForagingtime();
		disttravelled=(time/tf+0.5)*2;

		(*ind).DecreaseForagingTime(time);(*ind).IncrERequired(disttravelled*ef);

		if ((GPStag)&&((*ind).GetStatus()>0)) {(*ind).GPSdata.Addflighttime(time);
			(*ind).GPSdata.AddAvailableEnergy(0.0); (*ind).GPSdata.AddConsumedEnergy(0.0);}
	}

	//if paired synchronise partner's energy levels (because they share feeding of nestlings):
	if ((*ind).IsPaired()==true) 
	{
		list<Stork>::iterator partner;
		for (partner = Birdlist.begin(); partner != Birdlist.end(); ++partner)
		{
			if ((*partner).GetId()==(*ind).GetMateID()) break;
		}
		(*partner).SyncEs((*ind).GetERequired(),(*ind).GetEConsumed(),(*ind).GetEGained());
	}

	if ((GPStag)&&((*ind).GetStatus()>0)) {(*ind).GPSdata.AddGPSx(xnew);(*ind).GPSdata.AddGPSy(ynew);
		(*ind).GPSdata.AddGPStime((*ind).GetForagingtimeperday()-(*ind).GetForagingtime());}

}

//----------------------------------------------------------------------------------
//optimal patch selection
void Colony::OptimalPatch(list<Stork>::iterator ind, int x, int y)
{
	int dist,distmax_opt,xopt=x,yopt=y,xnew,ynew,memory_denom=0;
	double epot_pertime,epot,eopt=0.0,eopt_pertime=0.0,flightx, flighttime,flighttime_opt=0.0,tau, 
		tau_pot,tau_opt=0.0,resource_lost, weight=1.0, weight_opt=1.0;
	bool better=false;

	if ((MemoryMap)&&((*ind).GetStatus()>0))
	{
	for (int xi=0;xi<max_x;xi++)
		for (int yi=0;yi<max_y;yi++)
		{
			memory_denom+=(*ind).memorymap[xi][yi];
		}
	}

	Parameter *param=0;
	if (DistmaxSetting==1) distmax_opt=distmax;
	else if (DistmaxSetting==2) distmax_opt=min( min(max_x/2,max_y/2), (int) ceil(min(((*ind).GetForagingtime()/tf+0.5)*2,tau_max_good*emax/ef)));
	
	//------------------------------------------------
	//search home range for optimal foraging patch
	for (int xnewt=(x-distmax_opt);xnewt<=(x+distmax_opt);xnewt++)
		for (int ynewt=(y-distmax_opt);ynewt<=(y+distmax_opt);ynewt++)
		{
			//cyclic boundary conditions
			xnew=xnewt;
			ynew=ynewt;
			if (xnew<0) xnew=(xnew+max_x)%max_x;
			if (xnew>=max_x) xnew=xnew%max_x;
			if (ynew<0) ynew=(ynew+max_y)%max_y;
			if (ynew>=max_y) ynew=ynew%max_y;

			dist=(int) floor(sqrt((pow((double) (x-xnewt),2)+pow((double) (y-ynewt),2))));
			if (dist<=distmax_opt)
			{
				if ((MemoryMap)&&((*ind).GetStatus()>0)) weight=(double) (*ind).memorymap[xnew][ynew]/memory_denom;

				if (dist<=1) flighttime=0; else {flightx=(dist/2)-0.5;flighttime = tf* flightx;}

				tau_pot = (*ind).GetForagingtime()-flighttime; 
				if (tau_pot>0)	//enough foraging time left for flight?
				{
					//Foraging time tau:
					if (Grid[xnew][ynew]->GetResource()>0) tau = dist * ef / Grid[xnew][ynew]->GetResource();

					if ((Grid[xnew][ynew]->GetResource()>0)&&(tau<tau_pot))	//enough foraging time left?
					{
					//set tau_max
					param = new Parameter;
					param->Set_tau_max(Grid[xnew][ynew]->GetMaxResource()/emax,(*ind).GetNLowquality());

					//potential energy gain is zero after accounting for flight costs
					epot = 0.0;
					//increase tau until tau_max (resp. tau_pot) or until energy requirements of adults and nestlings covered
					while ((tau<=(tau_max-0.01))&&(tau<=(tau_pot-0.01))&&(epot < ((*ind).GetERequired() - (*ind).GetEGained())) )
					{
						tau+=0.01;
						epot = Grid[xnew][ynew]->GetResource() * tau - (dist * ef);
					}

					epot_pertime = epot / (tau + flighttime);
					
					if ((economy==1)&&(epot_pertime*weight>eopt_pertime*weight_opt)) better=true; //time-minimising
					else if ((economy==2)&&(epot*weight>eopt*weight_opt)) better=true; //energy-maximising
					else better=false;

					if (better==true) {eopt=epot; eopt_pertime=epot_pertime; xopt=xnew; yopt=ynew; 
						flighttime_opt=flighttime;tau_opt=tau;weight_opt=weight;}

					delete param;
					param=0;
					}
				}
			}
		}	//end loop over home range cells

		(*ind).DecreaseForagingTime(flighttime_opt);	//substract flighttime from available foraging time of day
		(*ind).DecreaseForagingTime(tau_opt);
		(*ind).IncrEGained(eopt);

		//if paired also increase partner's net energy gain (because they share feeding of nestlings):
		if ((*ind).IsPaired()==true) 
		{
			list<Stork>::iterator partner;
			for (partner = Birdlist.begin(); partner != Birdlist.end(); ++partner)
			{
				if ((*partner).GetId()==(*ind).GetMateID()) break;
			}
			(*partner).IncrEGained(eopt);
		}

		resource_lost=Grid[xopt][yopt]->GetResource()*tau_opt/hours_prey;	//resource taken by bird averaged over daytime hours
		Grid[xopt][yopt]->DecrResource((double) resource_lost);

		//remember if low quality patch was visited
		if (Grid[xopt][yopt]->GetMaxResource()/emax<0.65) (*ind).IncrNLowquality(); else (*ind).InitNLowquality();

		if ((GPStag)&&((*ind).GetStatus()>0)) {(*ind).GPSdata.AddGPSx(xopt);(*ind).GPSdata.AddGPSy(yopt);
			(*ind).GPSdata.AddGPStime((*ind).GetForagingtimeperday()-(*ind).GetForagingtime());(*ind).GPSdata.Addflighttime(flighttime_opt);
			(*ind).GPSdata.AddAvailableEnergy(Grid[xopt][yopt]->GetResource()); (*ind).GPSdata.AddConsumedEnergy(eopt);}

		if ((MemoryMap)&&((*ind).GetStatus()>0)) (*ind).memorymap[xopt][yopt]+=1*memoryweight;
}

//----------------------------------------------------------------------------------
//lottery patch selection
void Colony::LotteryPatch(list<Stork>::iterator ind, int x, int y)
{
	int dist,distmax_opt,xopt=x,yopt=y,xnew,ynew,cell_opt=0, memory_denom=0;
	double flightx, flighttime,epot,tau, tau_pot,resource_lost,all_e=0.0, weight=1.0;
	bool better=false;

	Parameter *param=0;
	if (DistmaxSetting==1) distmax_opt=distmax;
	else if (DistmaxSetting==2) distmax_opt=min( min(max_x/2,max_y/2),(int) ceil(min(((*ind).GetForagingtime()/tf+0.5)*2,tau_max_good*emax/ef)));

	int total_cellnumb=(distmax_opt*2+1)*(distmax_opt*2+1);

	vector<double> Prob_PatchSelect(total_cellnumb,0.0),epot_pertime(total_cellnumb,0.0),epot_abs(total_cellnumb,0.0),
		flighttime_all(total_cellnumb,0.0),tau_all(total_cellnumb,0.0);
	vector<int> coord_x,coord_y;

	for (int xnewt=(x-distmax_opt);xnewt<=(x+distmax_opt);xnewt++)
		for (int ynewt=(y-distmax_opt);ynewt<=(y+distmax_opt);ynewt++)
		{		
			coord_x.push_back(xnewt);
			coord_y.push_back(ynewt);
		}
	
	//------------------------------------------------
	//search home range for optimal foraging patch
	for (int cell=0; cell < total_cellnumb; cell++)
	{	
		//cyclic boundary conditions
		xnew=coord_x[cell];
		ynew=coord_y[cell];
		if (xnew<0) xnew=(xnew+max_x)%max_x;
		if (xnew>=max_x) xnew=xnew%max_x;
		if (ynew<0) ynew=(ynew+max_y)%max_y;
		if (ynew>=max_y) ynew=ynew%max_y;

		dist=(int) floor(sqrt((pow((double) (x-coord_x[cell]),2)+pow((double) (y-coord_y[cell]),2))));
		if (dist<=distmax_opt)
		{
			if (dist<=1) flighttime=0; else {flightx=(dist/2)-0.5;flighttime = tf* flightx;}

			tau_pot = (*ind).GetForagingtime()-flighttime; 
			if (tau_pot>0)	//enough foraging time left for flight?
			{
				//Foraging time tau:
				if (Grid[xnew][ynew]->GetResource()>0) tau = dist * ef / Grid[xnew][ynew]->GetResource();

				if ((Grid[xnew][ynew]->GetResource()>0)&&(tau<tau_pot))	//enough foraging time left?
				{
				//set tau_max
				param = new Parameter;
				param->Set_tau_max(Grid[xnew][ynew]->GetMaxResource()/emax,(*ind).GetNLowquality());

				//potential energy gain is zero after accounting for flight costs
				epot = 0.0;
				//increase tau until tau_max (resp. tau_pot) or until energy requirements of adults and nestlings covered
				while ((tau<=(tau_max-0.01))&&(tau<=(tau_pot-0.01))&&(epot < ((*ind).GetERequired() - (*ind).GetEGained())) )
				{
					tau+=0.01;
					epot = Grid[xnew][ynew]->GetResource() * tau - (dist * ef);
				}

				epot_abs[cell] = epot;
				epot_pertime[cell] = epot / (tau + flighttime);
				flighttime_all[cell] = flighttime;
				tau_all[cell] = tau;

				delete param;
				param=0;
				}
			}
		}
	}	//end loop over home range cells

	if ((MemoryMap) &&((*ind).GetStatus()>0))
	{
		for (int xi=0;xi<max_x;xi++)
			for (int yi=0;yi<max_y;yi++)
			{
				memory_denom+=(*ind).memorymap[xi][yi];
			}
	}

	// update coordinates wrt cyclic boundary conditions
	for (int cell=0; cell < total_cellnumb; cell++)
	{
		xnew=coord_x[cell];
		ynew=coord_y[cell];
		if (xnew<0) xnew=(xnew+max_x)%max_x;
		if (xnew>=max_x) xnew=xnew%max_x;
		if (ynew<0) ynew=(ynew+max_y)%max_y;
		if (ynew>=max_y) ynew=ynew%max_y;

		coord_x[cell]=xnew;
		coord_y[cell]=ynew;
	}

	// count all potential energy
	for (int cell=0; cell < total_cellnumb; cell++)
	{
		if ((MemoryMap)&&((*ind).GetStatus()>0)) weight=(double) (*ind).memorymap[coord_x[cell]][coord_y[cell]]/memory_denom;

		if (economy==1) all_e += epot_pertime[cell]*weight;	//time-minimising
		else if (economy==2) all_e += epot_abs[cell]*weight;	//energy-maximising
	}

	// calculate cumulative probabilities
	if ((MemoryMap)&&((*ind).GetStatus()>0)) weight=(double) (*ind).memorymap[coord_x[0]][coord_y[0]]/memory_denom;
	if (economy==1) Prob_PatchSelect[0] = epot_pertime[0]*weight/all_e;	//time-minimising
	else if (economy==2) Prob_PatchSelect[0] = epot_abs[0]*weight/all_e;	//energy-maximising

	for (int cell=1; cell < total_cellnumb; cell++)
	{
		if ((MemoryMap)&&((*ind).GetStatus()>0)) weight=(double) (*ind).memorymap[coord_x[cell]][coord_y[cell]]/memory_denom;

		if (economy==1) Prob_PatchSelect[cell] = epot_pertime[cell]*weight/all_e + Prob_PatchSelect[cell-1];	//time-minimising
		else if (economy==2) Prob_PatchSelect[cell] = epot_abs[cell]*weight/all_e + Prob_PatchSelect[cell-1];	//energy-maximising
	}

	// to make sure that all probs really sum up to 1.0
	for (int cell=1; cell < total_cellnumb; cell++)
	{
		Prob_PatchSelect[cell]=Prob_PatchSelect[cell]/Prob_PatchSelect[total_cellnumb-1];
	}

	//draw from probability distribution
	double randnumb = rand01();
	for (int cell=0; cell < total_cellnumb; cell++)
	{
		if (randnumb < Prob_PatchSelect[cell]) {cell_opt=cell; break;}
	}

	(*ind).DecreaseForagingTime(flighttime_all[cell_opt]);	//substract flighttime from available foraging time of day
	(*ind).DecreaseForagingTime(tau_all[cell_opt]);
	(*ind).IncrEGained(epot_abs[cell_opt]);

	//if paired also increase partner's net energy gain (because they share feeding of nestlings):
	if ((*ind).IsPaired()==true) 
	{
		list<Stork>::iterator partner;
		for (partner = Birdlist.begin(); partner != Birdlist.end(); ++partner)
		{
			if ((*partner).GetId()==(*ind).GetMateID()) break;
		}
		(*partner).IncrEGained(epot_abs[cell_opt]);
	}

	xopt=coord_x[cell_opt];
	yopt=coord_y[cell_opt];

	resource_lost=Grid[xopt][yopt]->GetResource()*tau_all[cell_opt]/hours_prey;	//resource taken by bird averaged over daytime hours
	Grid[xopt][yopt]->DecrResource((double) resource_lost);

	//remember if low quality patch was visited
	if (Grid[xopt][yopt]->GetMaxResource()/emax<0.65) (*ind).IncrNLowquality(); else (*ind).InitNLowquality();

	if ((GPStag)&&((*ind).GetStatus()>0)) {(*ind).GPSdata.AddGPSx(xopt);(*ind).GPSdata.AddGPSy(yopt);
		(*ind).GPSdata.AddGPStime((*ind).GetForagingtimeperday()-(*ind).GetForagingtime());(*ind).GPSdata.Addflighttime(flighttime_all[cell_opt]);
		(*ind).GPSdata.AddAvailableEnergy(Grid[xopt][yopt]->GetResource()); (*ind).GPSdata.AddConsumedEnergy(epot_abs[cell_opt]);}

	if ((MemoryMap)&&((*ind).GetStatus()>0)) (*ind).memorymap[xopt][yopt]+=1*memoryweight;
}


//---------------------------------------------------------------------------------
//use optimal patch search algorithm to test whether any foraging patch left within home range
bool Colony::IsPatchAvailable(list<Stork>::iterator ind, int x, int y)
{
	int dist,distmax_opt,xnew,ynew;
	double epot,flighttime,tau, tau_pot;
	bool better=false, patch_left=false;

	Parameter *param=0;
	if (DistmaxSetting==1) distmax_opt=distmax;
	else if (DistmaxSetting==2) distmax_opt=min( min(max_x/2,max_y/2),(int) ceil(min(((*ind).GetForagingtime()/tf+0.5)*2,tau_max_good*emax/ef)));

	//------------------------------------------------
	//search home range for optimal foraging patch
	for (int xnewt=(x-distmax_opt);xnewt<=(x+distmax_opt);xnewt++)
		for (int ynewt=(y-distmax_opt);ynewt<=(y+distmax_opt);ynewt++)
		{
			if (patch_left==false)
			{
			//cyclic boundary conditions
			xnew=xnewt;
			ynew=ynewt;
			if (xnewt<0) xnew=(xnewt+max_x)%max_x;
			if (xnewt>=max_x) xnew=xnewt%max_x;
			if (ynewt<0) ynew=(ynewt+max_y)%max_y;
			if (ynewt>=max_y) ynew=ynewt%max_y;


			dist=(int) floor(sqrt((pow((double) (x-xnewt),2)+pow((double) (y-ynewt),2))));
			if (dist<=distmax_opt)
			{
				if (dist<=1) flighttime=0; else flighttime = tf*(dist/2-0.5);

				tau_pot = (*ind).GetForagingtime()-flighttime; 
				if (tau_pot>0)	//enough foraging time left for flight?
				{
					//Foraging time tau:
					if (Grid[xnew][ynew]->GetResource()>0) tau = (dist * ef) / (Grid[xnew][ynew]->GetResource());

					if ((Grid[xnew][ynew]->GetResource()>0)&&(tau<tau_pot))	//enough foraging time left?
					{
					//set tau_max
					param = new Parameter;
					param->Set_tau_max(Grid[xnew][ynew]->GetMaxResource()/emax,(*ind).GetNLowquality());

					//potential energy gain is zero after accounting for flight costs
					epot = 0.0;
					//increase tau until tau_max (resp. tau_pot) or until energy requirements of adults and nestlings covered
					while ((tau<=(tau_max-0.01))&&(tau<=(tau_pot-0.01))&&(epot < ((*ind).GetERequired() - (*ind).GetEGained())) )
					{
						tau+=0.01;
						epot = Grid[xnew][ynew]->GetResource() * tau - (dist * ef);
					}

					if (epot>0) patch_left=true;
					delete param;
					param=0;
					}
				}
			}
			}}	//end loop over home range cells

		return(patch_left);
}