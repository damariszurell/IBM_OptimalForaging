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


#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "Simulation.h"
#include "parameter.h"
#include "Cell.h"
#include "Stork.h"
#include "Colony.h"
#include "functions.h"

using namespace std;

//-----------------------------------------------------------------
//Constructor
Simulation::Simulation(void)
{
}

//-----------------------------------------------------------------
//Destructor
Simulation::~Simulation(void)
{
}

//-----------------------------------------------------------------
void Simulation::SingleHomerange(void)
{
	Parameter *param=0;
	param = new Parameter; 
	
	Colony *pop=0;
	pop = new Colony;

	list<Stork>::iterator singlebird;

	//set maximum number of birds to one!
	param->Set_max_birds(1);

	pop->InitBirds(); 
	singlebird = pop->Birdlist.begin();
	std::cout << (*singlebird).GetX() <<" " << (*singlebird).GetY() <<std::endl;

	for (int t=0;t<5;t++)
	{
	pop->Forage(t);
	pop->GoodNight();
	}
	
	int num=(*singlebird).GPSdata.GetGPSsize();

	for (int i=0;i<num;i++)
		std::cout << (*singlebird).GPSdata.GetGPSday()[i] << " " << (*singlebird).GPSdata.GetGPSx()[i] << " " 
		<< (*singlebird).GPSdata.GetGPSy()[i] << " " << (*singlebird).GPSdata.GetGPStime()[i] << " " 
		<< (*singlebird).GPSdata.Getflighttime()[i] << " " 
		<< (*singlebird).GPSdata.GetAvailableEnergy()[i] << " "
		<< (*singlebird).GPSdata.GetConsumedEnergy()[i] << std::endl;



	//-------------------
	delete pop;
	pop=0;

	delete param;
	param=0;
}

//-----------------------------------------------------------------
//Run simulation for homerange establishment only
void Simulation::HomerangeEstablishment(void)
{
	//what would be interesting to know?
	//how many birds and respective home range centers
	//how many pairs and respective home range centers
	//
	
	Colony *pop=0;
	pop = new Colony;

	int pairs=0;

	list<Stork>::iterator ind;
	//------------------

	pop->InitBirds(); 
	std::cout << pop->Birdlist.size()<<std::endl;

	for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
	{
		std::cout << (*ind).GetId() <<" x="<< (*ind).GetX() <<" y="<< (*ind).GetY() <<" sex="<< (*ind).GetSex() 
			<<" paired="<< (*ind).IsPaired() <<" mate="<< (*ind).GetMateID()<<std::endl;

		if ((*ind).IsPaired()==true) pairs++;
	}

	
	std::cout << pairs/2 <<" pairs" << std::endl;

	//-------------------
	delete pop;
	pop=0;
}

//-----------------------------------------------------------------
//Run simulation for homerange establishment and colony development over breeding season
void Simulation::BreedingPeriod(void)
{
	//std::ofstream outfile ("e:/Projects/DIP_WhiteStork/breeding/CiconiaBreeding/OUTPUT/protocol.txt");
	std::ofstream outfile ("OUTPUT/protocol_K.txt");
	std::ofstream xyfile ("OUTPUT/protocol_K_resources.txt");

	outfile<<"time\t"<<"x\t"<<"y\t"<<"id\t"<<"mate\t"<<"pair\t"<<"JZG\n";
	xyfile<<"time\t"<<"x\t"<<"y\t"<<"max_resource\t"<<"resource\n";

	Colony *pop=0;
	pop = new Colony;

	int pairs=0;

	// define list iterator
	list<Stork>::iterator ind;
	//------------------

	pop->InitBirds(); 

	for (int t=0;t<spinup;t++)
	{
		pop->Forage(t);

		//------------------------
		// write resources to file
		if (t==(spinup-1))
		{
			//loop over all grid cells
			for (int x=0;x<max_x;x++)
				for (int y=0;y<max_y;y++)
				{
					xyfile<<"prebreeding\t"<< x <<"\t"<< y <<"\t"<<
					pop->Grid[x][y]->GetMaxResource() << "\t"<< 
					pop->Grid[x][y]->GetResource() <<"\n";
				}
		}
		//------------------------

		pop->GoodNight();
	}

	// write to file
	for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
		{
			if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID()))
			{
			pairs++; 
			outfile <<"prebreeding\t"<<  (*ind).GetX() <<"\t"<< (*ind).GetY() <<"\t"<< (*ind).GetId() <<"\t"<<
				(*ind).GetMateID()<< "\t" << pairs <<"\t"<< (*ind).GetChicks() <<std::endl;
			}
			else if ((*ind).IsPaired()==false)
			{
				outfile <<"prebreeding\t"<< (*ind).GetX() <<"\t"<< (*ind).GetY() <<"\t"<< (*ind).GetId() <<"\t"<< 
				"0\t" << "NA\t"<< "NA" <<std::endl;
			}
		}
		pairs=0;

	//-------------------
	pop->InitBreeding();

	for (int t=0;t<max_days;t++)
	{
		pop->Forage(t);

		//------------------------
		// write resources to file
		if (t==(spinup-1))
		{
			//loop over all grid cells
			for (int x=0;x<max_x;x++)
				for (int y=0;y<max_y;y++)
				{
					xyfile<<"postbreeding\t"<< x <<"\t"<< y <<"\t"<<
					pop->Grid[x][y]->GetMaxResource() << "\t"<< 
					pop->Grid[x][y]->GetResource() <<"\n";
				}
		}
		//------------------------

		pop->GoodNight();
	}

	// write to file
	for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
		{
			if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID()))
			{
			pairs++; 
			outfile <<"postbreeding\t"<< (*ind).GetX() <<"\t"<< (*ind).GetY() <<"\t"<< (*ind).GetId() <<"\t"<<
				(*ind).GetMateID()<< "\t" << pairs <<"\t"<< (*ind).GetChicks() <<std::endl;
			}
			else if ((*ind).IsPaired()==false)
			{
				outfile <<"postbreeding\t"<< (*ind).GetX() <<"\t"<< (*ind).GetY() <<"\t"<< (*ind).GetId() <<"\t"<<
				"0\t" << "NA\t"<< "NA" <<std::endl;
			}
		}
		pairs=0;


	//-------------------
	delete pop;
	pop=0;

	outfile.close();
}

//---------------------------------------------------------------------------
// run replicate simulations
void Simulation::Sim_Breeding(void)
{
	int nsims=100;
	Colony *pop=0;

	//--------------------------------------
	//prepare vectors for storing results [nsims][ndays]
	array2Dint n_ind;	
	array2Dint HPa;
	array2Dint HPm; 
	array2Dint HPo;
	array2Ddouble JZa;
	array2Ddouble JZm;

	SetArraySize2d(n_ind,nsims,1+spinup+max_days);
	SetArraySize2d(HPa,nsims,1+spinup+max_days);
	SetArraySize2d(HPm,nsims,1+spinup+max_days);
	SetArraySize2d(HPo,nsims,1+spinup+max_days);
	SetArraySize2d(JZa,nsims,1+spinup+max_days);
	SetArraySize2d(JZm,nsims,1+spinup+max_days);

	//--------------------------------------
	//prepare file output

	//std::ofstream simfile ("e:/Projects/DIP_WhiteStork/breeding/CiconiaBreeding/OUTPUT/SimBreeding.txt");
	std::ofstream simfile ("OUTPUT/SimBreeding.txt");

	simfile<<"max_birds\t"<<"day\t"<<"ind\t"<<"sd_ind\t"<<"HPa\t"<<"sd_HPa\t"<<"HPo\t"
		<<"sd_HPo\t"<<"HPm\t"<<"sd_HPm\t"<<"JZa\t"<<"sd_JZa\t"<<"JZm\t"<<"sd_JZm\n";

	simfile.close();

	//--------------------------------------
	//start simulations
	for (int sim=0;sim<nsims;sim++)
	{
		//variables
		pop = new Colony;
		int hpa=0,hpo=0,hpm=0,jz=0;
		list<Stork>::iterator ind;
		
		//--------------------------------------
		//home range search

		pop->InitBirds();

		//count numbers
		for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
		{
			if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
			{
				hpa++;
				if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
				else if ((*ind).GetChicks()==0) hpo++;
			}
		}

		n_ind[sim][0]=(int) pop->Birdlist.size();
		HPa[sim][0]=hpa;
		HPm[sim][0]=hpm;
		HPo[sim][0]=hpo;
		JZa[sim][0]=(double) jz/hpa;
		JZm[sim][0]=(double) jz/hpm;
		hpa=0;hpo=0;hpm=0;jz=0;

		//--------------------------------
		// spinup / burnin

		for (int t=0;t<spinup;t++)
		{
			pop->Forage(t);
			pop->GoodNight();

			//count numbers
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
				{
					hpa++;
					if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
					else if ((*ind).GetChicks()==0) hpo++;
				}
			}

			n_ind[sim][t+1]=(int) pop->Birdlist.size();
			HPa[sim][t+1]=hpa;
			HPm[sim][t+1]=hpm;
			HPo[sim][t+1]=hpo;
			JZa[sim][t+1]=(double) jz/hpa;
			JZm[sim][t+1]=(double) jz/hpm;
			hpa=0;hpo=0;hpm=0;jz=0;
		}

		//--------------------------------
		// initiate breeding and simulate breeding period

		pop->InitBreeding();

		for (int t=0;t<max_days;t++)
		{
			pop->Forage(t);
			pop->GoodNight();

			//count numbers
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
				{
					hpa++;
					if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
					else if ((*ind).GetChicks()==0) hpo++;
				}
			}

			n_ind[sim][t+1+spinup]=(int) pop->Birdlist.size();
			HPa[sim][t+1+spinup]=hpa;
			HPm[sim][t+1+spinup]=hpm;
			HPo[sim][t+1+spinup]=hpo;
			JZa[sim][t+1+spinup]=(double) jz/hpa;
			JZm[sim][t+1+spinup]=(double) jz/hpm;
			hpa=0;hpo=0;hpm=0;jz=0;
		}

		//------------------------------------
		// free dynamic memory
		delete pop;
		pop=0;
	}
	
	//--------------------------------------
	//write to output file
	simfile.open("OUTPUT/SimBreeding.txt", ofstream::app);

	for (int t=0;t<66;t++)
	{
		simfile << 100000 <<"\t"<< t-6 <<"\t"<< mean(n_ind,nsims,t) <<"\t"<< sd(n_ind,nsims,t) <<"\t"<<
			mean(HPa,nsims,t) <<"\t"<< sd(HPa,nsims,t) <<"\t"<<
			mean(HPo,nsims,t) <<"\t"<< sd(HPo,nsims,t) <<"\t"<<
			mean(HPm,nsims,t) <<"\t"<< sd(HPm,nsims,t) <<"\t"<<
			max(0.0,mean(JZa,nsims,t)) <<"\t"<< max(0.0,sd(JZa,nsims,t)) <<"\t"<<
			max(0.0,mean(JZm,nsims,t)) <<"\t"<< max(0.0,sd(JZm,nsims,t)) <<"\n";
	}

	simfile.close();
}

//---------------------------------------------------------------------------
// run replicate simulations
void Simulation::Sim_Density(void)
{
	bool run_dens=true;
	int densquant=0; //for density manipulation: 0=all quantiles;8= N/K=0.1
	int nsims=100;
	int max_ncell=max_x*max_y;
	Colony *pop=0;
	Parameter *param=0;
	std::string DD_file="OUTPUT/SimDDpt_PS.txt";	
	std::string DD_xyfile="OUTPUT/SimDDpt_xy_PS.txt";
	std::string DD_xyfile_init="OUTPUT/SimDDpt_xyInit_PS.txt";

	//--------------------------------------
	//prepare vectors for storing results [nsims][ndays]
	array2Dint n_ind;
	array2Dint HPa;
	array2Dint HPm;
	array2Dint HPo;
	array2Dint HE;
	array2Dint JZG;
	array2Ddouble JZa;
	array2Ddouble JZm;

	SetArraySize2d(n_ind,nsims,1+spinup+max_days);
	SetArraySize2d(HPa,nsims,1+spinup+max_days);
	SetArraySize2d(HPm,nsims,1+spinup+max_days);
	SetArraySize2d(HPo,nsims,1+spinup+max_days);
	SetArraySize2d(HE,nsims,1+spinup+max_days);
	SetArraySize2d(JZG,nsims,1+spinup+max_days);
	SetArraySize2d(JZa,nsims,1+spinup+max_days);
	SetArraySize2d(JZm,nsims,1+spinup+max_days);

	//-----------------
	//vectors for storing spatial results xy_init
	vecdouble n_ind_xy_init;
	vecdouble HPa_xy_init;
	vecdouble HPm_xy_init;
	vecdouble HPo_xy_init;
	vecdouble HE_xy_init;
	vecdouble JZG_xy_init;
	vecdouble JZa_xy_init;
	vecdouble JZm_xy_init;
	vecdouble n_occ_xy_init;

	vecint min_ind_xy_init;
	vecint min_HPa_xy_init;
	vecint min_HPm_xy_init;
	vecint min_HPo_xy_init;
	vecint min_HE_xy_init;
	vecint min_JZG_xy_init;
	vecdouble min_JZa_xy_init;
	vecdouble min_JZm_xy_init;

	vecint max_ind_xy_init;
	vecint max_HPa_xy_init;
	vecint max_HPm_xy_init;
	vecint max_HPo_xy_init;
	vecint max_HE_xy_init;
	vecint max_JZG_xy_init;
	vecdouble max_JZa_xy_init;
	vecdouble max_JZm_xy_init;

	SetVectorSize(n_ind_xy_init,max_ncell);
	SetVectorSize(HPa_xy_init,max_ncell);
	SetVectorSize(HPm_xy_init,max_ncell);
	SetVectorSize(HPo_xy_init,max_ncell);
	SetVectorSize(HE_xy_init,max_ncell);
	SetVectorSize(JZG_xy_init,max_ncell);
	SetVectorSize(JZa_xy_init,max_ncell);
	SetVectorSize(JZm_xy_init,max_ncell);
	SetVectorSize(n_occ_xy_init,max_ncell);

	SetVectorSize(min_ind_xy_init,max_ncell);
	SetVectorSize(min_HPa_xy_init,max_ncell);
	SetVectorSize(min_HPm_xy_init,max_ncell);
	SetVectorSize(min_HPo_xy_init,max_ncell);
	SetVectorSize(min_HE_xy_init,max_ncell);
	SetVectorSize(min_JZG_xy_init,max_ncell);
	SetVectorSize(min_JZa_xy_init,max_ncell);
	SetVectorSize(min_JZm_xy_init,max_ncell);

	SetVectorSize(max_ind_xy_init,max_ncell);
	SetVectorSize(max_HPa_xy_init,max_ncell);
	SetVectorSize(max_HPm_xy_init,max_ncell);
	SetVectorSize(max_HPo_xy_init,max_ncell);
	SetVectorSize(max_HE_xy_init,max_ncell);
	SetVectorSize(max_JZG_xy_init,max_ncell);
	SetVectorSize(max_JZa_xy_init,max_ncell);
	SetVectorSize(max_JZm_xy_init,max_ncell);

	// vectors for storing spatial results xy
	vecdouble n_ind_xy;
	vecdouble HPa_xy;
	vecdouble HPm_xy;
	vecdouble HPo_xy;
	vecdouble HE_xy;
	vecdouble JZG_xy;
	vecdouble JZa_xy;
	vecdouble JZm_xy;
	vecdouble n_occ_xy;

	vecint min_ind_xy;
	vecint min_HPa_xy;
	vecint min_HPm_xy;
	vecint min_HPo_xy;
	vecint min_HE_xy;
	vecint min_JZG_xy;
	vecdouble min_JZa_xy;
	vecdouble min_JZm_xy;

	vecint max_ind_xy;
	vecint max_HPa_xy;
	vecint max_HPm_xy;
	vecint max_HPo_xy;
	vecint max_HE_xy;
	vecint max_JZG_xy;
	vecdouble max_JZa_xy;
	vecdouble max_JZm_xy;

	SetVectorSize(n_ind_xy,max_ncell);
	SetVectorSize(HPa_xy,max_ncell);
	SetVectorSize(HPm_xy,max_ncell);
	SetVectorSize(HPo_xy,max_ncell);
	SetVectorSize(HE_xy,max_ncell);
	SetVectorSize(JZG_xy,max_ncell);
	SetVectorSize(JZa_xy,max_ncell);
	SetVectorSize(JZm_xy,max_ncell);
	SetVectorSize(n_occ_xy,max_ncell);

	SetVectorSize(min_ind_xy,max_ncell);
	SetVectorSize(min_HPa_xy,max_ncell);
	SetVectorSize(min_HPm_xy,max_ncell);
	SetVectorSize(min_HPo_xy,max_ncell);
	SetVectorSize(min_HE_xy,max_ncell);
	SetVectorSize(min_JZG_xy,max_ncell);
	SetVectorSize(min_JZa_xy,max_ncell);
	SetVectorSize(min_JZm_xy,max_ncell);

	SetVectorSize(max_ind_xy,max_ncell);
	SetVectorSize(max_HPa_xy,max_ncell);
	SetVectorSize(max_HPm_xy,max_ncell);
	SetVectorSize(max_HPo_xy,max_ncell);
	SetVectorSize(max_HE_xy,max_ncell);
	SetVectorSize(max_JZG_xy,max_ncell);
	SetVectorSize(max_JZa_xy,max_ncell);
	SetVectorSize(max_JZm_xy,max_ncell);

	//--------------------------------------
	//prepare file output

	std::ofstream simfile (DD_file.c_str());

	simfile<<"max_birds\t"<<"day\t"<<"ind\t"<<"sd_ind\t"<<"HPa\t"<<"sd_HPa\t"<<"HPo\t"
		<<"sd_HPo\t"<<"HPm\t"<<"sd_HPm\t"<<"HE\t"<<"sd_HE\t"<<"JZG\t"<<"sd_JZG\t"<<"JZa\t"<<"sd_JZa\t"<<"JZm\t"<<"sd_JZm\t"
		<<"max_ind\t"<<"min_ind\t"<<"max_HPa\t"<<"min_HPa\t"<<"max_HPo\t"
		<<"min_HPo\t"<<"max_HPm\t"<<"min_HPm\t"<<"min_HE\t"<<"max_HE\t"<<"min_JZG\t"<<"max_JZG\t"
		<<"max_JZa\t"<<"min_JZa\t"<<"max_JZm\t"<<"min_JZm\n";

	simfile.close();

	//---------
	std::ofstream xy_init_simfile (DD_xyfile_init.c_str());
	xy_init_simfile<<"max_birds\t"<<"x\t"<<"y\t"<<"n_occ\t"<<"ind\t"<<"HPa\t"<<"HPo\t"
		<<"HPm\t"<<"HE\t"<<"JZG\t"<<"JZa\t"<<"JZm\t"
		<<"max_ind\t"<<"min_ind\t"<<"max_HPa\t"<<"min_HPa\t"<<"max_HPo\t"
		<<"min_HPo\t"<<"max_HPm\t"<<"min_HPm\t"<<"min_HE\t"<<"max_HE\t"<<"min_JZG\t"<<"max_JZG\t"
		<<"max_JZa\t"<<"min_JZa\t"<<"max_JZm\t"<<"min_JZm\n";

	xy_init_simfile.close();

	//--------------
	std::ofstream xy_simfile (DD_xyfile.c_str());
	xy_simfile<<"max_birds\t"<<"x\t"<<"y\t"<<"n_occ\t"<<"ind\t"<<"HPa\t"<<"HPo\t"
		<<"HPm\t"<<"HE\t"<<"JZG\t"<<"JZa\t"<<"JZm\t"
		<<"max_ind\t"<<"min_ind\t"<<"max_HPa\t"<<"min_HPa\t"<<"max_HPo\t"
		<<"min_HPo\t"<<"max_HPm\t"<<"min_HPm\t"<<"min_HE\t"<<"max_HE\t"<<"min_JZG\t"<<"max_JZG\t"
		<<"max_JZa\t"<<"min_JZa\t"<<"max_JZm\t"<<"min_JZm\n";

	xy_simfile.close();
	
	//--------------------------------------
	//start simulations - full capacity /max. density
	for (int sim=0;sim<nsims;sim++)
	{
		//variables
		pop = new Colony;
		int hpa=0,hpo=0,hpm=0,jz=0;
		list<Stork>::iterator ind;
		
		//--------------------------------------
		//home range search

		pop->InitBirds();

		//count numbers
		for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
		{
			if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
			{
				hpa++;
				if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
				else if ((*ind).GetChicks()==0) hpo++;
			}
		}

		n_ind[sim][0]=(int) pop->Birdlist.size();
		HPa[sim][0]=hpa;
		HPm[sim][0]=hpm;
		HPo[sim][0]=hpo;
		HE[sim][0]=(int) pop->Birdlist.size()-(hpa*2);
		JZG[sim][0]=jz;
		JZa[sim][0]=max(0.0,(double) jz/hpa);
		JZm[sim][0]=max(0.0,(double) jz/hpm);
		hpa=0;hpo=0;hpm=0;jz=0;

		//--------------------------------
		// spinup / burnin

		for (int t=0;t<spinup;t++)
		{
			pop->Forage(t);
			pop->GoodNight();

			//count numbers
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
				{
					hpa++;
					if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
					else if ((*ind).GetChicks()==0) hpo++;
				}
			}

			n_ind[sim][t+1]=(int) pop->Birdlist.size();
			HPa[sim][t+1]=hpa;
			HPm[sim][t+1]=hpm;
			HPo[sim][t+1]=hpo;
			HE[sim][t+1]=(int) pop->Birdlist.size()-(hpa*2);
			JZG[sim][t+1]=jz;
			JZa[sim][t+1]=max(0.0,(double) jz/hpa);
			JZm[sim][t+1]=max(0.0,(double) jz/hpm);
			hpa=0;hpo=0;hpm=0;jz=0;
		}
		//--------------------------------
		// output initial spatial distribution after home range formation
		for (int xi=0;xi<max_x;xi++)
			for (int yi=0;yi<max_y;yi++)
			{
				int n=0;
				list<Stork>::iterator ind;

				//count numbers
				for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
				{
					if (((*ind).GetX()==xi)&&((*ind).GetY()==yi)) 
					{
						n++;

						if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
						{
							hpa++;
							if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
							else if ((*ind).GetChicks()==0) hpo++;
						}
					}
				}

				if (n>0) n_occ_xy_init[xi*max_x+yi]++;
				
				// store results
				if (sim==0)
				{
					n_ind_xy_init[xi*max_x+yi]=n;min_ind_xy_init[xi*max_x+yi]=n;max_ind_xy_init[xi*max_x+yi]=n;
					HPa_xy_init[xi*max_x+yi]=hpa;min_HPa_xy_init[xi*max_x+yi]=hpa;max_HPa_xy_init[xi*max_x+yi]=hpa;
					HPm_xy_init[xi*max_x+yi]=hpm;min_HPm_xy_init[xi*max_x+yi]=hpm;max_HPm_xy_init[xi*max_x+yi]=hpm;
					HPo_xy_init[xi*max_x+yi]=hpo;min_HPo_xy_init[xi*max_x+yi]=hpo;max_HPo_xy_init[xi*max_x+yi]=hpo;
					HE_xy_init[xi*max_x+yi]=n-(hpa*2);min_HE_xy_init[xi*max_x+yi]=n-(hpa*2);max_HE_xy_init[xi*max_x+yi]=n-(hpa*2);
					JZG_xy_init[xi*max_x+yi]=jz;min_JZG_xy_init[xi*max_x+yi]=jz;max_JZG_xy_init[xi*max_x+yi]=jz;
					JZa_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpa);min_JZa_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpa);max_JZa_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpa);
					JZm_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpm);min_JZm_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpm);max_JZm_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpm);
				} else
				{
					n_ind_xy_init[xi*max_x+yi]=(n_ind_xy_init[xi*max_x+yi]+n)/2;
					HPa_xy_init[xi*max_x+yi]=(HPa_xy_init[xi*max_x+yi]+hpa)/2;
					HPm_xy_init[xi*max_x+yi]=(HPm_xy_init[xi*max_x+yi]+hpm)/2;
					HPo_xy_init[xi*max_x+yi]=(HPo_xy_init[xi*max_x+yi]+hpo)/2;
					HE_xy_init[xi*max_x+yi]=(HE_xy_init[xi*max_x+yi]+(n-(hpa*2)))/2;
					JZG_xy_init[xi*max_x+yi]=(JZG_xy_init[xi*max_x+yi]+jz)/2;
					JZa_xy_init[xi*max_x+yi]=(JZa_xy_init[xi*max_x+yi]+max(0.0,(double) jz/hpa))/2;
					JZm_xy_init[xi*max_x+yi]=(JZm_xy_init[xi*max_x+yi]+max(0.0,(double) jz/hpm))/2;	

					min_ind_xy_init[xi*max_x+yi]=min(min_ind_xy_init[xi*max_x+yi],n);
					min_HPa_xy_init[xi*max_x+yi]=min(min_HPa_xy_init[xi*max_x+yi],hpa);
					min_HPm_xy_init[xi*max_x+yi]=min(hpm,min_HPm_xy_init[xi*max_x+yi]);
					min_HPo_xy_init[xi*max_x+yi]=min(min_HPo_xy_init[xi*max_x+yi],hpo);
					min_HE_xy_init[xi*max_x+yi]=min(min_HE_xy_init[xi*max_x+yi],(n-(hpa*2)));
					min_JZG_xy_init[xi*max_x+yi]=min(min_JZG_xy_init[xi*max_x+yi],jz);
					min_JZa_xy_init[xi*max_x+yi]=min(min_JZa_xy_init[xi*max_x+yi],max(0.0,(double) jz/hpa));
					min_JZm_xy_init[xi*max_x+yi]=min(min_JZm_xy_init[xi*max_x+yi],max(0.0,(double) jz/hpm));

					max_ind_xy_init[xi*max_x+yi]=max(max_ind_xy_init[xi*max_x+yi],n);
					max_HPa_xy_init[xi*max_x+yi]=max(max_HPa_xy_init[xi*max_x+yi],hpa);
					max_HPm_xy_init[xi*max_x+yi]=max(hpm,max_HPm_xy_init[xi*max_x+yi]);
					max_HPo_xy_init[xi*max_x+yi]=max(max_HPo_xy_init[xi*max_x+yi],hpo);
					max_HE_xy_init[xi*max_x+yi]=max(max_HE_xy_init[xi*max_x+yi],(n-(hpa*2)));
					max_JZG_xy_init[xi*max_x+yi]=max(max_JZG_xy_init[xi*max_x+yi],jz);
					max_JZa_xy_init[xi*max_x+yi]=max(max_JZa_xy_init[xi*max_x+yi],max(0.0,(double) jz/hpa));
					max_JZm_xy_init[xi*max_x+yi]=max(max_JZm_xy_init[xi*max_x+yi],max(0.0,(double) jz/hpm));
				}
				n=0;hpa=0;hpo=0;hpm=0;jz=0;
			}

		//--------------------------------
		// initiate breeding and simulate breeding period

		pop->InitBreeding();

		for (int t=0;t<max_days;t++)
		{
			pop->Forage(t);
			pop->GoodNight();

			//count numbers
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
				{
					hpa++;
					if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
					else if ((*ind).GetChicks()==0) hpo++;
				}
			}

			n_ind[sim][t+1+spinup]=(int) pop->Birdlist.size();
			HPa[sim][t+1+spinup]=hpa;
			HPm[sim][t+1+spinup]=hpm;
			HPo[sim][t+1+spinup]=hpo;
			HE[sim][t+1+spinup]=(int) pop->Birdlist.size()-(hpa*2);
			JZG[sim][t+1+spinup]=jz;
			JZa[sim][t+1+spinup]=max(0.0,(double) jz/hpa);
			JZm[sim][t+1+spinup]=max(0.0,(double) jz/hpm);
			hpa=0;hpo=0;hpm=0;jz=0;
		}
		
		//----------------------------------------------------
		// loop over all grid cells to count and store spatial results
		for (int xi=0;xi<max_x;xi++)
			for (int yi=0;yi<max_y;yi++)
			{
				int n=0;
				list<Stork>::iterator ind;

				//count numbers
				for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
				{
					if (((*ind).GetX()==xi)&&((*ind).GetY()==yi)) 
					{
						n++;

						if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
						{
							hpa++;
							if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
							else if ((*ind).GetChicks()==0) hpo++;
						}
					}
				}

				if (n>0) n_occ_xy[xi*max_x+yi]++;
				
				// store results
				if (sim==0)
				{
					n_ind_xy[xi*max_x+yi]=n;min_ind_xy[xi*max_x+yi]=n;max_ind_xy[xi*max_x+yi]=n;
					HPa_xy[xi*max_x+yi]=hpa;min_HPa_xy[xi*max_x+yi]=hpa;max_HPa_xy[xi*max_x+yi]=hpa;
					HPm_xy[xi*max_x+yi]=hpm;min_HPm_xy[xi*max_x+yi]=hpm;max_HPm_xy[xi*max_x+yi]=hpm;
					HPo_xy[xi*max_x+yi]=hpo;min_HPo_xy[xi*max_x+yi]=hpo;max_HPo_xy[xi*max_x+yi]=hpo;
					HE_xy[xi*max_x+yi]=n-(hpa*2);min_HE_xy[xi*max_x+yi]=n-(hpa*2);max_HE_xy[xi*max_x+yi]=n-(hpa*2);
					JZG_xy[xi*max_x+yi]=jz;min_JZG_xy[xi*max_x+yi]=jz;max_JZG_xy[xi*max_x+yi]=jz;
					JZa_xy[xi*max_x+yi]=max(0.0,(double) jz/hpa);min_JZa_xy[xi*max_x+yi]=max(0.0,(double) jz/hpa);max_JZa_xy[xi*max_x+yi]=max(0.0,(double) jz/hpa);
					JZm_xy[xi*max_x+yi]=max(0.0,(double) jz/hpm);min_JZm_xy[xi*max_x+yi]=max(0.0,(double) jz/hpm);max_JZm_xy[xi*max_x+yi]=max(0.0,(double) jz/hpm);
				} else
				{
					n_ind_xy[xi*max_x+yi]=(n_ind_xy[xi*max_x+yi]+n)/2;
					HPa_xy[xi*max_x+yi]=(HPa_xy[xi*max_x+yi]+hpa)/2;
					HPm_xy[xi*max_x+yi]=(HPm_xy[xi*max_x+yi]+hpm)/2;
					HPo_xy[xi*max_x+yi]=(HPo_xy[xi*max_x+yi]+hpo)/2;
					HE_xy[xi*max_x+yi]=(HE_xy[xi*max_x+yi]+(n-(hpa*2)))/2;
					JZG_xy[xi*max_x+yi]=(JZG_xy[xi*max_x+yi]+jz)/2;
					JZa_xy[xi*max_x+yi]=(JZa_xy[xi*max_x+yi]+max(0.0,(double) jz/hpa))/2;
					JZm_xy[xi*max_x+yi]=(JZm_xy[xi*max_x+yi]+max(0.0,(double) jz/hpm))/2;	

					min_ind_xy[xi*max_x+yi]=min(min_ind_xy[xi*max_x+yi],n);
					min_HPa_xy[xi*max_x+yi]=min(min_HPa_xy[xi*max_x+yi],hpa);
					min_HPm_xy[xi*max_x+yi]=min(hpm,min_HPm_xy[xi*max_x+yi]);
					min_HPo_xy[xi*max_x+yi]=min(min_HPo_xy[xi*max_x+yi],hpo);
					min_HE_xy[xi*max_x+yi]=min(min_HE_xy[xi*max_x+yi],(n-(hpa*2)));
					min_JZG_xy[xi*max_x+yi]=min(min_JZG_xy[xi*max_x+yi],jz);
					min_JZa_xy[xi*max_x+yi]=min(min_JZa_xy[xi*max_x+yi],max(0.0,(double) jz/hpa));
					min_JZm_xy[xi*max_x+yi]=min(min_JZm_xy[xi*max_x+yi],max(0.0,(double) jz/hpm));

					max_ind_xy[xi*max_x+yi]=max(max_ind_xy[xi*max_x+yi],n);
					max_HPa_xy[xi*max_x+yi]=max(max_HPa_xy[xi*max_x+yi],hpa);
					max_HPm_xy[xi*max_x+yi]=max(hpm,max_HPm_xy[xi*max_x+yi]);
					max_HPo_xy[xi*max_x+yi]=max(max_HPo_xy[xi*max_x+yi],hpo);
					max_HE_xy[xi*max_x+yi]=max(max_HE_xy[xi*max_x+yi],(n-(hpa*2)));
					max_JZG_xy[xi*max_x+yi]=max(max_JZG_xy[xi*max_x+yi],jz);
					max_JZa_xy[xi*max_x+yi]=max(max_JZa_xy[xi*max_x+yi],max(0.0,(double) jz/hpa));
					max_JZm_xy[xi*max_x+yi]=max(max_JZm_xy[xi*max_x+yi],max(0.0,(double) jz/hpm));
				}
				n=0;hpa=0;hpo=0;hpm=0;jz=0;
			}
			
		//------------------------------------
		// free dynamic memory
		delete pop;
		pop=0;
	}
	
	//--------------------------------------
	//write to output file
	simfile.open(DD_file.c_str(), ofstream::app);
	for (int t=0;t<(max_days+spinup+1);t++)
	{
		simfile << "maxdensity" <<"\t"<< t-6 <<"\t"<< mean(n_ind,nsims,t) <<"\t"<< sd(n_ind,nsims,t) <<"\t"<<
			mean(HPa,nsims,t) <<"\t"<< sd(HPa,nsims,t) <<"\t"<<
			mean(HPo,nsims,t) <<"\t"<< sd(HPo,nsims,t) <<"\t"<<
			mean(HPm,nsims,t) <<"\t"<< sd(HPm,nsims,t) <<"\t"<<
			mean(HE,nsims,t) <<"\t"<< sd(HE,nsims,t) <<"\t"<<
			mean(JZG,nsims,t) <<"\t"<< sd(JZG,nsims,t) <<"\t"<<
			mean(JZa,nsims,t) <<"\t"<< sd(JZa,nsims,t) <<"\t"<<
			mean(JZm,nsims,t) <<"\t"<< sd(JZm,nsims,t) <<"\t"<<
			max_array(n_ind,nsims,t) <<"\t"<< min_array(n_ind,nsims,t) <<"\t"<<
			max_array(HPa,nsims,t) <<"\t"<< min_array(HPa,nsims,t) <<"\t"<<
			max_array(HPo,nsims,t) <<"\t"<< min_array(HPo,nsims,t) <<"\t"<<
			max_array(HPm,nsims,t) <<"\t"<< min_array(HPm,nsims,t) <<"\t"<<
			max_array(HE,nsims,t) <<"\t"<< min_array(HE,nsims,t) <<"\t"<<
			max_array(JZG,nsims,t) <<"\t"<< min_array(JZG,nsims,t) <<"\t"<<
			max_array(JZa,nsims,t) <<"\t"<< min_array(JZa,nsims,t) <<"\t"<<
			max_array(JZm,nsims,t) <<"\t"<< min_array(JZm,nsims,t) <<"\n";
	}
	simfile.close();
	
	xy_simfile.open(DD_xyfile.c_str(), ofstream::app);
	for (int xi=0;xi<max_x;xi++)
		for (int yi=0;yi<max_y;yi++)
		{
		xy_simfile << "maxdensity" <<"\t"<< xi <<"\t"<< yi <<"\t"<< 
			n_occ_xy[xi*max_x+yi] <<"\t"<<
			n_ind_xy[xi*max_x+yi] <<"\t"<<
			HPa_xy[xi*max_x+yi] <<"\t"<<
			HPo_xy[xi*max_x+yi] <<"\t"<< 
			HPm_xy[xi*max_x+yi] <<"\t"<< 
			HE_xy[xi*max_x+yi] <<"\t"<<
			JZG_xy[xi*max_x+yi] <<"\t"<< 
			JZa_xy[xi*max_x+yi] <<"\t"<<
			JZm_xy[xi*max_x+yi] <<"\t"<< 
			max_ind_xy[xi*max_x+yi] <<"\t"<< min_ind_xy[xi*max_x+yi] <<"\t"<<
			max_HPa_xy[xi*max_x+yi] <<"\t"<< min_HPa_xy[xi*max_x+yi] <<"\t"<<
			max_HPo_xy[xi*max_x+yi] <<"\t"<< min_HPo_xy[xi*max_x+yi] <<"\t"<<
			max_HPm_xy[xi*max_x+yi] <<"\t"<< min_HPm_xy[xi*max_x+yi] <<"\t"<<
			max_HE_xy[xi*max_x+yi] <<"\t"<< min_HE_xy[xi*max_x+yi] <<"\t"<<
			max_JZG_xy[xi*max_x+yi] <<"\t"<< min_JZG_xy[xi*max_x+yi] <<"\t"<<
			max_JZa_xy[xi*max_x+yi] <<"\t"<< min_JZa_xy[xi*max_x+yi] <<"\t"<<
			max_JZm_xy[xi*max_x+yi] <<"\t"<< min_JZm_xy[xi*max_x+yi] <<"\n";
		}
	xy_simfile.close();

	xy_init_simfile.open(DD_xyfile_init.c_str(), ofstream::app);
	for (int xi=0;xi<max_x;xi++)
		for (int yi=0;yi<max_y;yi++)
		{
		xy_init_simfile << "maxdensity" <<"\t"<< xi <<"\t"<< yi <<"\t"<< 
			n_occ_xy_init[xi*max_x+yi] <<"\t"<<
			n_ind_xy_init[xi*max_x+yi] <<"\t"<<
			HPa_xy_init[xi*max_x+yi] <<"\t"<<
			HPo_xy_init[xi*max_x+yi] <<"\t"<< 
			HPm_xy_init[xi*max_x+yi] <<"\t"<< 
			HE_xy_init[xi*max_x+yi] <<"\t"<<
			JZG_xy_init[xi*max_x+yi] <<"\t"<< 
			JZa_xy_init[xi*max_x+yi] <<"\t"<<
			JZm_xy_init[xi*max_x+yi] <<"\t"<< 
			max_ind_xy_init[xi*max_x+yi] <<"\t"<< min_ind_xy_init[xi*max_x+yi] <<"\t"<<
			max_HPa_xy_init[xi*max_x+yi] <<"\t"<< min_HPa_xy_init[xi*max_x+yi] <<"\t"<<
			max_HPo_xy_init[xi*max_x+yi] <<"\t"<< min_HPo_xy_init[xi*max_x+yi] <<"\t"<<
			max_HPm_xy_init[xi*max_x+yi] <<"\t"<< min_HPm_xy_init[xi*max_x+yi] <<"\t"<<
			max_HE_xy_init[xi*max_x+yi] <<"\t"<< min_HE_xy_init[xi*max_x+yi] <<"\t"<<
			max_JZG_xy_init[xi*max_x+yi] <<"\t"<< min_JZG_xy_init[xi*max_x+yi] <<"\t"<<
			max_JZa_xy_init[xi*max_x+yi] <<"\t"<< min_JZa_xy_init[xi*max_x+yi] <<"\t"<<
			max_JZm_xy_init[xi*max_x+yi] <<"\t"<< min_JZm_xy_init[xi*max_x+yi] <<"\n";
		}
	xy_init_simfile.close();


	//-----------------------------------------------------------------------
	//--------------------------------------
	//manipulate density

	if (run_dens)
	{
	int max_dens=max_array(n_ind,nsims,0);

	double quantiles[]={0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1};
	int n_quantiles=9;

	for (int quant=densquant;quant<n_quantiles;quant++)
	{

	int max_ind=(int) ceil(max_dens*quantiles[quant]);
	param = new Parameter;
	param->Set_max_birds(max_ind);

	//--------------------------------------
	//reset vectors for storing results [nsims][ndays] to zero
	for (int s=0;s<nsims;s++)
		for (int z=0;z<(1+spinup+max_days);z++)
		{
			n_ind[s][z]=0;	
			HPa[s][z]=0;
			HPm[s][z]=0;
			HPo[s][z]=0;
			HE[s][z]=0;
			JZG[s][z]=0;
			JZa[s][z]=0;
			JZm[s][z]=0;
		}
		
	for (int sz=0;sz<max_ncell;sz++)
	{
		n_ind_xy[sz]=0;	
		HPa_xy[sz]=0;
		HPm_xy[sz]=0;
		HPo_xy[sz]=0;
		HE_xy[sz]=0;
		JZG_xy[sz]=0;
		JZa_xy[sz]=0;
		JZm_xy[sz]=0;
		n_occ_xy[sz]=0;

		min_ind_xy[sz]=0;	
		min_HPa_xy[sz]=0;
		min_HPm_xy[sz]=0;
		min_HPo_xy[sz]=0;
		min_HE_xy[sz]=0;
		min_JZG_xy[sz]=0;
		min_JZa_xy[sz]=0;
		min_JZm_xy[sz]=0;

		max_ind_xy[sz]=0;	
		max_HPa_xy[sz]=0;
		max_HPm_xy[sz]=0;
		max_HPo_xy[sz]=0;
		max_HE_xy[sz]=0;
		max_JZG_xy[sz]=0;
		max_JZa_xy[sz]=0;
		max_JZm_xy[sz]=0;

		n_ind_xy_init[sz]=0;	
		HPa_xy_init[sz]=0;
		HPm_xy_init[sz]=0;
		HPo_xy_init[sz]=0;
		HE_xy_init[sz]=0;
		JZG_xy_init[sz]=0;
		JZa_xy_init[sz]=0;
		JZm_xy_init[sz]=0;
		n_occ_xy_init[sz]=0;

		min_ind_xy_init[sz]=0;	
		min_HPa_xy_init[sz]=0;
		min_HPm_xy_init[sz]=0;
		min_HPo_xy_init[sz]=0;
		min_HE_xy_init[sz]=0;
		min_JZG_xy_init[sz]=0;
		min_JZa_xy_init[sz]=0;
		min_JZm_xy_init[sz]=0;

		max_ind_xy_init[sz]=0;	
		max_HPa_xy_init[sz]=0;
		max_HPm_xy_init[sz]=0;
		max_HPo_xy_init[sz]=0;
		max_HE_xy_init[sz]=0;
		max_JZG_xy_init[sz]=0;
		max_JZa_xy_init[sz]=0;
		max_JZm_xy_init[sz]=0;
	}
	
	//------------------------------
	// start simulations
	for (int sim=0;sim<nsims;sim++)
	{
		//variables
		pop = new Colony;
		int hpa=0,hpo=0,hpm=0,jz=0;
		list<Stork>::iterator ind;
		
		//--------------------------------------
		//home range search

		pop->InitBirds();

		//count numbers
		for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
		{
			if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
			{
				hpa++;
				if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
				else if ((*ind).GetChicks()==0) hpo++;
			}
		}

		n_ind[sim][0]=(int) pop->Birdlist.size();
		HPa[sim][0]=hpa;
		HPm[sim][0]=hpm;
		HPo[sim][0]=hpo;
		HE[sim][0]=(int) pop->Birdlist.size()-(hpa*2);
		JZG[sim][0]=jz;
		JZa[sim][0]=max(0.0,(double) jz/hpa);
		JZm[sim][0]=max(0.0,(double) jz/hpm);
		hpa=0;hpo=0;hpm=0;jz=0;

		//--------------------------------
		// spinup / burnin

		for (int t=0;t<spinup;t++)
		{
			pop->Forage(t);
			pop->GoodNight();

			//count numbers
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
				{
					hpa++;
					if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
					else if ((*ind).GetChicks()==0) hpo++;
				}
			}

			n_ind[sim][t+1]=(int) pop->Birdlist.size();
			HPa[sim][t+1]=hpa;
			HPm[sim][t+1]=hpm;
			HPo[sim][t+1]=hpo;
			HE[sim][t+1]=(int) pop->Birdlist.size()-(hpa*2);
			JZG[sim][t+1]=jz;
			JZa[sim][t+1]=max(0.0,(double) jz/hpa);
			JZm[sim][t+1]=max(0.0,(double) jz/hpm);
			hpa=0;hpo=0;hpm=0;jz=0;
		}

		//--------------------------------
		// output initial spatial distribution after home range formation
		for (int xi=0;xi<max_x;xi++)
			for (int yi=0;yi<max_y;yi++)
			{
				int n=0;
				list<Stork>::iterator ind;

				//count numbers
				for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
				{
					if (((*ind).GetX()==xi)&&((*ind).GetY()==yi)) 
					{
						n++;

						if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
						{
							hpa++;
							if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
							else if ((*ind).GetChicks()==0) hpo++;
						}
					}
				}

				if (n>0) n_occ_xy_init[xi*max_x+yi]++;
				
				// store results
				if (sim==0)
				{
					n_ind_xy_init[xi*max_x+yi]=n;min_ind_xy_init[xi*max_x+yi]=n;max_ind_xy_init[xi*max_x+yi]=n;
					HPa_xy_init[xi*max_x+yi]=hpa;min_HPa_xy_init[xi*max_x+yi]=hpa;max_HPa_xy_init[xi*max_x+yi]=hpa;
					HPm_xy_init[xi*max_x+yi]=hpm;min_HPm_xy_init[xi*max_x+yi]=hpm;max_HPm_xy_init[xi*max_x+yi]=hpm;
					HPo_xy_init[xi*max_x+yi]=hpo;min_HPo_xy_init[xi*max_x+yi]=hpo;max_HPo_xy_init[xi*max_x+yi]=hpo;
					HE_xy_init[xi*max_x+yi]=n-(hpa*2);min_HE_xy_init[xi*max_x+yi]=n-(hpa*2);max_HE_xy_init[xi*max_x+yi]=n-(hpa*2);
					JZG_xy_init[xi*max_x+yi]=jz;min_JZG_xy_init[xi*max_x+yi]=jz;max_JZG_xy_init[xi*max_x+yi]=jz;
					JZa_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpa);min_JZa_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpa);max_JZa_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpa);
					JZm_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpm);min_JZm_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpm);max_JZm_xy_init[xi*max_x+yi]=max(0.0,(double) jz/hpm);
				} else
				{
					n_ind_xy_init[xi*max_x+yi]=(n_ind_xy_init[xi*max_x+yi]+n)/2;
					HPa_xy_init[xi*max_x+yi]=(HPa_xy_init[xi*max_x+yi]+hpa)/2;
					HPm_xy_init[xi*max_x+yi]=(HPm_xy_init[xi*max_x+yi]+hpm)/2;
					HPo_xy_init[xi*max_x+yi]=(HPo_xy_init[xi*max_x+yi]+hpo)/2;
					HE_xy_init[xi*max_x+yi]=(HE_xy_init[xi*max_x+yi]+(n-(hpa*2)))/2;
					JZG_xy_init[xi*max_x+yi]=(JZG_xy_init[xi*max_x+yi]+jz)/2;
					JZa_xy_init[xi*max_x+yi]=(JZa_xy_init[xi*max_x+yi]+max(0.0,(double) jz/hpa))/2;
					JZm_xy_init[xi*max_x+yi]=(JZm_xy_init[xi*max_x+yi]+max(0.0,(double) jz/hpm))/2;	

					min_ind_xy_init[xi*max_x+yi]=min(min_ind_xy_init[xi*max_x+yi],n);
					min_HPa_xy_init[xi*max_x+yi]=min(min_HPa_xy_init[xi*max_x+yi],hpa);
					min_HPm_xy_init[xi*max_x+yi]=min(hpm,min_HPm_xy_init[xi*max_x+yi]);
					min_HPo_xy_init[xi*max_x+yi]=min(min_HPo_xy_init[xi*max_x+yi],hpo);
					min_HE_xy_init[xi*max_x+yi]=min(min_HE_xy_init[xi*max_x+yi],(n-(hpa*2)));
					min_JZG_xy_init[xi*max_x+yi]=min(min_JZG_xy_init[xi*max_x+yi],jz);
					min_JZa_xy_init[xi*max_x+yi]=min(min_JZa_xy_init[xi*max_x+yi],max(0.0,(double) jz/hpa));
					min_JZm_xy_init[xi*max_x+yi]=min(min_JZm_xy_init[xi*max_x+yi],max(0.0,(double) jz/hpm));

					max_ind_xy_init[xi*max_x+yi]=max(max_ind_xy_init[xi*max_x+yi],n);
					max_HPa_xy_init[xi*max_x+yi]=max(max_HPa_xy_init[xi*max_x+yi],hpa);
					max_HPm_xy_init[xi*max_x+yi]=max(hpm,max_HPm_xy_init[xi*max_x+yi]);
					max_HPo_xy_init[xi*max_x+yi]=max(max_HPo_xy_init[xi*max_x+yi],hpo);
					max_HE_xy_init[xi*max_x+yi]=max(max_HE_xy_init[xi*max_x+yi],(n-(hpa*2)));
					max_JZG_xy_init[xi*max_x+yi]=max(max_JZG_xy_init[xi*max_x+yi],jz);
					max_JZa_xy_init[xi*max_x+yi]=max(max_JZa_xy_init[xi*max_x+yi],max(0.0,(double) jz/hpa));
					max_JZm_xy_init[xi*max_x+yi]=max(max_JZm_xy_init[xi*max_x+yi],max(0.0,(double) jz/hpm));
				}
				n=0;hpa=0;hpo=0;hpm=0;jz=0;
			}

		//--------------------------------
		// initiate breeding and simulate breeding period

		pop->InitBreeding();

		for (int t=0;t<max_days;t++)
		{
			pop->Forage(t);
			pop->GoodNight();

			//count numbers
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
				{
					hpa++;
					if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
					else if ((*ind).GetChicks()==0) hpo++;
				}
			}

			n_ind[sim][t+1+spinup]=(int) pop->Birdlist.size();
			HPa[sim][t+1+spinup]=hpa;
			HPm[sim][t+1+spinup]=hpm;
			HPo[sim][t+1+spinup]=hpo;
			HE[sim][t+1+spinup]=(int) pop->Birdlist.size()-(hpa*2);
			JZG[sim][t+1+spinup]=jz;
			JZa[sim][t+1+spinup]=max(0.0,(double) jz/hpa);
			JZm[sim][t+1+spinup]=max(0.0,(double) jz/hpm);
			hpa=0;hpo=0;hpm=0;jz=0;
		}
		
		//----------------------------------------------------
		// loop over all grid cells to count and store spatial results
		for (int xi=0;xi<max_x;xi++)
			for (int yi=0;yi<max_y;yi++)
			{
				int n=0;
				list<Stork>::iterator ind;

				//count numbers
				for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
				{
					if (((*ind).GetX()==xi)&&((*ind).GetY()==yi)) 
					{
						n++;

						if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
						{
							hpa++;
							if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
							else if ((*ind).GetChicks()==0) hpo++;
						}
					}
				}

				if (n>0) n_occ_xy[xi*max_x+yi]++;
				
				// store results
				if (sim==0)
				{
					n_ind_xy[xi*max_x+yi]=n;min_ind_xy[xi*max_x+yi]=n;max_ind_xy[xi*max_x+yi]=n;
					HPa_xy[xi*max_x+yi]=hpa;min_HPa_xy[xi*max_x+yi]=hpa;max_HPa_xy[xi*max_x+yi]=hpa;
					HPm_xy[xi*max_x+yi]=hpm;min_HPm_xy[xi*max_x+yi]=hpm;max_HPm_xy[xi*max_x+yi]=hpm;
					HPo_xy[xi*max_x+yi]=hpo;min_HPo_xy[xi*max_x+yi]=hpo;max_HPo_xy[xi*max_x+yi]=hpo;
					HE_xy[xi*max_x+yi]=n-(hpa*2);min_HE_xy[xi*max_x+yi]=n-(hpa*2);max_HE_xy[xi*max_x+yi]=n-(hpa*2);
					JZG_xy[xi*max_x+yi]=jz;min_JZG_xy[xi*max_x+yi]=jz;max_JZG_xy[xi*max_x+yi]=jz;
					JZa_xy[xi*max_x+yi]=max(0.0,(double) jz/hpa);min_JZa_xy[xi*max_x+yi]=max(0.0,(double) jz/hpa);max_JZa_xy[xi*max_x+yi]=max(0.0,(double) jz/hpa);
					JZm_xy[xi*max_x+yi]=max(0.0,(double) jz/hpm);min_JZm_xy[xi*max_x+yi]=max(0.0,(double) jz/hpm);max_JZm_xy[xi*max_x+yi]=max(0.0,(double) jz/hpm);
				} else
				{
					n_ind_xy[xi*max_x+yi]=(n_ind_xy[xi*max_x+yi]+n)/2;
					HPa_xy[xi*max_x+yi]=(HPa_xy[xi*max_x+yi]+hpa)/2;
					HPm_xy[xi*max_x+yi]=(HPm_xy[xi*max_x+yi]+hpm)/2;
					HPo_xy[xi*max_x+yi]=(HPo_xy[xi*max_x+yi]+hpo)/2;
					HE_xy[xi*max_x+yi]=(HE_xy[xi*max_x+yi]+(n-(hpa*2)))/2;
					JZG_xy[xi*max_x+yi]=(JZG_xy[xi*max_x+yi]+jz)/2;
					JZa_xy[xi*max_x+yi]=(JZa_xy[xi*max_x+yi]+max(0.0,(double) jz/hpa))/2;
					JZm_xy[xi*max_x+yi]=(JZm_xy[xi*max_x+yi]+max(0.0,(double) jz/hpm))/2;	

					min_ind_xy[xi*max_x+yi]=min(min_ind_xy[xi*max_x+yi],n);
					min_HPa_xy[xi*max_x+yi]=min(min_HPa_xy[xi*max_x+yi],hpa);
					min_HPm_xy[xi*max_x+yi]=min(hpm,min_HPm_xy[xi*max_x+yi]);
					min_HPo_xy[xi*max_x+yi]=min(min_HPo_xy[xi*max_x+yi],hpo);
					min_HE_xy[xi*max_x+yi]=min(min_HE_xy[xi*max_x+yi],(n-(hpa*2)));
					min_JZG_xy[xi*max_x+yi]=min(min_JZG_xy[xi*max_x+yi],jz);
					min_JZa_xy[xi*max_x+yi]=min(min_JZa_xy[xi*max_x+yi],max(0.0,(double) jz/hpa));
					min_JZm_xy[xi*max_x+yi]=min(min_JZm_xy[xi*max_x+yi],max(0.0,(double) jz/hpm));

					max_ind_xy[xi*max_x+yi]=max(max_ind_xy[xi*max_x+yi],n);
					max_HPa_xy[xi*max_x+yi]=max(max_HPa_xy[xi*max_x+yi],hpa);
					max_HPm_xy[xi*max_x+yi]=max(hpm,max_HPm_xy[xi*max_x+yi]);
					max_HPo_xy[xi*max_x+yi]=max(max_HPo_xy[xi*max_x+yi],hpo);
					max_HE_xy[xi*max_x+yi]=max(max_HE_xy[xi*max_x+yi],(n-(hpa*2)));
					max_JZG_xy[xi*max_x+yi]=max(max_JZG_xy[xi*max_x+yi],jz);
					max_JZa_xy[xi*max_x+yi]=max(max_JZa_xy[xi*max_x+yi],max(0.0,(double) jz/hpa));
					max_JZm_xy[xi*max_x+yi]=max(max_JZm_xy[xi*max_x+yi],max(0.0,(double) jz/hpm));
				}
				n=0;hpa=0;hpo=0;hpm=0;jz=0;
			}
			
		//------------------------------------
		// free dynamic memory
		delete pop;
		pop=0;
	}
	
	//--------------------------------------
	//write to output file
	simfile.open(DD_file.c_str(), ofstream::app);

	for (int t=0;t<(max_days+spinup+1);t++)
	{
		simfile << max_ind <<"\t"<< t-6 <<"\t"<< mean(n_ind,nsims,t) <<"\t"<< sd(n_ind,nsims,t) <<"\t"<<
			mean(HPa,nsims,t) <<"\t"<< sd(HPa,nsims,t) <<"\t"<<
			mean(HPo,nsims,t) <<"\t"<< sd(HPo,nsims,t) <<"\t"<<
			mean(HPm,nsims,t) <<"\t"<< sd(HPm,nsims,t) <<"\t"<<
			mean(HE,nsims,t) <<"\t"<< sd(HE,nsims,t) <<"\t"<<
			mean(JZG,nsims,t) <<"\t"<< sd(JZG,nsims,t) <<"\t"<<
			mean(JZa,nsims,t) <<"\t"<< sd(JZa,nsims,t) <<"\t"<<
			mean(JZm,nsims,t) <<"\t"<< sd(JZm,nsims,t) <<"\t"<<
			max_array(n_ind,nsims,t) <<"\t"<< min_array(n_ind,nsims,t) <<"\t"<<
			max_array(HPa,nsims,t) <<"\t"<< min_array(HPa,nsims,t) <<"\t"<<
			max_array(HPo,nsims,t) <<"\t"<< min_array(HPo,nsims,t) <<"\t"<<
			max_array(HPm,nsims,t) <<"\t"<< min_array(HPm,nsims,t) <<"\t"<<
			max_array(HE,nsims,t) <<"\t"<< min_array(HE,nsims,t) <<"\t"<<
			max_array(JZG,nsims,t) <<"\t"<< min_array(JZG,nsims,t) <<"\t"<<
			max_array(JZa,nsims,t) <<"\t"<< min_array(JZa,nsims,t) <<"\t"<<
			max_array(JZm,nsims,t) <<"\t"<< min_array(JZm,nsims,t) <<"\n";
	}

	simfile.close();

	
	xy_simfile.open(DD_xyfile.c_str(), ofstream::app);
	for (int xi=0;xi<max_x;xi++)
		for (int yi=0;yi<max_y;yi++)
		{
		xy_simfile << max_ind <<"\t"<< xi <<"\t"<< yi <<"\t"<< 
			n_occ_xy[xi*max_x+yi] <<"\t"<<
			n_ind_xy[xi*max_x+yi] <<"\t"<<
			HPa_xy[xi*max_x+yi] <<"\t"<<
			HPo_xy[xi*max_x+yi] <<"\t"<< 
			HPm_xy[xi*max_x+yi] <<"\t"<< 
			HE_xy[xi*max_x+yi] <<"\t"<<
			JZG_xy[xi*max_x+yi] <<"\t"<< 
			JZa_xy[xi*max_x+yi] <<"\t"<<
			JZm_xy[xi*max_x+yi] <<"\t"<< 
			max_ind_xy[xi*max_x+yi] <<"\t"<< min_ind_xy[xi*max_x+yi] <<"\t"<<
			max_HPa_xy[xi*max_x+yi] <<"\t"<< min_HPa_xy[xi*max_x+yi] <<"\t"<<
			max_HPo_xy[xi*max_x+yi] <<"\t"<< min_HPo_xy[xi*max_x+yi] <<"\t"<<
			max_HPm_xy[xi*max_x+yi] <<"\t"<< min_HPm_xy[xi*max_x+yi] <<"\t"<<
			max_HE_xy[xi*max_x+yi] <<"\t"<< min_HE_xy[xi*max_x+yi] <<"\t"<<
			max_JZG_xy[xi*max_x+yi] <<"\t"<< min_JZG_xy[xi*max_x+yi] <<"\t"<<
			max_JZa_xy[xi*max_x+yi] <<"\t"<< min_JZa_xy[xi*max_x+yi] <<"\t"<<
			max_JZm_xy[xi*max_x+yi] <<"\t"<< min_JZm_xy[xi*max_x+yi] <<"\n";
		}
	xy_simfile.close();

	xy_init_simfile.open(DD_xyfile_init.c_str(), ofstream::app);
	for (int xi=0;xi<max_x;xi++)
		for (int yi=0;yi<max_y;yi++)
		{
		xy_init_simfile << "maxdensity" <<"\t"<< xi <<"\t"<< yi <<"\t"<< 
			n_occ_xy_init[xi*max_x+yi] <<"\t"<<
			n_ind_xy_init[xi*max_x+yi] <<"\t"<<
			HPa_xy_init[xi*max_x+yi] <<"\t"<<
			HPo_xy_init[xi*max_x+yi] <<"\t"<< 
			HPm_xy_init[xi*max_x+yi] <<"\t"<< 
			HE_xy_init[xi*max_x+yi] <<"\t"<<
			JZG_xy_init[xi*max_x+yi] <<"\t"<< 
			JZa_xy_init[xi*max_x+yi] <<"\t"<<
			JZm_xy_init[xi*max_x+yi] <<"\t"<< 
			max_ind_xy_init[xi*max_x+yi] <<"\t"<< min_ind_xy_init[xi*max_x+yi] <<"\t"<<
			max_HPa_xy_init[xi*max_x+yi] <<"\t"<< min_HPa_xy_init[xi*max_x+yi] <<"\t"<<
			max_HPo_xy_init[xi*max_x+yi] <<"\t"<< min_HPo_xy_init[xi*max_x+yi] <<"\t"<<
			max_HPm_xy_init[xi*max_x+yi] <<"\t"<< min_HPm_xy_init[xi*max_x+yi] <<"\t"<<
			max_HE_xy_init[xi*max_x+yi] <<"\t"<< min_HE_xy_init[xi*max_x+yi] <<"\t"<<
			max_JZG_xy_init[xi*max_x+yi] <<"\t"<< min_JZG_xy_init[xi*max_x+yi] <<"\t"<<
			max_JZa_xy_init[xi*max_x+yi] <<"\t"<< min_JZa_xy_init[xi*max_x+yi] <<"\t"<<
			max_JZm_xy_init[xi*max_x+yi] <<"\t"<< min_JZm_xy_init[xi*max_x+yi] <<"\n";
		}
	xy_init_simfile.close();
	
	//-----------------------------
	delete param;
	param=0;

	} }// end loop over density simulations
}

//---------------------------------------------------------------------------
// run simulations for GPS-tagged storks and different density levels
void Simulation::SimDensity_GPS(void)
{
	int nsims=1;
	int max_ncell=max_x*max_y;
	Colony *pop=0;
	Parameter *param=0;

	std::string ss;
	std::string filename="OUTPUT/SimDDGPS_ptLotteryPS0Memoryw1";
	ss = "_meta.txt";
	std::string meta_file=filename+ss;
	ss = "_GPS.txt";
	std::string GPS_file=filename+ss;
	ss = "_pop.txt";
	std::string pop_file=filename+ss;

	//---------------------------------------
	// prepare vectors for storing results
	vecint n_ind;
	SetVectorSize(n_ind,nsims);

	//--------------------------------------
	//prepare file output

	std::ofstream meta_simfile (meta_file.c_str());
	meta_simfile<< "economy=\t" << economy << "\n"
		<<"forageChoice=\t"<<forageChoice<<"\n"
		<<"DistmaxSetting=\t"<<DistmaxSetting<<"\n"
		<<"GPStag=\t"<<GPStag<<"\n"
		<<"MemoryMap=\t"<<MemoryMap<<"\n"
		<<"memoryweight=\t"<<memoryweight<<"\n"
		<<"pseudocolonizer=\t"<<pseudocolonizer<<"\n";

	meta_simfile.close();
	//-------------

	std::ofstream pop_simfile (pop_file.c_str());
	pop_simfile<< "sim\t" << "max_birds\t"<<"day\t"<<"n_ind\t"<<"HPa\t"<<"HPo\t"<<"HPm\t"<<"HE\t"<<"JZG\t"<<"JZa\t"<<"JZm\n";
	pop_simfile.close();
	//-------------

	std::ofstream GPS_simfile (GPS_file.c_str());
	GPS_simfile<< "sim\t" << "max_birds\t"<<"day\t"<<"ID\t"<<"PartnerID\t"<<"n_chicks\t"<<"AgeChicks\t"<<"NestX\t"<<"NestY\t"
		<< "GPSx\t" << "GPSy\t" << "time\t" << "flighttime\t" << "AvailableEnergy\t" << "ConsumedEnergy\n";
	GPS_simfile.close();

	//--------------------------------------
	//start simulations - full capacity /max. density
	for (int sim=0;sim<nsims;sim++)
	{
		//variables
		pop = new Colony;
		int hpa=0,hpo=0,hpm=0,jz=0;
		list<Stork>::iterator ind;
		
		//--------------------------------------
		//home range search

		pop->InitBirds();

		//count max. numbers
		n_ind[sim]=(int) pop->Birdlist.size();


		//--------------------------------
		// spinup / burnin
		for (int t=0;t<spinup;t++)
		{
			pop->Forage(t-spinup);


			//------------------------------------------------------
			//count numbers - NOTE! numbers correspond to those experienced over this day - breeding status etc. updated after outputting results
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) hpa++;
			}

			//--------------------------------------
			//write to output file
			pop_simfile.open(pop_file.c_str(), ofstream::app);
			// "sim\t" << "max_birds\t"<<"day\t"<< "n_ind\t" << "HPa\t"<<"HPo\t"<<"HPm\t"<<"HE\t"<<"JZG\t"<<"JZa\t"<<"JZm\n";
		
			pop_simfile << sim << "\t" << max_birds <<"\t"<< t-spinup <<"\t"<< pop->Birdlist.size() <<"\t"<< hpa <<"\t"<<"0\t"<<"0\t"<<
					(int) pop->Birdlist.size()-(hpa*2) <<"\t" <<"0\t"<<"0\t"<<"0\n";
			pop_simfile.close();

			hpa=0;

			//--------------------------------------
			//write to output file
			GPS_simfile.open(GPS_file.c_str(), ofstream::app);
			// "sim\t" << "max_birds\t"<<"day\t"<<"ID\t"<<"PartnerID\t"<<"n_chicks\t"<<"AgeChicks\t"<<"NestX\t"<<"NestY\t"
			//<< "GPSx\t" << "GPSy\t" << "time\t" << "flighttime\t" << "AvailableEnergy\t" << "ConsumedEnergy\n";
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				int num=(*ind).GPSdata.GetGPSsize();

				for (int i=0;i<num;i++)
				{
					GPS_simfile << sim <<"\t"<< max_birds <<"\t"<< t-spinup <<"\t"<< (*ind).GetId() <<"\t"<< (*ind).GetMateID() <<"\t"<<
						(*ind).GetChicks() <<"\t"<< (*ind).GetAgeChicks() <<"\t"<< (*ind).GetX() <<"\t"<< (*ind).GetY() <<"\t"<<
						(*ind).GPSdata.GetGPSx()[i] <<"\t"<< (*ind).GPSdata.GetGPSy()[i] <<"\t"<< (*ind).GPSdata.GetGPStime()[i] <<"\t"<<
						(*ind).GPSdata.Getflighttime()[i] <<"\t"<< (*ind).GPSdata.GetAvailableEnergy()[i] <<"\t"<<
						(*ind).GPSdata.GetConsumedEnergy()[i] <<"\n";			
				}
			}
			GPS_simfile.close();

			//---------------------------------------
			// update breeding status, initialise variables
			pop->GoodNight();

		} //end loop over spinup

		//-------------------------------------------------------------
		// initiate breeding and simulate breeding period

		pop->InitBreeding();

		for (int t=0;t<max_days;t++)
		{
			pop->Forage(t);

			hpa=0;hpm=0;hpo=0;jz=0;
			//------------------
			//count numbers: NOTE! numbers correspond to those experienced over this day - breeding status etc. updated after outputting results
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
				{					
					hpa++;
					if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
					else if ((*ind).GetChicks()==0) hpo++;
				}
			}

			//--------------------------------------
			//write to output file
			pop_simfile.open(pop_file.c_str(), ofstream::app);
			// "sim\t" << "max_birds\t"<<"day\t"<< "n_ind\t" << "HPa\t"<<"HPo\t"<<"HPm\t"<<"HE\t"<<"JZG\t"<<"JZa\t"<<"JZm\n";
		
			pop_simfile << sim << "\t" << max_birds <<"\t"<< t <<"\t"<< pop->Birdlist.size() <<"\t"<< hpa <<"\t"<< hpo << "\t"<< hpm <<"\t"<<
					(int) pop->Birdlist.size()-(hpa*2) <<"\t" << jz <<"\t"<< max(0.0,(double) jz/hpa) <<"\t"<< max(0.0,(double) jz/hpm) <<"\n";
			pop_simfile.close();

			hpa=0;hpm=0;hpo=0;jz=0;

			//--------------------------------------
			//write to output file
			GPS_simfile.open(GPS_file.c_str(), ofstream::app);
			// "sim\t" << "max_birds\t"<<"day\t"<<"ID\t"<<"PartnerID\t"<<"n_chicks\t"<<"AgeChicks\t"<<"NestX\t"<<"NestY\t"
			//<< "GPSx\t" << "GPSy\t" << "time\t" << "flighttime\t" << "AvailableEnergy\t" << "ConsumedEnergy\n";
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				int num=(*ind).GPSdata.GetGPSsize();

				for (int i=0;i<num;i++)
				{
					GPS_simfile << sim <<"\t"<< max_birds <<"\t"<< t <<"\t"<< (*ind).GetId() <<"\t"<< (*ind).GetMateID() <<"\t"<<
						(*ind).GetChicks() <<"\t"<< (*ind).GetAgeChicks() <<"\t"<< (*ind).GetX() <<"\t"<< (*ind).GetY() <<"\t"<<
						(*ind).GPSdata.GetGPSx()[i] <<"\t"<< (*ind).GPSdata.GetGPSy()[i] <<"\t"<< (*ind).GPSdata.GetGPStime()[i] <<"\t"<<
						(*ind).GPSdata.Getflighttime()[i] <<"\t"<< (*ind).GPSdata.GetAvailableEnergy()[i] <<"\t"<<
						(*ind).GPSdata.GetConsumedEnergy()[i] <<"\n";			
				}
			}
			GPS_simfile.close();

			//-------------------
			// update breeding status, initialise variables
			pop->GoodNight();

		} //end loop over breeding period

	} // end loop over nsims

	//-----------------------------------------------------------------------
	//--------------------------------------
	//manipulate density

	int max_dens=max_vec(n_ind);

	double quantiles[]={0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1};
	int n_quantiles=9;

	for (int quant=0;quant<n_quantiles;quant++)
	{

	int max_ind=(int) ceil(max_dens*quantiles[quant]);
	param = new Parameter;
	param->Set_max_birds(max_ind);

	//--------------------------------------
	//start simulations - for different density levels
	for (int sim=0;sim<nsims;sim++)
	{
		//variables
		pop = new Colony;
		int hpa=0,hpo=0,hpm=0,jz=0;
		list<Stork>::iterator ind;
		
		//--------------------------------------
		//home range search

		pop->InitBirds();

		//--------------------------------
		// spinup / burnin
		for (int t=0;t<spinup;t++)
		{
			pop->Forage(t-spinup);


			//------------------------------------------------------
			//count numbers - NOTE! numbers correspond to those experienced over this day - breeding status etc. updated after outputting results
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) hpa++;
			}

			//--------------------------------------
			//write to output file
			pop_simfile.open(pop_file.c_str(), ofstream::app);
			// "sim\t" << "max_birds\t"<<"day\t"<< "n_ind\t" << "HPa\t"<<"HPo\t"<<"HPm\t"<<"HE\t"<<"JZG\t"<<"JZa\t"<<"JZm\n";
		
			pop_simfile << sim << "\t" << max_birds <<"\t"<< t-spinup <<"\t"<< pop->Birdlist.size() <<"\t"<< hpa <<"\t"<<"0\t"<<"0\t"<<
					(int) pop->Birdlist.size()-(hpa*2) <<"\t" <<"0\t"<<"0\t"<<"0\n";
			pop_simfile.close();

			hpa=0;

			//--------------------------------------
			//write to output file
			GPS_simfile.open(GPS_file.c_str(), ofstream::app);
			// "sim\t" << "max_birds\t"<<"day\t"<<"ID\t"<<"PartnerID\t"<<"n_chicks\t"<<"AgeChicks\t"<<"NestX\t"<<"NestY\t"
			//<< "GPSx\t" << "GPSy\t" << "time\t" << "flighttime\t" << "AvailableEnergy\t" << "ConsumedEnergy\n";
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				int num=(*ind).GPSdata.GetGPSsize();

				for (int i=0;i<num;i++)
				{
					GPS_simfile << sim <<"\t"<< max_birds <<"\t"<< t-spinup <<"\t"<< (*ind).GetId() <<"\t"<< (*ind).GetMateID() <<"\t"<<
						(*ind).GetChicks() <<"\t"<< (*ind).GetAgeChicks() <<"\t"<< (*ind).GetX() <<"\t"<< (*ind).GetY() <<"\t"<<
						(*ind).GPSdata.GetGPSx()[i] <<"\t"<< (*ind).GPSdata.GetGPSy()[i] <<"\t"<< (*ind).GPSdata.GetGPStime()[i] <<"\t"<<
						(*ind).GPSdata.Getflighttime()[i] <<"\t"<< (*ind).GPSdata.GetAvailableEnergy()[i] <<"\t"<<
						(*ind).GPSdata.GetConsumedEnergy()[i] <<"\n";			
				}
			}
			GPS_simfile.close();

			//---------------------------------------
			// update breeding status, initialise variables
			pop->GoodNight();

		} //end loop over spinup

		//-------------------------------------------------------------
		// initiate breeding and simulate breeding period

		pop->InitBreeding();

		for (int t=0;t<max_days;t++)
		{
			pop->Forage(t);

			hpa=0;hpm=0;hpo=0;jz=0;
			//------------------
			//count numbers: NOTE! numbers correspond to those experienced over this day - breeding status etc. updated after outputting results
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				if (((*ind).IsPaired()==true)&&((*ind).GetId()<(*ind).GetMateID())) 
				{					
					hpa++;
					if ((*ind).GetChicks()>0) {hpm++; jz+=(*ind).GetChicks();}
					else if ((*ind).GetChicks()==0) hpo++;
				}
			}

			//--------------------------------------
			//write to output file
			pop_simfile.open(pop_file.c_str(), ofstream::app);
			// "sim\t" << "max_birds\t"<<"day\t"<< "n_ind\t" << "HPa\t"<<"HPo\t"<<"HPm\t"<<"HE\t"<<"JZG\t"<<"JZa\t"<<"JZm\n";
		
			pop_simfile << sim << "\t" << max_birds <<"\t"<< t <<"\t"<< pop->Birdlist.size() <<"\t"<< hpa <<"\t"<< hpo << "\t"<< hpm <<"\t"<<
					(int) pop->Birdlist.size()-(hpa*2) <<"\t" << jz <<"\t"<< max(0.0,(double) jz/hpa) <<"\t"<< max(0.0,(double) jz/hpm) <<"\n";
			pop_simfile.close();

			hpa=0;hpm=0;hpo=0;jz=0;

			//--------------------------------------
			//write to output file
			GPS_simfile.open(GPS_file.c_str(), ofstream::app);
			// "sim\t" << "max_birds\t"<<"day\t"<<"ID\t"<<"PartnerID\t"<<"n_chicks\t"<<"AgeChicks\t"<<"NestX\t"<<"NestY\t"
			//<< "GPSx\t" << "GPSy\t" << "time\t" << "flighttime\t" << "AvailableEnergy\t" << "ConsumedEnergy\n";
			for(ind=pop->Birdlist.begin(); ind != pop->Birdlist.end(); ++ind)
			{
				int num=(*ind).GPSdata.GetGPSsize();

				for (int i=0;i<num;i++)
				{
					GPS_simfile << sim <<"\t"<< max_birds <<"\t"<< t <<"\t"<< (*ind).GetId() <<"\t"<< (*ind).GetMateID() <<"\t"<<
						(*ind).GetChicks() <<"\t"<< (*ind).GetAgeChicks() <<"\t"<< (*ind).GetX() <<"\t"<< (*ind).GetY() <<"\t"<<
						(*ind).GPSdata.GetGPSx()[i] <<"\t"<< (*ind).GPSdata.GetGPSy()[i] <<"\t"<< (*ind).GPSdata.GetGPStime()[i] <<"\t"<<
						(*ind).GPSdata.Getflighttime()[i] <<"\t"<< (*ind).GPSdata.GetAvailableEnergy()[i] <<"\t"<<
						(*ind).GPSdata.GetConsumedEnergy()[i] <<"\n";			
				}
			}
			GPS_simfile.close();

			//-------------------
			// update breeding status, initialise variables
			pop->GoodNight();

		} //end loop over breeding period

	} // end loop over nsims

	//-----------------------------
	delete param;
	param=0;

	}
}
