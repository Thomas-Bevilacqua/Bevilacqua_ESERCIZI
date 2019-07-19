/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() { 
	Input(); //Inizialization
	
	for(double t=temp; t<=temp_fin; t=t+dT) {
		Load_spin();	//Reset spin configuration

		cout << "Equilibration..." << endl;
		for(int i=0; i<neq; i++)	// Equilibrate
			Move(metro);
	
		cout << "Simulation at temperature = " << t << endl;
		beta = 1./t;

		for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation
			Reset(iblk);   //Reset block averages
			
			for(int istep=1; istep <= nstep; ++istep) {
				Move(metro);
				Measure();
				Accumulate(); //Update block averages
			}
		
			Averages(iblk);   //Print results for current block
		}

		ConfFinal(); //Write final configuration
		
		Temp_print(t);
	}

	return 0;
}


void Input(void) {
	ifstream ReadInput;

	cout << "Classic 1D Ising model             " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Nearest neighbour interaction      " << endl << endl;
	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
	cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();

//Read input informations
	ReadInput.open("input.dat");

	ReadInput >> temp;
	beta = 1.0/temp;
	cout << "Initial temperature = " << temp << endl;

	ReadInput >> temp_fin;
	ReadInput >> dT;
		if(temp_fin == temp)
			cout << "Fixed temperature simulation" << endl;
		
		else {
			cout << "Final temperature = " << temp_fin << endl;
			cout << "Temperature step = " << dT << endl;
		}
	
	ReadInput >> neq;
	cout << "Equilibration steps = " << neq << endl;
	ReadInput >> nspin;
	cout << "Number of spins = " << nspin << endl;

	ReadInput >> J;
	cout << "Exchange interaction = " << J << endl;

	ReadInput >> h;
	cout << "External field = " << h << endl << endl;

	ReadInput >> metro; // if=1 Metropolis else Gibbs

	ReadInput >> nblk;

	ReadInput >> nstep;

	if(metro==1) 
		cout << "The program perform Metropolis moves" << endl;
	else 
		cout << "The program perform Gibbs moves" << endl;

	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;
	ReadInput.close();

//Prepare arrays for measurements
	iu = 0; //Energy
	ic = 1; //Heat capacity
	im = 2; //Magnetization
	ix = 3; //Magnetic susceptibility
	iu2 = 4; //Quadratic energy

	n_props = 5; //Number of observables

//initial configuration
//then save it
	for (int i=0; i<nspin; ++i) {
		if(rnd.Rannyu() >= 0.5) 
			s[i] = 1;
		else 
			s[i] = -1;
	}

	Save_spin();

//Evaluate energy etc. of the initial configuration
	Measure();

//Print initial values for the potential energy and virial and other observables
	cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
	cout << "Initial heat capacity = " << walker[ic] << endl;
	cout << "Initial magnetization = " << walker[im] << endl;
	cout << "Initial susceptibility = " << walker[ix] << endl << endl;
}


void Move(int metro) {
	int o;
	double p;
	double dE, r;

	for(int i=0; i<nspin; ++i) {
		//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
		o = (int)(rnd.Rannyu()*nspin);
		
		if(metro==1) { //Metropolis
			//Evaluate the difference in energy when flipping s[o]
			dE = 2. * J * s[o] * (s[Pbc(o-1)] + s[Pbc(o+1)]) + 2. * h * s[o];
			p = exp(-beta * dE);
			
			//cout << dE << "	" << p << endl;

			//Flip when satisy the Metropolis conditions
			if( dE < 0. ) {
				s[o] = -s[o];
				accepted++;
			}

			else {
				r = rnd.Rannyu();
				if( r <= p ) {
					s[o] = -s[o];
					accepted++;
				}
			}
			
			attempted++;
		}

		else { //Gibbs sampling
			//Assume the flip of a spin up (s[o]=+ -)
			dE = 2. * J * beta * (s[Pbc(o-1)] + s[Pbc(o+1)]) + 2. * h * beta;
			p = 1./(1. + exp( -dE ) );	
			r = rnd.Rannyu();

			if( r < p ) {
				s[o] = 1.;
				accepted++;
			}

			else
				s[o] = -1.;
			
			attempted++;
		}
	}
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure() {
	double u = 0.0, m = 0.0;

//cycle over spins
	for (int i=0; i<nspin; ++i) {
		u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		m += s[i];
	}

	walker[iu] = u;
	walker[iu2] = u * u;	//For heat capacity
	walker[ix] = beta * m * m/(double)nspin;
	walker[im] = m/(double)nspin;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) { //Print results for current block
    
	ofstream Ene, Heat, Mag, Chi;
	//const int wd=12;
	
	if(nblk%iblk == 0) {
		cout << "Block number " << iblk << endl;
		cout << "Acceptance rate " << accepted/attempted * 100. << " %" << endl << endl;
		cout << "----------------------------" << endl << endl;
	}
	
	if(metro == 1) {
		Ene.open("Risultati/ene_ave.met",ios::app);
		Heat.open("Risultati/heat_ave.met",ios::app);
		Chi.open("Risultati/chi_ave.met",ios::app);
		Mag.open("Risultati/mag_ave.met",ios::app);
	}

	else {
		Ene.open("Risultati/ene_ave.gib",ios::app);
		Heat.open("Risultati/heat_ave.gib",ios::app);
		Chi.open("Risultati/chi_ave.gib",ios::app);
		Mag.open("Risultati/mag_ave.gib",ios::app);
	}

	stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
	glob_av[iu]  += stima_u;
	double uu = glob_av[iu] * (double)nspin/(double)iblk;

	stima_u2 = blk_av[iu2]/blk_norm/(double)nspin;
	glob_av[iu2]  += stima_u2;
	double uu2 = glob_av[iu2] * (double)nspin/(double)iblk;

	glob_av2[iu] += stima_u*stima_u;
	err_u=Error(glob_av[iu],glob_av2[iu],iblk);

	Ene << iblk << "   " << stima_u << "   " << glob_av[iu]/(double)iblk << "   " << err_u << endl;

	stima_c = pow(beta, 2) * (uu2 - pow(uu, 2) )/(double)nspin; //Heat capacity
	glob_av[ic]  += stima_c;
	glob_av2[ic] += stima_c*stima_c;
	err_c=Error(glob_av[ic],glob_av2[ic],iblk);

	Heat << iblk << "   " << stima_c << "   " << glob_av[ic]/(double)iblk << "   " << err_c << endl;
	
	stima_x = blk_av[ix]/blk_norm; //Susceptibility
	glob_av[ix]  += stima_x;
	glob_av2[ix] += stima_x*stima_x;
	err_x=Error(glob_av[ix],glob_av2[ix],iblk);

	Chi << iblk << "   " << stima_x << "   " << glob_av[ix]/(double)iblk << "   " << err_x << endl;

	stima_m = blk_av[im]/blk_norm; //Magnetization
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m=Error(glob_av[im],glob_av2[im],iblk);

	Mag << iblk << "   " << stima_m << "   " << glob_av[im]/(double)iblk << "   " << err_m << endl;

	Ene.close();
	Heat.close();
	Chi.close();
	Mag.close();
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Save_spin() {
	for(int i=0; i<nspin; i++)
		deps[i] = s[i];
}

void Load_spin() {
	for(int i=0; i<nspin; i++)
		s[i] = deps[i];
}

void Temp_print(double temp) {
	ofstream Ene, Heat, Mag, Chi;

	if(metro == 1) {
		Ene.open("Risultati/ene.met",ios::app);
		Heat.open("Risultati/heat.met",ios::app);
		Chi.open("Risultati/chi.met",ios::app);
		Mag.open("Risultati/mag.met",ios::app);
	}

	else {
		Ene.open("Risultati/ene.gib",ios::app);
		Heat.open("Risultati/heat.gib",ios::app);
		Chi.open("Risultati/chi.gib",ios::app);
		Mag.open("Risultati/mag.gib",ios::app);
	}

	Ene << temp << "   " << glob_av[iu]/(double)nblk << "   " << err_u << endl;
	Heat << temp << "   " << glob_av[ic]/(double)nblk << "   " << err_c << endl;
	Chi << temp << "   " << glob_av[ix]/(double)nblk << "   " << err_x << endl;
	Mag << temp << "   " << glob_av[im]/(double)nblk << "   " << err_m << endl;

	Ene.close();
	Heat.close();
	Chi.close();
	Mag.close();
}
	

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
