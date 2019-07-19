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
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "func.h"

using namespace std;
 
int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");

	if (Primes.is_open()) {
		Primes >> p1 >> p2 ;
	}
	else cerr << "PROBLEM: Unable to open Primes" << endl;
	
	Primes.close();

	ifstream input("seed.in");
	string property;

	if (input.is_open()) {
		while ( !input.eof() ) {
			input >> property;
			if( property == "RANDOMSEED" ) {
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	}
	else cerr << "PROBLEM: Unable to open seed.in" << endl;

	//for(int i=0; i<20; i++) {
	//	cout << rnd.Rannyu() << endl;
	//}

	double a = 1.;
	double va;
	int rep = 10000;
	int nstep = 100;
	double x, y, z;
	vector<double> r(rep);
	ofstream risd, risc;

	risd.open("Risultati/RW.discr");
	risc.open("Risultati/RW.cont");

// 1. discrete RW
	cout << "Discrete RW starts..." << endl;
	cout << "Repetitions = " << rep << endl;
	cout << "Steps from 1 to " << nstep << endl << endl;
	for(int k=0; k<=nstep; k++) {
		for(int i=0; i<rep; i++) {
			x = 0.;
			y = 0.;
			z = 0.;
// Make nsteps: choose the direction and move
			for(int j=0; j<k; j++) {
				va = rnd.Rannyu(-3., 3.);

				if(va <= -2.)
					x--;

				else if(va >= -2. && va <= -1.)
					y--;

				else if(va >= -1. && va <= 0.)
					z--;

				else if(va >= 0. && va <= 1.)
					x++;

				else if(va >= 1. && va <= 2.)
					y++;

				else
					z++;
			}
// Compute the path length squared
			r[i] = pow(x, 2)
				+ pow(y, 2) + pow(z, 2);
		}

// Average and error
		risd << k << "    " << sqrt(Mean(r)) << "    " << sqrt(Sigma(r)) << endl;
	}
	
	risd.close();

// 2. Continuum RW		
	cout << "Continuum RW starts..." << endl;
	cout << "Repetitions = " << rep << endl;
	cout << "Steps from 1 to " << nstep << endl;
	double th, phi;

	for(int k=0; k<=nstep; k++) {
		for(int i=0; i<rep; i++) {
			x = 0.;
			y = 0.;
			z = 0.;
// Make nsteps in spherical coordinates
			for(int j=0; j<k; j++) {
				th = rnd.Rannyu(0., M_PI);
				phi = rnd.Rannyu(0., 2.*M_PI);

				x += a * sin(th) * cos(phi);
				y += a * sin(th) * sin(phi);
				z += a * cos(th);
			}
// Compute the path length
			r[i] = pow(x, 2)
				+ pow(y, 2) + pow(z, 2);
		}

// Average and error
		risc << k << "    " << sqrt(Mean(r)) << "    " << sqrt(Sigma(r)) << endl;
	}
	
	risc.close();

	rnd.SaveSeed();

	return 0;
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
