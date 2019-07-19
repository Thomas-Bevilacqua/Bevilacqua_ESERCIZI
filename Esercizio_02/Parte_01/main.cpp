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

	int thr = 100000;
	int bl = 100;
	vector<double> f(thr);
	ofstream of;

	of.open("Risultati/sampled.points");

// 1. Sampling a uniform distribution
	for(int i=0; i<thr; i++)
		f[i] = M_PI/2. * cos(M_PI * rnd.Rannyu()/2.);

	block_unc(f, bl, "Risultati/int.uniform");

// 2. Importance sampling with p = 2(1-x)
	for(int i=0; i<thr; i++) {
		double x = rnd.Cos();
		of << x << endl;
		f[i] = M_PI * cos(M_PI*x/2.) / (4.*(1.-x));
	}
		
	block_unc(f, bl, "Risultati/int.imp");
	of.close();

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
