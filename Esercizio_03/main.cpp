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

// Input parameters
	double t = 0.;
	double S0 = 100.;
	double s;
	double T = 1.;
	double K = 100.;
	double r = 0.1;
	double sigma = 0.25;

	int M = 100000;
	int bl = 100;
	vector<double> C, P, S;

// 1. Direct sampling
	for(int i=0; i<M; i++) {
		s = S0 * exp( (r - 0.5 * pow(sigma, 2)) * T + sigma * rnd.Gauss(0., T) );
		C.push_back( exp(-r*T) * max(s-K, 0.) );
		P.push_back( exp(-r*T) * max(K-s, 0.) );
		//cout << s << endl;
	}

	block_unc(C, bl, "Risultati/Call.dir");
	block_unc(P, bl, "Risultati/Put.dir");

// 2. Discrete sampling (each interval is long dt)
	int step = 100;
	double dt = (T-t)/(double)step;

	for(int i=0; i<M; i++) {
		s = S0;
		for(t=0.; t<=T; t=t+dt)
			s = s * exp( (r - 0.5 * pow(sigma, 2)) * dt + sigma * rnd.Gauss(0., 1.) * sqrt(dt) );
	
		//cout << t << endl;
		C[i] = exp(-r*T) * max(s-K, 0.);
		P[i] = exp(-r*T) * max(K-s, 0.);
		//cout << s << endl;
	}

	block_unc(C, bl, "Risultati/Call.dis");
	block_unc(P, bl, "Risultati/Put.dis");
	
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
