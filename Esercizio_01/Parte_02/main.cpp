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
#include "random.h"

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

	int n = 10000;
	int N[4] = {1, 2, 10, 100};
	double Ss, Se, Sl;
	ofstream stan, ex, lor;
	
// Open the files for the results and compute
	for(int c=0; c<4; c++) {
		stan.open("Risultati/stand_" + to_string(N[c]) + ".dice");
		ex.open("Risultati/exp_" + to_string(N[c]) + ".dice");
		lor.open("Risultati/lorentz_" + to_string(N[c]) + ".dice");

		for(int j=0; j<n; j++) {
			Ss = 0.;
			Se = 0.;
			Sl = 0.;
			for(int k=0; k<N[c]; k++) {
				Ss += (int)rnd.Rannyu(1., 7.); 			// Normal dice
				Se += rnd.Exp(1.);						// Exp dice
				Sl += rnd.Lor(1., 0.);					// Lorentz dice
			}
			
			stan << (double)Ss/(double)N[c] << endl;
			ex << (double)Se/(double)N[c] << endl;
			lor << (double)Sl/(double)N[c] << endl;
		}
	
		stan.close();
		ex.close();
		lor.close();
	}
	
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
