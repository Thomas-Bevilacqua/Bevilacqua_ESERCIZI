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
#include <vector>
#include <cmath>
#include <iomanip>
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

	double L = 1.;
	double d = 1.5;
	double x, y, l;
	int count, n;

	int thr = 10000;
	int bl = 100;
	vector<double> pi;
	
// Make thr throws evaluations of pi
	for(int i=0; i<thr; i++) {
		count = 0;
		n = 0;
		while(n <= thr) {					// Make thr throws ..
			l = rnd.Rannyu(0., d);			// Distance center-line
			x = rnd.Rannyu(-1., 1.);
			y = rnd.Rannyu();				

			if(x*x + y*y < 1.) {			// Compute the sin
				if(L*Sin(x,y)/2. >= d-l || L*Sin(x,y)/2. >= l )	// If cross a line ..
					count++;		
				
				n++;		
			}
		}
		
		pi.push_back( 2.*L*(double)thr/( d*(double)count) );
	}
	
	//cout << pi.size() << endl;
	//for(double elem : pi)
		//cout << elem << endl;

// Print the results
	block_unc(pi, bl, "Risultati/Pi.eval");

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
