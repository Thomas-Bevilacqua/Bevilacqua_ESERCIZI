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
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include "genetic.h"
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

// Construct 30 cities on a circle with r=1 without repetitions
// and check
	cout << "---SALESMAN PROBLEM---" << endl;
	cout << "Performing a simulated annealing algorithm" << endl << endl;

	int ncity = 30;		// Number of cities
	int nmut = 14;	// Number of mutations
	int r = 1.;
	vector<double> x(ncity), y(ncity);
	int i = 1;	
	
	cout << "Number of cities = " << ncity << endl;
	cout << "Number of mutations (probability 10%) = " << nmut << endl;
	cout << "Cities on a circle of radius = " << r << endl << endl;

	x[0] = rnd.Rannyu(-1., 1);
	while(i < ncity) {
		double city = rnd.Rannyu(-1., 1);
		double dep;

		for(int j=0; j<i; j++) {
			if(x[j] != city)
				dep = city;

			else {
				dep = 0;
				break;
			}
		}

		if(dep != 0) {
			x[i] = city;
			i++;
		}
	}
	
	Circle(x, y, r, rnd);

	cout << "Vector x coordinates" << endl;
	Check(x);
	cout << "Vector y coordinates" << endl;
	Check(y);
	Print(x, y, "Risultati/circ.0");
	cout << "Initial configuration printed" << endl << endl;
	//cout << L2(x, y) << endl;

// Find the minimum
	double bi = 0., bf = 70., db = 0.0005;
	int thr, n = 1000;
	vector<double> depx(ncity), depy(ncity);
	ofstream bp;

	bp.open("Risultati/circ.shortest");		// File with best path
	
	cout << "Beta in range = [ " << bi << ":" << bf << " ] with step = " << db << endl;
	cout << "Number of MC steps = " << n << endl << endl;
	
	for(double beta=bi; beta<bf; beta=beta+db) {
		thr = 0;
		for(int i=0; i<n; i++) {
			depx = x;
			depy = y;
			// Find a new configuration
			for(int j=0; j<nmut; j++) {
				double prob = rnd.Rannyu();
				Make_Mutation(depx, depy, rnd, prob);
			}

			double prob = rnd.Rannyu();
			// Substitute if better
			SA(x, y, depx, depy, prob, beta, thr);
		}

		bp << beta << "    " << L2(x, y) << endl;

		if((int)beta % 5 == 0) {
			cout << "Beta = " << beta << "/" << bf << endl;
			cout << "Acceptance ratio = " << (double)thr/(double)n * 100. << endl;
		}
	}
		
	Print(x, y, "Risultati/circ.fin");
	cout << "Optimized configuration printed" << endl << endl;
	bp.close();

// Repeat with a SQUARE
	rnd.SetRandom(seed,p1,p2);
	double l = 1.;
	cout << "Performing the same algorithm with the same input" << endl;
	cout << "Cities in a square with edge = " << 2.*l << endl << endl;
	
	// Generate random cities in a square and check
	for(int i=0; i<ncity; i++) {
		x[i] = rnd.Rannyu(-l, l);
		y[i] = rnd.Rannyu(-l, l);
	}
	
	Check(x, y);
	Print(x, y, "Risultati/sq.0");
	cout << "Initial configuration printed" << endl << endl;

	bp.open("Risultati/sq.shortest");		// File with best path

	for(double beta=bi; beta<bf; beta=beta+db) {
		thr = 0;
		for(int i=0; i<n; i++) {
			depx = x;
			depy = y;
			// Find a new configuration
			for(int j=0; j<nmut; j++) {
				double prob = rnd.Rannyu();
				Make_Mutation(depx, depy, rnd, prob);
			}

			double prob = rnd.Rannyu();
			// Substitute if better
			SA(x, y, depx, depy, prob, beta, thr);
		}

		bp << beta << "    " << L2(x, y) << endl;

		if((int)beta % 5 == 0) {
			cout << "Beta = " << beta << "/" << bf << endl;
			cout << "Acceptance ratio = " << (double)thr/(double)n * 100. << endl;
		}
	}
		
	Print(x, y, "Risultati/sq.fin");
	cout << "Optimized configuration printed" << endl << endl;
	bp.close();
	

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
