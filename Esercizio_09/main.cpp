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
	cout << "Performing a genetic algorithm" << endl << endl;

	int ncity = 30;		// Number of cities
	int nmut = 14;	// Number of mutations
	int ngen = 2000;	// Number of generations
	int r = 1.;
	vector<double> x(ncity), y(ncity);
	int i = 1;	
	
	cout << "Number of cities = " << ncity << endl;
	cout << "Number of mutations (probability 10%) = " << nmut << endl;
	cout << "Number of generations (crossover prob. 50%) = " << ngen << endl;
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

// Starting from this vector make a population using
// using some mutations
	int pop = 1000;
	vector< vector<double> > xpop(pop), ypop(pop);

	for(int i=0; i<pop; i++) {
		for(int j=0; j<1000; j++) {
			double prob = rnd.Rannyu();
			Make_Mutation(x, y, rnd, prob);
		}

		xpop[i] = x;
		ypop[i] = y;
		//Print(xpop[i]);
		//Print(ypop[i]);
	}
	cout << "Population of " << pop << " people created" << endl << endl;

// Make new generations
	ofstream bp, al;
	bp.open("Risultati/circ.shortest");		// File with best path for generation
	al.open("Risultati/circ.ave");		// File with average length
			
	for(int i=0; i<ngen; i++) {
		if( i%500 == 0)
			cout << "Generation " << i << "/" << ngen << endl;

		Sons(xpop, ypop, rnd, nmut);
		bp << i+1 << "    " << Shortest(xpop, ypop) << endl;
		al << i+1 << "    " << (Stat(xpop, ypop))[0] << "    " << (Stat(xpop, ypop))[1] << endl;
	}

	bp.close();
	al.close();
	cout << endl;

	Print_best_path(xpop, ypop, "Risultati/circ.fin");
	cout << "Optimized configuration printed" << endl << endl << endl;

// SQUARE
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

	// Generate the population
	for(int i=0; i<pop; i++) {
		for(int j=0; j<1000; j++) {
			double prob = rnd.Rannyu();
			Make_Mutation(x, y, rnd, prob);
		}

		xpop[i] = x;
		ypop[i] = y;
		//Print(xpop[i]);
		//Print(ypop[i]);
	}
	cout << "Population of " << pop << " people created" << endl << endl;
	
	// Make new generations
	bp.open("Risultati/sq.shortest");		// File with best path for generation
	al.open("Risultati/sq.ave");		// File with average length
			
	for(int i=0; i<ngen; i++) {
		if( i%500 == 0)
			cout << "Generation " << i << "/" << ngen << endl;

		Sons(xpop, ypop, rnd, nmut);
		bp << i+1 << "    " << Shortest(xpop, ypop) << endl;
		al << i+1 << "    " << (Stat(xpop, ypop))[0] << "    " << (Stat(xpop, ypop))[1] << endl;
	}

	bp.close();
	al.close();
	cout << endl;

	Print_best_path(xpop, ypop, "Risultati/sq.fin");
	cout << "Optimized configuration printed" << endl << endl << endl;

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
