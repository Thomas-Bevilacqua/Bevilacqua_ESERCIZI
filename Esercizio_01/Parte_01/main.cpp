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

	//for(int i=0; i<20; i++)
	//	cout << rnd.Rannyu() << endl;

	int throws = 100000;
	int blocks = 100;
	int rapp = throws/blocks;
	vector<double> nrs;

// 1. Evaluate the first integral
	for(int i=0; i<throws; i++)
		nrs.push_back( rnd.Rannyu() );

	block_unc(nrs, blocks, "Risultati/r.ave");


// 2. Evaluate the second integral
	for(int i=0; i<blocks; i++)
		for(int j=0; j<rapp; j++)
			nrs[j + i*rapp] = pow( nrs[j + i*rapp] - 0.5, 2 );

	//cout << nrs.size() << endl;
	//for(double elem : nrs)
	//	cout << elem << endl;

	block_unc(nrs, blocks, "Risultati/r.sigma");


// 3. Chi test
	int M = 100;	// Subintervals
	int n = 10000;
	double exp = (double)n/(double)M;	// Expected
	int dim = 1E6;	// Total throws

	double step = 1./(double)M;
	vector<double> num;
	int count;
	double chi;
	ofstream chi2;

	chi2.open("Risultati/Chi2.test");

// Generate dim randoms
	for(int i=0; i<dim; i++)
		num.push_back( rnd.Rannyu() );

	//cout << num.size() << endl;

	for(int k=0; k<M; k++) {	// Take the n-th subinterval
		chi = 0.;
		for(int i=0; i<M; i++) {	// In the i-th step...
			count = 0;
			for(int j=0; j<n; j++)	// ... count how many fall in
				if( num[j + k*n] >= step*(double)i && num[j + k*n] < step*(double)(i+1) )
					count++;
			
			chi += pow( (double)count - exp, 2)/exp;
		}
		
		chi2 << k+1 << "    " << chi << endl;
	}

	chi2.close();

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
