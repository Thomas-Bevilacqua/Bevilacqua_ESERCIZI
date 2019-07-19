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
#include "mpi.h"

using namespace std;
 
int main (int argc, char *argv[]){

// Parallel computing
	MPI::Init(argc,argv);					// initialization
	int size = MPI::COMM_WORLD.Get_size();	// how many processes
	int rank = MPI::COMM_WORLD.Get_rank();	// who am I	

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");

	if (Primes.is_open()) {
		Primes >> p1 >> p2 ;
	}
	else cerr << "PROBLEM: Unable to open Primes" << endl;
	
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
	int ncity = 30;		// Number of cities
	int nmut = 14;	// Number of mutations
	int r = 1.;
	vector<double> x(ncity), y(ncity);
	int i = 1;

// Prepare the vectors for the results
	double x_irecv[size][ncity];
	double y_irecv[size][ncity];
	vector<double> best_irecv(size);

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

	//cout << "Vector x coordinates" << endl;
	//Check(x);
	//cout << "Vector y coordinates" << endl;
	//Check(y);
	//Print(x, y, "circ.0");
	//cout << "Initial configuration printed" << endl << endl;
	//cout << L2(x, y) << endl;

// Find the minimum
// Set different random for each node
	int pa[size], pb[size];
	for(int i=0; i<size; i++)
		Primes >> pa[i] >> pb[i];

	rnd.SetRandom(seed, pa[rank], pb[rank]);

	double bi = 0., bf = 70., db = 0.005;
	int thr, n = 1000;
	vector<double> depx(ncity), depy(ncity);
		
	//cout << "Beta in range = [ " << bi << ":" << bf << " ] with step = " << db << endl;
	//cout << "Number of MC steps = " << n << endl << endl;
	
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
	}
	
// Prepare the results for each node
	double l_isend = L2(x, y);
	double * x_isend = new double[ncity];
	double * y_isend = new double[ncity];

	for(int i=0; i<ncity; i++) {
		x_isend[i] = x[i];
		y_isend[i] = y[i];
	}
		
// Back to node 0
	MPI_Gather(x_isend, ncity, MPI_DOUBLE, x_irecv[rank], ncity, MPI_DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(y_isend, ncity, MPI_DOUBLE, y_irecv[rank], ncity, MPI_DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(&l_isend, 1, MPI_DOUBLE, &best_irecv[rank], 1, MPI_DOUBLE, 0, MPI::COMM_WORLD);
	
// Print results: best path and correspondent configuration
	if(rank==0) {
		//for(int i=0; i<4; i++)
		//	cout << best_irecv[i] << endl;
		double min = 10000.;
		int ib = 0;
		ofstream ris;

		ris.open("Risultati/Circle.ris");
		
		for(int i=0; i<size; i++)
			if( best_irecv[i] < min ) {
				min = best_irecv[i];
				ib = i;
			}

		//ris << min << endl;
		for(int i=0; i<ncity; i++)
			ris << x_irecv[ib][i] << "    " << y_irecv[ib][i] << endl;
		ris << x_irecv[ib][0] << "    " << y_irecv[ib][0] << endl;

		ris.close();
		cout << "Circle best path = " << min << endl;
		cout << "Circle: results are printed" << endl;
	}

	rnd.SaveSeed();

	MPI::Finalize();	// finish

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
