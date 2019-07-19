#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "random.h"

using namespace std;

// Statisctical uncertainties: Blocking method
void block_unc(vector<double>& data, int blocks, string file) {
	int throws = data.size();
	double sum, sum2;
	int ratio = throws/blocks;
	vector<double> mean, mean2;

	ofstream outfile;
	outfile.open(file);

	if( !outfile.is_open() ) {
		cerr << "Error: can't open " << file << " ..." << endl;
		return;
	}

// Average in each block
	for(int i=0; i<blocks; i++) {
		sum = 0.;
		for(int j=0; j<ratio; j++)
			sum += data[j + i*ratio];

		mean.push_back(sum/(double)ratio);
		mean2.push_back( pow(mean[i], 2) );
	}

	//cout << mean.size() << endl;
	//for(double elem : mean)
	//	cout << elem << endl;

// Global average and print
	for(int i=0; i<blocks; i++) {
		sum = 0.;
		sum2 = 0.;
		for(int j=0; j<i+1; j++) {
			sum += mean[j];
			sum2 += mean2[j];
		}

		outfile << (i+1)*ratio << "    " << sum/(double)(i+1) << "    "
		<< sqrt(sum2/(double)(i+1) - pow(sum/(double)(i+1), 2) )/sqrt( (double)(i+1) ) 
		<< endl;

	}

	outfile.close();
}

// Wave functions (probability distributions)
// Lengths in a0 units
double psi_100(vector<double>& coord) {
	//double a0 = 5.292E-11;
	double r = 0.;

	for(double elem : coord)
		r += elem*elem;
	
	r = sqrt(r);
	//cout << r << endl;

	return pow( exp(-r)/sqrt(M_PI), 2 );
}

double psi_210(vector<double>& coord) {
	//double a0 = 5.292E-11;
	double r = 0.;

	for(double elem : coord)
		r += elem*elem;
	
	r = sqrt(r);
	double cos_th = coord[2]/r;
	
	return pow( 1./8. * sqrt(2./M_PI) * r * exp( -r/2. ) * cos_th, 2 );
}

// Metropolis algorithm
void Metropolis( double (*distr)(vector<double>&), vector<double>& pos_in, int thr, int scar, string prob, double step, string file) {
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
	
	ofstream of;
	of.open(file);

	int dim = pos_in.size();
	vector<double> passo(dim), nuova_pos(dim);
	double alpha, r;
	int n = 0;
	
// Initial position (control)
	//for(int i=0; i<dim; i++)
	//	of << pos_in[i] << "    ";
	//of << endl;

	if(prob == "Unif")
		cout << "Metropolis algorithm with uniform probability transition" << endl;

	else if(prob == "Gauss")
		cout << "Metropolis algorithm with uniform probability transition" << endl;

	else {
		cerr << "Error: you must choose 'Unif' or 'Gauss'" << endl;
		return;
	}

// Make thr steps
	for(int k=0; k<thr; k++) {
	// Make a step in each direction (uniform or gaussian)
		for(int i=0; i<dim; i++) {
			if(prob == "Unif")
				passo[i] = rnd.Rannyu(-step, step);

			else
				passo[i] = rnd.Gauss(0., step);

			nuova_pos[i] = pos_in[i] + passo[i];
			//cout << nuova_pos[i] << endl;
		}
		
	// 'Transition probability'
		//cout << distr(nuova_pos) << endl;
		alpha = distr(nuova_pos)/distr(pos_in);

	// Make the step or reject
	// If make print on file
		if(alpha >= 1.) {
			for(int i=0; i<dim; i++) {
				pos_in[i] = nuova_pos[i];
				if( k>=scar )
					of << pos_in[i] << "    ";
			}
			if( k>=scar )
				of << endl;
			n++;
		}

		else {
			r = rnd.Rannyu();
			if(r <= alpha) {
				for(int i=0; i<dim; i++) {
					pos_in[i] = nuova_pos[i];
					if( k>=scar )
						of << pos_in[i] << "    ";
				}
				if( k>=scar )
					of << endl;
				n++;
			}

			else {
				for(int i=0; i<dim; i++)
					if( k>=scar )
						of << pos_in[i] << "    ";
				if( k>=scar )
					of << endl;
			}
		}
	}
	
	cout << "Sampling points accepted = " << (double)n/(double)thr * 100. << " %" << endl;
	cout << "The first " << scar << " have been rejected" << endl;

	of.close();
	rnd.SaveSeed();
}

// Take a file with coordinates and compute the distances from the origin
vector<double> Raggi(string infile) {
	ifstream inf;
	inf.open(infile);

	double x, y, z;
	vector<double> r;

	while(inf >> x >> y >> z)
		r.push_back( sqrt(x*x + y*y + z*z) );

	inf.close();
	
	return r;
}
