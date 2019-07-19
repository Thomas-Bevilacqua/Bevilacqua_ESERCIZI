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

// Potential
double V(double x) {
	return pow(x, 4) - 2.5 * pow(x,2);
}


// Trial wave function
// Usage: x, mu, sigma
double psi_T(double x, vector<double>& param) {

	double mu = param[0];
	double sigma = param[1] * param[1];
	
	double n1 = exp( -pow(x-mu, 2)/(2. * sigma) );
	double n2 = exp( -pow(x+mu, 2)/(2. * sigma) );

	return n1 + n2;
}

double psi2_T(double x, vector<double>& param) {
	return pow(psi_T(x, param), 2);
}


// Energy
double En(double x, vector<double>& param) {
	double mu = param[0];
	double sigma = param[1] * param[1];

	double e1 = pow(x-mu, 2)/sigma;
	double e2 = pow(x+mu, 2)/sigma;
	
	double psi_sec = ( -psi_T(x, param) + e1 * exp(-e1/2.) + e2 * exp(-e2/2.) )/sigma;
	
	return -0.5 * psi_sec/psi_T(x, param) + V(x);
}


double mean(vector<double>& data) {
	double sum = 0.;
	for(double elem : data)
		sum += elem;

	return sum/(double)data.size();
}
	

// Metropolis algorithm
vector<double> Metropolis( double (*distr)(double, vector<double>&), double pos_in, vector<double>& param, int thr, int scar, double step, string file ) {
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

	vector<double> ris;
	double passo, nuova_pos;
	double alpha, r;
	int n = 0;

// Make thr steps
	for(int k=0; k<thr; k++) {
	// Make a step in each direction
		passo = rnd.Rannyu(-step, step);
		nuova_pos = pos_in + passo;
		
	// 'Transition probability'
		//cout << distr(nuova_pos) << endl;
		alpha = distr(nuova_pos, param)/distr(pos_in, param);

	// Make the step or reject
	// If make print on file and save
		if(alpha >= 1.) {
			pos_in = nuova_pos;
			if( k>=scar ) {
				of << pos_in << endl;
				ris.push_back(pos_in);
			}
			n++;
		}

		else {
			r = rnd.Rannyu();
			if(r <= alpha) {
				pos_in = nuova_pos;
				if( k>=scar ) {
					of << pos_in << endl;
					ris.push_back(pos_in);
				}
				n++;
			}

			else
				if( k>=scar ) {
					of << pos_in << endl;
					ris.push_back(pos_in);
				}
		}
	}
	
	cout << "Sampling points accepted = " << (double)n/(double)thr * 100. << " %" << endl;
	cout << "The first " << scar << " have been rejected" << endl << endl;

	of.close();
	rnd.SaveSeed();

	return ris;
}
