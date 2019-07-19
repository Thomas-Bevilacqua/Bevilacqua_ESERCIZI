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
	
// Initial position and mu, sigma
	vector<double> param = {1., 0.5};

	int thr = 11000;
	int scar = 1000;

// Sampling test
	cout << "Sampling test of |psi|^2 and <H>" << endl;
	cout << "using mu = " << param[0] << " and sigma = " << param[1] << endl << endl;

	cout << "|psi|^2:" << endl;
	vector<double> psi = Metropolis(psi2_T, 0., param, thr, scar, 2.1, "Risultati/psi2_cam.prova");

	vector<double> ene;

	for(double elem : psi)
		ene.push_back( En(elem, param) );

	block_unc(ene, 100, "En.prova");

// Minimization of the parameters
// Fix the initial parameters and then evaluate
// the function as a grid
	double sigma, mu;
	double mu_m = 0.7, mu_M = 1., mu_s = 0.01;
	double sigma_m = 0.4, sigma_M = 0.8, sigma_s = 0.01;

	double E_min = 10000., sigma_min = sigma_m, mu_min = mu_m;
	double sum;
	
	ofstream of;
	of.open("Risultati/matrix.dat");

	for(mu=mu_m; mu<mu_M; mu=mu+mu_s) {
		param[0] = mu;
		for(sigma=sigma_m; sigma<sigma_M; sigma=sigma+sigma_s) {
			sum = 0.;
			param[1] = sigma;
			
			//cout << param[0] << "	" << param[1] << endl;
			psi = Metropolis(psi2_T, 0., param, thr, scar, 2.3, "Risultati/Psi.min");

			for(double elem : psi)
				sum += En(elem, param);
			
			sum = sum/(double)psi.size();

			if(sum < E_min) {
				E_min = sum;
				sigma_min = sigma;
				mu_min = mu;
			}

			of << mu << "    " << sigma << "    " << sum << endl;
		}
	}

	of.close();
			
	cout << "Minimized parameters:" << endl;
	cout << "E = " << E_min << " mu = " << mu_min << " sigma = " << sigma_min << endl << endl;

// Assign the parameters
	param[0] = mu_min;
	param[1] = sigma_min;
	
	psi = Metropolis(psi2_T, 0., param, thr, scar, 2.3, "Matrix/Psi.min");

	for(int i=0; i<thr-scar; i++)
		ene[i] = En(psi[i], param);

	block_unc(ene, 100, "Matrix/En.min");
			
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
