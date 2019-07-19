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
	
// Starting parameters
	vector<double> pos_in = {0., 0., 0.};
	int thr = 1E6;
	ifstream p100, p210;
	int bl = 100;
	int scar = 1000;

// Sampling and print
	cout << "psi_100" << endl;
	Metropolis(*psi_100, pos_in, thr, scar, "Unif", 1.2, "Risultati/psi_100.uniform");
	cout << endl;
	cout << "psi_210" << endl;
	Metropolis(*psi_210, pos_in, thr, scar, "Unif", 3., "Risultati/psi_210.uniform");
	
	cout << endl;

	cout << "psi_100" << endl;
	Metropolis(*psi_100, pos_in, thr, scar, "Gauss", 0.8, "Risultati/psi_100.norm");
	cout << endl;
	cout << "psi_210" << endl;
	Metropolis(*psi_210, pos_in, thr, scar, "Gauss", 1.8, "Risultati/psi_210.norm");

// Compute the average distances and their errors
	vector<double> r1 = Raggi("Risultati/psi_100.uniform");
	vector<double> r2 = Raggi("Risultati/psi_210.uniform");
	vector<double> r1G = Raggi("Risultati/psi_100.norm");
	vector<double> r2G = Raggi("Risultati/psi_210.norm");

	block_unc(r1, bl, "Risultati/r100.uniform");
	block_unc(r2, bl, "Risultati/r210.uniform");
	block_unc(r1G, bl, "Risultati/r100.norm");
	block_unc(r2G, bl, "Risultati/r210.norm");
	
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
