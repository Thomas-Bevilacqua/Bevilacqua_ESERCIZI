#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "random.h"

using namespace std;

void block_unc(vector<double>&, int, string);

double psi_100(vector<double>&);
double psi_210(vector<double>&);

void Metropolis(double (*distr)(vector<double>&), vector<double>&, int, int, string, double, string);

vector<double> Raggi(string);
