#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include "random.h"

using namespace std;

void Check(vector<double> vett);
void Check(vector<double> x, vector<double> y);

void Print(vector<double> vett);
void Print(vector <double> x, vector <double> y, string file);
void Print_best_path(vector< vector<double> > xpop, vector< vector<double> > ypop, string file);
vector<double> Stat(vector< vector<double> > xpop, vector< vector<double> > ypop);

void Circle(vector<double>& x, vector<double>& y, double r, Random rnd);

double L2(vector<double> x, vector<double> y);
double Shortest(vector< vector<double> > xpop, vector< vector<double> > ypop);


// Mutations
void Switch(vector<double>& x, int in, int fin);
void Shift(vector<double>& x, int n, int in);
void Permutation(vector<double>& x, int m, int in, int fin);
void Inversion(vector<double>& x, int in, int fin);

void Make_Mutation(vector<double>& x, vector<double>& y, Random rnd, double prob);

// Crossover
void Crossover(vector<double>& x, vector<double>& y, vector<double>& xm, vector<double>& yd, double prob, int pos);

// Make a new generation
void Sons(vector< vector<double> >& xpop, vector< vector<double> >& ypop, Random rnd, int nmut);


int BC(int i, int dim);
void Write(vector<double>& x, vector<double>& y, int n);
