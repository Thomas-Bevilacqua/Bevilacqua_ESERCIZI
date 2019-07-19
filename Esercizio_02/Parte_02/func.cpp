#include "func.h"

// Statistics
double Mean(vector<double> data) {
	double sum = 0.;
	int size = data.size();

	for(int i=0; i<size; i++)
		sum += data[i];

	return sum/(double)size;
}

double Sigma(vector<double> data) {
	double sum = 0.;
	double mean = Mean(data);
	int size = data.size();

	for(int i=0; i<size; i++)
		sum += pow(data[i] - mean, 2);

	return sqrt(sum/(double)size);
}
