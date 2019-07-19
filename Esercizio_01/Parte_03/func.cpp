#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

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

//funzione seno
double Sin(double x, double y) {
	return y/sqrt(x*x + y*y);
}
