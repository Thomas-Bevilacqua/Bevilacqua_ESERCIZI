#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

//Funzione per il calcolo dell'incertezza statistica col
//blocking method. Lo metto qui per comodità
void block_unc(vector<double>& data, int blocks, string file) {
	int throws = data.size();
	double sum, sum2;
	int ratio = throws/blocks;
	vector<double> mean, mean2;

	ofstream outfile;
	outfile.open(file);

//calcolo la media per ogni block
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

//calcolo la media e il relativo errore al crescere
//del numero di block compresi poi stampo su file i risultati
	for(int i=0; i<blocks; i++) {
		sum = 0.;
		sum2 = 0.;
		for(int j=0; j<i+1; j++) {
			sum += mean[j];
			sum2 += mean2[j];
		}

		if( outfile.is_open() )
			outfile << (i+1)*ratio << "    " << sum/(double)(i+1) << "    "
			<< sqrt(sum2/(double)(i+1) - pow(sum/(double)(i+1), 2) )/sqrt( (double)(i+1) ) 
			<< endl;

		else
			cerr << "Error: can't open " << file << "..." << endl;
	}

	outfile.close();
}
