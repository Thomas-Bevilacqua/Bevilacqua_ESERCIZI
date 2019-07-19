#include "genetic.h"

// Print
void Print(vector<double> vett) {
	int n = vett.size();

	for(int i=0; i<n; i++)
		cout << vett[i] << "    ";
	cout << endl;
}

// Print coordinates on a file
void Print(vector<double> x, vector<double> y, string file) {
	if(x.size() != y.size() ) {
		cerr << "Error: you should put two vector with the same dimension" << endl;
		cout << x.size() << " =/= " << y.size() << endl;
		return;
	}

	int n = x.size();
	ofstream of;

	of.open(file);

	for(int i=0; i<n; i++)
		of << x[i] << "    " << y[i] << endl;
	of << x[0] << "    " << y[0] << endl;
}

// Print the ordered cities
void Print_best_path(vector< vector<double> > xpop, vector< vector<double> > ypop, string file) {
	if(xpop.size() != ypop.size()) {
		cerr << "Error: the vectors have different size" << endl;
		return;
	}
	
	int pop = xpop.size();
	int b = 0;
	double min = 10000.;
	
	for(int i=0; i<pop; i++)
		if( L2(xpop[i], ypop[i]) < min ) {
			min = L2(xpop[i], ypop[i]);
			b = i;
		}
	
	Print(xpop[b], ypop[b], file);
}

// Print average path length (on the best half of the population)
vector<double> Stat(vector< vector<double> > xpop, vector< vector<double> > ypop) {
	if(xpop.size() != ypop.size()) {
		cerr << "Error: the vectors have different size" << endl;
		return xpop[0];
	}
	
	int pop = xpop.size();
	double sum = 0., sum2 = 0.;
	vector<double> ris(2);
	vector<double> d(pop);
	
	for(int i=0; i<pop; i++)
		d[i] = L2(xpop[i], ypop[i]);
	
	sort( d.begin(), d.end() );
	
	pop = pop/2;
	for(int i=0; i<pop; i++) {
		sum += d[i];
		sum2 += pow(d[i], 2);
	}
	
	ris[0] = sum/(double)pop;
	ris[1] = sqrt( sum2/(double)pop - ris[0]*ris[0] );

	return ris;
}


// Check
void Check(vector<double> vett) {
	int n = vett.size();

	for(int i=0; i<n; i++)
		for(int j=i+1; j<n; j++)
			if( vett[i] == vett[j] ) {
				cerr << "Error: same cities in positions " << i << " , " << j << endl;
				Print(vett);
				return;
			}
	
	cout << "This vector has all different elements" << endl;
}

void Check(vector<double> x, vector<double> y) {
	if(x.size() != y.size() ) {
		cerr << "Error: you should put two vector with the same dimension" << endl;
		cout << x.size() << " =/= " << y.size() << endl;
		return;
	}

	int n = x.size();
	
	for(int i=0; i<n; i++)
		for(int j=i+1; j<n; j++)
			if( x[i] == x[j] && y[i] == y[j]) {
				cerr << "Error: same cities in positions " << i << " , " << j << endl;
				Print(x);
				Print(y);
				return;
			}
	
	cout << "All points are different" << endl;
}	


// Given x construct the y on a circle with radius = r
void Circle(vector<double>& x, vector<double>& y, double r, Random rnd) {
	if(x.size() != y.size() ) {
		cerr << "Error: you should put two vector with the same dimension" << endl;
		cout << x.size() << " =/= " << y.size() << endl;
		return;
	}

	int n = x.size();

	for(int i=0; i<n; i++) {
		double a = rnd.Rannyu(0., 1.);
		if(a < 0.5)
			y[i] = sqrt(pow(r, 2) - pow(x[i], 2) );

		else
			y[i] = -sqrt(pow(r, 2) - pow(x[i], 2) );
	}
}


// MUTATIONS: don't change the first city
// Pair permutation
void Switch(vector<double>& x, int in, int fin) {
	if( fin >= (int)x.size() ) {
		cerr << "Error: the city " << fin << " doesn't exist" << endl;
		return;
	}

	swap( x[in], x[fin] );
}


// Shift of every city from a selected position
void Shift(vector<double>& x, int n, int in) {
	int dim = x.size() - in;
	vector<double> dep(dim);

	for(int i=0; i<dim; i++)
		dep[i] = x[i + in];

	for(int i=0; i<dim; i++)
		x[i + in] = dep[ BC(i-n, dim) ];
}


// Permutation of m cities, [in, in+m] and [fin, fin+m]
void Permutation(vector<double>& x, int m, int in, int fin) {
	int dim = x.size();

	if(m > dim/2) {
		cerr << "Error: There are not enough cities to permute" << endl;
		return;
	}

	if(fin < in + m || fin > dim - m)
		{
		cerr << "Error: There are not enough cities to permute" << endl;
		return;
	}

	vector<double> dep(m);

	for(int i=0; i<m; i++) {
		dep[i] = x[i + in];
		x[i + in] = x[i + fin];
		x[i + fin] = dep[i];
	}
}


// Inversion of m cities
void Inversion(vector<double>& x, int in, int fin) {
	int dim = x.size();

	if(fin > dim) {
		cerr << "Error: Can't reverse the order, city " << fin << " doesn't exist" << endl;
		return;
	}

	reverse(x.begin() + in, x.begin() + fin);
}


// Make a randomly chosen mutation: prob in [0:1]
void Make_Mutation(vector<double>& x, vector<double>& y, Random rnd, double prob) {
	if(x.size() != y.size() ) {
		cerr << "Error: the vectors have different size" << endl;
		return;
	}
	
	double dim = x.size();
	if(prob < 0.1) {
		
		if( prob <= 0.025 ) {
			int in = rnd.Rannyu(1., dim);
			int fin = rnd.Rannyu(1., dim);

			// Don't switch the same city
			while(in == fin)
				fin = rnd.Rannyu(1., dim);

			Switch(x, in, fin);
			Switch(y, in, fin);
		}

		else if(prob > 0.025 && prob <= 0.050) {
			int in = rnd.Rannyu(1., dim);
			int sh = rnd.Rannyu( -(double)(dim-in) - 1., (double)(dim-in) + 1. );

			while(sh == 0)
				sh = rnd.Rannyu( -(double)(dim-in) - 1., (double)(dim-in) + 1. );
			
			Shift(x, sh, in);
			Shift(y, sh, in);
		}

		else if(prob > 0.050 && prob <= 0.075) {
			int N = 0.5 * (double)dim;

			int m = rnd.Rannyu(2., (double)N);
			int in = rnd.Rannyu(1., (double)(N - m));
			int fin = rnd.Rannyu(N, (double)(dim-N));

			Permutation(x, m, in, fin);
			Permutation(y, m, in, fin);
		}
		
		else if(prob > 0.075) {
			int in = rnd.Rannyu(1., (double)dim);
			int fin = rnd.Rannyu((double)in, (double)(dim+1));

			Inversion(x, in, fin);
			Inversion(y, in, fin);
		}
	}
	
	else{ }

}


// Path length
double L2(vector<double> x, vector<double> y) {
	if(x.size() != y.size()) {
		cerr << "Error: the vectors have different size" << endl;
		return -1.;
	}
	
	int dim = x.size();
	double sum = 0.;

	for(int i=0; i<dim; i++)
		sum += pow( x[ BC(i+1, dim) ]-x[i], 2) + pow( y[ BC(i+1, dim) ]-y[i], 2);

	return sum;
}

// Find the shortest path in a population
double Shortest(vector< vector<double> > xpop, vector< vector<double> > ypop) {
	if(xpop.size() != ypop.size()) {
		cerr << "Error: the vectors have different size" << endl;
		return -1;
	}
	
	int pop = xpop.size();
	double min = 10000.;
	
	for(int i=0; i<pop; i++)
		if( L2(xpop[i], ypop[i]) < min )
			min = L2(xpop[i], ypop[i]);


	return min;
}
	


// Crossover: writes sons over dad and mum (cuts from n+1)
// give a couple (x,y) to crossover to (xm,yd)
void Crossover(vector<double>& x, vector<double>& y, vector<double>& xx, vector<double>& yy, double prob, int pos) {
	if(x.size() != y.size()) {
		cerr << "Error: the vectors have different size" << endl;
		return;
	}

	if(prob <= 0.5) {
		vector<double> depx = x, depy = y, dx = xx, dy = yy;
		
		Write(x, dx, pos);
		Write(xx, depx, pos);
		Write(y, dy, pos);
		Write(yy, depy, pos);
	}
}


// New generation: start from (x,y) population and then overwrite
void Sons(vector< vector<double> >& xpop, vector< vector<double> >& ypop, Random rnd, int nmut) {
	if(xpop.size() != ypop.size()) {
		cerr << "Error: the vectors have different size" << endl;
		return;
	}
	
	int pop = xpop.size();
	int ncity = (xpop[0]).size();
	vector< vector<double> > depx, depy;
	vector<double> d(pop);
	//cout << ncity << endl;
	
	// Compute the path length and sort the population
	map<int, double> p;
	for(int i=0; i<pop; i++) {
		d[i] = L2(xpop[i], ypop[i]);
		p[i] = d[i];
	}
	
	sort( d.begin(), d.end() );

	//cout << d[0] << endl;
	//for(int i=0; i<pop; i++)
	//	cout << d[i] << "	" << endl;
	//cout << endl;
	
	// Make an optimized population
	double acc = (double)pop/2.;
	int c = 0;
	while(c < pop) {
		int k = acc * pow(rnd.Rannyu(), 2);
		for(int i=0; i<pop; i++)
			if(d[k] == p[i]) {
				depx.push_back(xpop[i]);
				depy.push_back(ypop[i]);
				c++;
				break;
			}
	}
	
	//cout << depx.size() << endl;

	for(int k=0; k<pop-1; k++) {
		// Extract two parents
		int im = rnd.Rannyu(0., (double)pop);
		int id = rnd.Rannyu(0., (double)pop);
		
		while(im == id)
			id = rnd.Rannyu(0., (double)pop);
		//cout << im << "		" << id << endl;
		
		// Crossover and mutate
		double prob = rnd.Rannyu();
		int pos = rnd.Rannyu(1., (double)(ncity-1));
		Crossover(depx[im], depy[im], depx[id], depy[id], prob, pos);

		for(int i=0; i<nmut; i++) {
			double prob1 = rnd.Rannyu();
			double prob2 = rnd.Rannyu();
			Make_Mutation(depx[im], depy[im], rnd, prob1);
			Make_Mutation(depx[id], depy[id], rnd, prob2);
		}

		// Overwrite the sons in the old population
			xpop[k] = depx[im];
			ypop[k] = depy[im];
			xpop[k+1] = depx[id];
			ypop[k+1] = depy[id];
	}
}


// Simulated annealing
void SA(vector<double>& x, vector<double>& y, vector<double> depx, vector<double> depy, double prob, double beta, int& thr) {
	if(x.size() != y.size() || depx.size() != depy.size() ) {
		cerr << "Error: the vectors have different size" << endl;
		return;
	}
	
	double ld, l, p;
	ld = L2(depx, depy);
	l = L2(x, y);
	
// If better substitute (Metropolis)
	if(ld - l < 0.) {	// Boltzmann p > 1
		x = depx;
		y = depy;
		thr++;
	}

	else {
		p = exp(-beta * (ld - l));		// Boltzmann
		if(prob <= p) {
			x = depx;
			y = depy;
			thr++;
		}
	}
}


// Boundary conditions
int BC(int i, int dim) {
	if(i >= dim)
		i = i - dim;

	else if(i < 0)
		i = i + dim;

	return i;
}


// Overwrite (for crossover) from nth position included
void Write(vector<double>& x, vector<double>& y, int n) {
	int dim = x.size();
	vector<double> dep = x;
	int c = n;
	
	for(int i=0; i<dim; i++)
			for(int j=n; j<dim; j++)
				if( dep[j] == y[i] ) {
					x[c] = y[i];
					c++;
				}
}

