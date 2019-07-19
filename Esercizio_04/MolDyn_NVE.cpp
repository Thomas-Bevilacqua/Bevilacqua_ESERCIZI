/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include <vector>
#include "MolDyn_NVE.h"
#include "random.h"
#include "func.h"

using namespace std;

int main() {
	Input();             //Inizialization
	int nconf = 1;
	
	Eq_Temp();		//Temperature equilibration

	for(int istep=1; istep <= nstep; ++istep) {
		Move();				//Move particles with Verlet algorithm

		if(istep%iprint == 0)
			cout << "Number of time-steps: " << istep << endl;

		if(istep%10 == 0){
			Measure();    	 //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
			
			nconf += 1;
		}
	}
	
	//cout << utot.size() << endl;
	Measure_ave();	//block method

	ConfFinal();       //Write final configuration to restart
	
	return 0;
}


void Input(void){ //Prepare all stuff for the simulation

	ifstream ReadInput,ReadConf;

	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;

	//seed = 1;    //Set seed for random numbers
	//srand(seed); //Initialize random number generator

	ReadInput.open("input.dat"); //Read input

	ReadInput >> temp;

	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	ReadInput >> iprint;
	ReadInput >> mode;
	ReadInput >> Tst;

	if(mode != "actual" && mode != "old") {
		cerr << "Error: you should choose 'old' or 'actual' in the "
			"reading file mode..." << endl << endl;
		return;
	}

	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl;
	cout << "Desired temperature = " << Tst << endl;
	cout << "You choose the '" << mode << "' mode for reading the input" << endl << endl;
	
	ReadInput.close();

	//Prepare array for measurements
	iv = 0; //Potential energy
	ik = 1; //Kinetic energy
	ie = 2; //Total energy
	it = 3; //Temperature
	n_props = 4; //Number of observables

	//Read initial configuration
	cout << "Read initial configuration from file config.0 " << endl << endl;
	ReadConf.open("config.0");
	
	for (int i=0; i<npart; ++i) {
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}

	ReadConf.close();

//If you want the old configuration too, read it
	if(mode == "old") {
		cout << "Read old configuration from file old.0 " << endl << endl;
		ReadConf.open("old.0");
		
		for (int i=0; i<npart; ++i){
			ReadConf >> xold[i] >> yold[i] >> zold[i];
			xold[i] = xold[i] * box;
			yold[i] = yold[i] * box;
			zold[i] = zold[i] * box;
		}

		ReadConf.close();

	return;
	}

	else if(mode == "actual") {
//Random generator
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

//Prepare initial velocities
		cout << "Prepare random velocities with center of mass velocity equal to zero " << 
				endl << endl;
		double sumv[3] = {0.0, 0.0, 0.0};

		for (int i=0; i<npart; ++i) {
			vx[i] = rnd.Rannyu() - 0.5;
			vy[i] = rnd.Rannyu() - 0.5;
			vz[i] = rnd.Rannyu() - 0.5;

			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}

		for (int idim=0; idim<3; ++idim)
			sumv[idim] /= (double)npart;

		double sumv2 = 0.0, fs;

		for (int i=0; i<npart; ++i) {
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];	
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}

		sumv2 /= (double)npart;
		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor

		for (int i=0; i<npart; ++i) {
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = x[i] - vx[i] * delta;	//from actual to
			yold[i] = y[i] - vy[i] * delta;	//the previous
			zold[i] = z[i] - vz[i] * delta;
		}
		
		rnd.SaveSeed();
		return;
	}
}

//Find the equilibrium temperature
void Eq_Temp(void) {	//dT is the difference condition for stability
	cout <<"Reaching the desired temperature..." << endl << endl;
	
	double sumv[3];
	double dT = 0.0001;		//temperature spread
	int n = 0;
	ofstream of;

	of.open("Eq.temp");

	do {
		for(int i=0; i<3; i++)
			sumv[i] = 0.;
		
	//Make few step to evaluate the velocity
		for(int i=0; i<10; i++)
			Move();

	//Mean square velocity
		for(int i=0; i<npart; i++) { //Verlet integration scheme
			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}

		for (int i=0; i<3; i++)
			sumv[i] = sumv[i]/(double)npart;

		double sumv2 = 0.;

		for (int i=0; i<npart; i++) {
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}

		sumv2 =  sumv2/(double)npart;

	//Evaluate temperature
		temp = sumv2/3.;
		of << n << "    " << temp << endl;

	//no loop condition
		n++;
		if(n == 5000) {
			cout << "Warning! Can't reach the desired temperature" << endl;
			cout << "The system oscillates near T = " << temp << endl;
			cout << "You should modify the first temperature in input" << endl << endl;
			return;
		}
		
	} while(abs(Tst - temp) >= dT);

	cout << "Equilbrium temperature: " << temp << endl;
	cout << "This value will be used in the simulation" << endl;
	cout << "Velocities are scaled and positions saved" << endl << endl;

	of.close();

	return;

}

void Move(void){ //Move particles with Verlet algorithm
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i) { //Force acting on particle i
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i) { //Verlet integration scheme
		dxold[i] = xold[i];		
		dyold[i] = yold[i];		//t-dt
		dzold[i] = zold[i];

		xnew = Pbc( 2.0 * x[i] - dxold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - dyold[i] + fy[i] * pow(delta,2) );	//t+dt
		znew = Pbc( 2.0 * z[i] - dzold[i] + fz[i] * pow(delta,2) );

		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);	
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);		//v
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

		xold[i] = x[i];
		yold[i] = y[i];		//t
		zold[i] = z[i];

		x[i] = xnew;
		y[i] = ynew;		//t+dt
		z[i] = znew;
	}

	return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r) (sulla ip)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){			//per ogni particella
    if(i != ip){						//se la i-esima non è sè stessa
      dvec[0] = Pbc( x[ip] - x[i] );  	// distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){					//se supera il cut-off non interagisce
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement (istantanee)
  double v, t, vij, wij, w;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  Pres.open("output_pres.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  w = 0.;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
	   wij = 48. * (1.0/pow(dr,12) - 0.5/pow(dr,6));

//Potential energy
       v += vij;
	   w += wij;
     }
    }        
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy
	stima_pres = rho * stima_temp + w/(3. * vol);

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
	Pres << stima_pres << endl;

//save the current values in a vector
//for evaluating the means with blocks method
//in the following function
	utot.push_back(stima_pot);
	ktot.push_back(stima_kin);
	Ttot.push_back(stima_temp);
	etot.push_back(stima_etot);
	ptot.push_back(stima_pres);

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
	Pres.close();

    return;
}

void Measure_ave() {
	block_unc(utot, blocks, "epot.ave");
	block_unc(ktot, blocks, "ekin.ave");
	block_unc(Ttot, blocks, "temp.ave");
	block_unc(etot, blocks, "etot.ave");
	block_unc(ptot, blocks, "pres.ave");
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
