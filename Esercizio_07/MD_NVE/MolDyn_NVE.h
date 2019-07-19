/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000;
int n_props, nbins;
int iv,ik,it,ie, igofr;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres, stima_gofr;
double glob_av[m_props], glob_av2[m_props];
double vtail,ptail,bin_size;
double walker[m_props];

// averages
double acc,att;
int blocks = 100;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double dxold[m_part], dyold[m_part], dzold[m_part];
std :: vector<double> etot, utot, ktot, Ttot, ptot;

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;
std :: string mode;
double Tst;
double teq = 100.;

//functions
void Input(void);
void Eq_Temp(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void Measure_ave(void);
void Measure_gofr(void);
double pressure(void);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
