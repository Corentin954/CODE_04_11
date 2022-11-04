// allocation dyn.
void freetab(void *ptr);
int **alloctab(int dim1, int dim2);
double **alloctabd(int dim1, int dim2);

// LOI D'ETAT
double fx(double V, int phase);
double THETA(double V, int phase);
double fEs(double V, int phase);
double fPs(double V, int phase);
double u(double V, double E, int phase);
double fS(double V, double E, int phase);
double fP(double V, double E, int phase);
double fT(double V, double E, int phase);
double fG(double V, double E, int phase);
double Es_prime(double V, int phase);
double Ps_prime(double V, int phase);
double dPdE(int phase);
double dPdV(double V, int phase);
double dTdE(int phase);
double dTdV(double V, int phase);
double d2PdV2(double V, int phase);
double d2TdV2(double V, int phase);

void coeff(int phase, double* K0, double* N0, double* gamma0, double* Crv, double* theta0, double* T0, double* P0, double* rho0, double* v0, double* E0, double* Sr, double* ur);

double fE_VP(double V, double P, int phase);

double Ps_prime2(double V, int phase);
double Es_prime2(double V, int phase);
double THETA_prime(double V, int phase);
double THETA_prime2(double V, int phase);