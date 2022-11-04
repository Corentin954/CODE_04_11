// allocation dyn.
void freetab(void *ptr);
int **alloctab(int dim1, int dim2);
double **alloctabd(int dim1, int dim2);

int print_simulation(int tst, int sch, int so, double CFL, int nx, int aff, int Nmax, int kin_fix, int dpi, 
                      int q, double Cq, double Cl, int cin_phase, int Ntau, double dt_first, double Ttau, 
                      int sch_cin_phase, int tst_multi, int EOS, double epsilon);

// LOI D'ETAT
double epsilonEOS(int tst, double tau, double p);
void EGS(int tst, double tau, double epsilon, double* Pres, double* Cres);
void G_EGS(int tst, double tau, double epsilon, double* pG);
void EGS_cin_phase(int tst, double tau, double epsilon, double* pP, double* pC, int CIN_PHASE, double dx, double dt, double Ttau, double* Un, double Pn, double Vn, double Cn );
double G_BIZ_cin_phase(double tau, double epsilon, double* pG, int CIN_PHASE, double dx, double dt, double Ttau, double* Un, double Pn, double Vn, double Cn);

// Montée en odre spatiale
double phi(int so, double* Z, double** Coef);
double phidiv(int so, double* RZ, double* RHO, double** Coef);
double fxc(int so,double* X, double** rk);
double delta(int so,int ind, double* Z, double** dk);
double YdZ(int so, int cas, double* P, double* Q, double* u, double** dk, double** Ckbar);
double phiQK(int so, double* DK, double** Q);
void initCoef(double** Ck, double** Ckbar, double** dk, double** Qbar, double** rk);

// intialisation 
int u0_multi(int tst, int ind, int tst_air, int tst_etain, double a, double b,double x, double* W);
int u0_bar_multi(int tst, int ind, int ind_air, int ind_etain, int tst_air, double a, double b, double xg, double xd, double* W);

// Table RUNGE - KUTTA
int initButcher(int sch, double** A, double* THETA, double* ALPHA);
int condlim(int ind_cond, int jdeb, int jfin, int nbghosts, double* Z, int iu);
int condlim_Gtot(int ind_cond, int jdeb,int jfin, int nbghosts, double* Z, int ind);
// The Minmod limiter
double MinMod(double r);
// The Superbee limiter
double Superbee(double r);
// The Christensen limiter
double Christensen(double r, double arg);
double vanLeer(double r_p, double r_m);
// pseusdo-viscosity
double qvis(int q, double Cq, double Cl, double* U, double tau, double c, double* DM);
// print
//int print_sol(int ind_u, int ideb, int ifin, double* X, double* Xc, double* TAU, double* U, double* E, double* P, double* EPS);
int print_sol(int ind_u, int ideb, int ifin, double* X, double* Xc, double* RHO, double* TAU, double* U, double* E, double* EPS, double* P, double* T, double* C, double* G, double* S, double** matFM);
int print_sol_RKint(int ind_u, int ideb, int ifin, double* X, double* Xc, double* RHO, double* TAU, double* U, double* E,  double* EPS, double* P, double* T, double* C, double* G, double* S, double* RTAU, double* RU, double* REPS);
int print_nrj(int n, double* Etot, double* IMPUL, double* Mtot, double* VOL);
double pastemps(double T, double t, double dt, double dt_old);
//schemes
int funcGtot_multi(int sch, int tst, double T, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int EOS,
	               int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);

int funcRKint_multi(int sch, int tst, double T, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int q, double Cq, double Cl, int z, int dpi, int so, int EOS,
	                int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                    double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                    double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);

int funcBBC_JCP2009_multi(int tst, double T, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int q, double Cq, double Cl, int EOS,
	        	          int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
             	          double* tabVx, double* tabEx, double* tabVy, double* tabEy,
        		          double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);

int funcBBC_PRED_CORR_multi(int tst, double T, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int q , double Cq, double Cl, int EOS,
	               			int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
      		                double* tabVx, double* tabEx, double* tabVy, double* tabEy,
              			    double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);

int funcBBC_RK2av_multi(int tst, double T, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int q, double Cq, double Cl, int EOS,
	            	    int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   		double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   		double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);

int funcvNR_multi(int tst, double T, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int q, double Cq, double Cl, int EOS,
	              int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                  double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                  double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);

int funcKOVA_multi(int tst, double T, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int z,
	               int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);

int funcKOVAo2_multi(int sch, int tst, double Tf, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int EOS, double zeta,
                          int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                          double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                          double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);
 
int funcGodGoHy_multi(int sch, int tst, double T, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int so, int EOS,
                   int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);

int funcGodGoHy3_multi(int sch, int tst, double Tf, double a, double b, int nx, double epsilon, int Nmax, double CFL, int aff, int ind_cond, int iu, int Ntau, double dt_first, double Ttau, int cin_phase, int sch_cin_phase, int tst_multi, int tst_air, int EOS,
                     int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, int** tabZONE, int* tab_Nb_ZONE, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);
    
// Equation d'etat
int lec_tabulation(int nb_points, int Nv_dis, int Ne_dis, int** tabZONE, int* tab_Nb_ZONE, double* tabVx, double* tabEx, double* tabVy, double* tabEy,
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC);

int fPTC(double V, double E, double* pP, double *pc, double* pT, double* px,
         int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
         double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
         double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE );

int epsilon_etain(double V, double P, double T, double* pE, double* lambda, int nb_points,
                  double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC );
  

int fVE(int phase, double P, double T, double* pV, double* pE);


void fnb_pts(int nb_points, int* pNTG1, int* pNPB1, int* pNPB2, int* pNTD2 , int* pNTD3, int* pNPH3, int* pNTG3, 
            int* pN12 ,int* pN13, int* pN23, int* pNA, int* pNB, int* pNC, int* pNAB, int* pNAC, int* pNBC );

int fPn1_etain_cin_phase(int tst, int EOS, int cin_phase, int sch_cin_phase, double epsilon, double TAUn, double TAUn1, double Pn, double Qn, double EPSn, double* pEPSn1, double* pPn1, double* lambda, int Ntau, double dX, double dt, double Ttau,
               int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
               double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
               double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE );

int fPTC_cin_phase(int affichage, int EOS, int CIN_PHASE, int sch_cin_phase, double epsilon, double V, double E, double* pP, double *pc, double* pT, double* pVAR, double* lambda_in, double* lambda_out, 
                   int Ntau, double dx, double dt, double Ttau, int* pzone, int* n_iter, double* pcritere, double* p_test_consis_thermo,
                   double dUdx, double Pn, double Vn, double Cn,// cinétique
                   int nb_points, int Nv_dis, int Ne_dis, double Vmin_dis, double Vmax_dis, double Emin_dis, double Emax_dis, 
                   double** SA, double** SB, double** SC, double** SAB, double** SAC, double** SBC, 
                   double* tabVx, double* tabEx, double* tabVy, double* tabEy, int** tabZONE, int* tab_Nb_ZONE );

int init_fluides(int tst_multi, int* indFLUIDE, int N, int ind_air, int ind_etain);

int c_celerite_son_init(double P, double V, double E, double* lambda, double *pc);

int G_S_cin_3_phase(double epsilon, double P, double V, double E, double* lambda, double *pc, double* pG, double* pS);

int c_celerite_son_cin_3_phase_all_dom(int Nmax, double epsilon, double P, double V, double E, double dV, double* lambda, double *pc);




// Schémas GoHy
double fpw_to_av(double* Ypw);
double fav_to_pw(double* Yav);
double f_grid(double* Xd);
double fdtu(double* RHO, double* P, double dx);
double fdtp(double* RHO, double* RHOC2, double* U, double dx);
double fdttu(double* RHO, double* RHOC2, double* U, double dx);
double fdttp(double* RHO, double* RHOC2, double* G, double* P, double* U, double dx);
double f_inter(int so, double* Yc, double** rk);

// Prise d'erreur 
int sol_err_exacte(int cas_test, int nx, int ind_u);
