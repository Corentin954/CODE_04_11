#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

// sh : compil_test_integration.sh

double f(double x){
    double res;
	if(x<=3./4){  res=-1.; } 
	else {  res=1.; } 
	return res;
}

int main(){
  double a=-1, b=1; 
  int res; 

  size_t n=10; 
  gsl_integration_glfixed_table * TTT = gsl_integration_glfixed_table_alloc(n);

  int gsl_integration_glfixed_point(double a, double b, size_t i, double *xi, double *wi, const gsl_integration_glfixed_table *t); 

  double xi,wi;   
  double *tabx=malloc(n*sizeof(double)); 
  double *tabw=malloc(n*sizeof(double));
  double *tabf=malloc(n*sizeof(double));

  double I=0;
  for(int i=0; i<n; i++){
    res= gsl_integration_glfixed_point(a, b, i, &xi, &wi, TTT); 
    if(res){return 1;}
    tabx[i]=xi; 
    tabw[i]=wi; 
    tabf[i]=f(xi);
    I+=wi*f(xi); 

    printf("i=%d | x=%e, w=%e, f=%e, I=%e\n",i, tabx[i], tabw[i],tabf[i],I );
  }
  
  printf("  I = %.10e\n",I );



}