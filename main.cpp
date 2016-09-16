#include <iostream>
#include <stdio.h>

#include "cubature.h"
#include <gsl/gsl_math.h>


using namespace std;

///GLOBALS///

const unsigned g_dim = 3;


///PROTOTYPES///

int integrand_function(unsigned, const double*, void*, unsigned, double*);
int phasespace_integral(const double[], const double[], double * const, double * const);
int calculate_sigma_jet(double * const, double * const);



///FUNCTIONS///

int integrand_function(unsigned ndim,     //Variable dimension
                       const double *p_x, //Pointer to variable
                       void *p_fdata,     //Pointer to additional arguments
                       unsigned fdim,     //Output dimension
                       double *p_fval) {  //Pointer to output

    const double sigma = *((double *) p_fdata);
    double sum = 0;
    unsigned i;

    for (i = 0; i < ndim; ++i) sum += p_x[i] * p_x[i];
    p_fval[0] = exp(-sigma * sum);

    return 0; // success
}

int phasespace_integral(const double upper_limits [3],
                        const double lower_limits [3],
                        double * const p_value,
                        double * const p_error){

    const unsigned fdim=1;
    double fdata= 0.5;

    int not_success = pcubature(fdim, //Integrand dimension
              integrand_function, //Integrand function
              &fdata,              //Pointer to additional arguments
              g_dim,              //Variable dimension
              lower_limits,       //Variables minimum
              upper_limits,       //Variables maximum
              0,                  //Max n:o of function evaluations
              0,                  //Required absolute error
              1e-4,               //Required relative error
              ERROR_INDIVIDUAL,   //Enumerate of which norm is used on errors
              p_value,             //Pointer to output
              p_error);            //Pointer to error output

    if(not_success) {cout<<"Problem with integration"<<endl; return 1;}
    else {
            printf("Computed integral = %0.10g +/- %g\n", *p_value, *p_error);
            return 0;
    }
}

int calculate_sigma_jet(double * const p_value, double * const p_error){
    int not_success;
    const double upper_limits [3] = {2,2,2}, lower_limits [3] = {-2,-2,-2};

    not_success = phasespace_integral(upper_limits, lower_limits, p_value, p_error);
    return not_success;
}


///MAIN///

int main()
{
    int not_success;
    double sigma_jet, error;

    not_success = calculate_sigma_jet(&sigma_jet, &error);

    if(not_success) {cout<<"Bad stuff happened"<<endl; return 1;}
    else {
            printf("Sigma_jet = %0.10g +/- %g\n", sigma_jet, error);
            return 0;
    }
}
