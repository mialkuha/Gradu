#include <iostream>
#include <stdio.h>

#include <complex>

#include "cubature.h"
#include <gsl/gsl_math.h>

#include "LHAPDF/LHAPDF.h"


using namespace std;
using namespace LHAPDF;

///GLOBALS///

const unsigned g_dim = 3;
const string g_pdfsetname = "HERAPDF20_LO_EIG";
const int g_pdfsetmember = 1;
const double g_s = 10;


///PROTOTYPES///

int integrand_function(unsigned, const double*, void*, unsigned, double*);
int phasespace_integral(const double[], const double[], double * const, double * const);
int calculate_sigma_jet(double * const, double * const);
double f_ses(const double * const, const double * const, const PDF*);
double s_hat(const double * const, const double * const);
double t_hat(const double * const, const double * const, const double * const);
double u_hat(const double * const, const double * const, const double * const);
double sigma_gg_gg(const double * const, const double * const, const double * const);


///FUNCTIONS///
double sigma_gg_gg(const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat){
    const double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return 4.5*(3-(u*t)/(s*s)-(u*s)/(t*t)-(s*t)/(u*u));
}

double u_hat(const double * const p_x1, const double * const p_x2, const double * const p_kt2){
    const double x1 = *p_x1, x2 = *p_x2, kt2 = *p_kt2;
    const double dummy = x1*x2*(x1*x2 -4*(kt2/g_s));
    return -0.5*g_s*(sqrt(dummy) +x1*x2);
}

double t_hat(const double * const p_x1, const double * const p_x2, const double * const p_kt2){
    const double x1 = *p_x1, x2 = *p_x2, kt2 = *p_kt2;
    const double dummy = x1*x2*(x1*x2 -4*(kt2/g_s));
    return 0.5*g_s*(sqrt(dummy) -x1*x2);
}

double s_hat(const double * const p_x1, const double * const p_x2){ return *p_x1 * *p_x2 * g_s; }

double f_ses(const double * const p_x, const double * const p_q2, const PDF* p_pdf){
    double sum = 0;
    for (int flavour=-6; flavour<=6; ++flavour){
        sum += (flavour == 0) ? p_pdf->xfxQ2(flavour, *p_x, *p_q2) : 4 * p_pdf->xfxQ2(flavour, *p_x, *p_q2) / 9;
    }
    //cout<<"x="<<*p_x<<", q2="<<*p_q2<<" SES="<<sum<<endl;
    return sum;
}


int integrand_function(unsigned ndim,     //Variable dimension
                       const double *p_x, //Pointer to variable
                       void *p_fdata,     //Pointer to additional arguments
                       unsigned fdim,     //Output dimension
                       double *p_fval) {  //Pointer to output

    const PDF* p_pdf = (PDF *) p_fdata;
    const auto kt2 = p_x[0], x1 = p_x[1], x2 = p_x[2];
    const auto s = s_hat(&x1, &x2);
    const auto t = t_hat(&x1, &x2, &kt2);
    const auto u = u_hat(&x1, &x2, &kt2);
    const auto subprocess_cs = sigma_gg_gg(&s, &t, &u);
    const auto jacobian = - 1/sqrt(x1*x2*(x1*x2 -4*(kt2/g_s)));

    p_fval[0] = 0.5 * x1 * f_ses(&x1, &kt2, p_pdf)
                    * x2 * f_ses(&x2, &kt2, p_pdf)
                    * subprocess_cs * jacobian;

    cout<<"kt2="<<kt2<<", x1="<<x1<<", x2="<<x2<<", s="<<s<<", t="<<t<<", u="<<u<<", cs="<<subprocess_cs<<", J="<<jacobian<<", I="<<p_fval[0]<<endl;
    return 0; // success
}

int phasespace_integral(const double upper_limits [3],
                        const double lower_limits [3],
                        double * const p_value,
                        double * const p_error){

    const unsigned fdim=1;
    PDF* p_pdf = mkPDF(g_pdfsetname, g_pdfsetmember);

    int not_success = pcubature(fdim, //Integrand dimension
              integrand_function, //Integrand function
              p_pdf,              //Pointer to additional arguments
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
    const double upper_limits [3] = {10,1,1}, lower_limits [3] = {1,1E-5,1E-5};

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
