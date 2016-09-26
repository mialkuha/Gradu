#include <iostream>
#include <stdio.h>

#include <complex>
#include <gsl/gsl_math.h>

#include "cubature.h"
#include "LHAPDF/LHAPDF.h"

using namespace std;
using namespace LHAPDF;



///GLOBALS///



const unsigned g_dim = 3;                        //Dimension fo the integrals
const string g_pdfsetname = "HERAPDF20_LO_EIG";  //Name of the used pdf-set from LHAPDF
const int g_pdfsetmember = 1;                    //Member ID of the pdf in the set
const double g_s = 10;                           //Mandelstam variable s
const double g_kt_upper_cutoff = 10;             //Upper cutoff kt used in integration
const double g_kt_lower_cutoff = 1;              //Lower cutoff kt used in integration
const double g_error_tolerance = 1e-4;           //Global error tolerance



///PROTOTYPES///



int integrand_function(unsigned, const double*, void*, unsigned, double*);
int phasespace_integral(const double[], const double[], double * const, double * const);
int calculate_sigma_jet(double * const, double * const);
double f_ses(const double * const, const double * const, const PDF*);
complex<double> s_hat(const double * const, const double * const);
complex<double> t_hat(const double * const, const double * const, const double * const);
complex<double> u_hat(const double * const, const double * const, const double * const);
complex<double> sigma_gg_gg(const complex<double> * const, const complex<double> * const, const complex<double> * const);



///MAIN///



int main()
{
    int not_success;
    double sigma_jet, error;

    not_success = calculate_sigma_jet(&sigma_jet, &error);

    if(not_success)
    {
        cout<<"Bad stuff happened"<<endl;
        return 1;
    }
    else
    {
        printf("Sigma_jet = %0.10g +/- %g\n", sigma_jet, error);
        return 0;
    }
}



///FUNCTIONS///



///calculate_sigma_jet
///Function that does the calculation of the sigma_jet. Calls phasespace_integral.
///
///\param p_value = pointer to the destination of the value of the sigma_jet
///\param p_error = pointer to the destination of the error output from cubature
///
///\return 0 on successful run
///
int calculate_sigma_jet(double * const p_value, double * const p_error)
{
    int not_success;
    const double upper_limits [3] = {g_kt_upper_cutoff, 1, 1};
    const double lower_limits [3] = {g_kt_lower_cutoff, 0, 0};

    not_success = phasespace_integral(upper_limits, lower_limits, p_value, p_error);
    return not_success;
}


///phasespace_integral
///Function that calculates the phasespace integral with cubature algorithm. Uses integrand_function.
///
///\param upper_limits = list of the upper limits of the integral. Order: kt2, x1, x2
///\param lower_limits = list of the lower limits of the integral. Order: kt2, x1, x2
///\param p_value = pointer to the destination of the value of the integral
///\param p_error = pointer to the destination of the error output from cubature
///
///\return 0 on successful run
///
int phasespace_integral(const double upper_limits [3], const double lower_limits [3], double * const p_value, double * const p_error)
{

    const unsigned fdim=1;
    PDF* p_pdf = mkPDF(g_pdfsetname, g_pdfsetmember);

    int not_success = hcubature(fdim,               //Integrand dimension
                                integrand_function, //Integrand function
                                p_pdf,              //Pointer to additional arguments
                                g_dim,              //Variable dimension
                                lower_limits,       //Variables minimum
                                upper_limits,       //Variables maximum
                                0,                  //Max n:o of function evaluations
                                0,                  //Required absolute error
                                g_error_tolerance,  //Required relative error
                                ERROR_INDIVIDUAL,   //Enumerate of which norm is used on errors
                                p_value,            //Pointer to output
                                p_error);           //Pointer to error output

    if(not_success)
    {
        cout<<"Problem with integration"<<endl;
        return 1;
    }
    else
    {
        printf("Computed integral = %0.10g +/- %g\n", *p_value, *p_error);
        return 0;
    }
}


///integrand_function
///Function to be integrated by cubature algorithm. Calculates the value of the integrand in a
///phasespace point according to single effective subprocess approximation.
///
///\param ndim = dimension of the integral variable
///\param p_x = pointer to the integral variable
///\param p_fdata = pointer to the additional arguments (constants in integration)
///\param fdim = dimension of the output
///\param p_fval = pointer to the output destination
///
///\return 0 on successful run
///
int integrand_function(unsigned ndim, const double *p_x, void *p_fdata, unsigned fdim, double *p_fval)
{

    const PDF* p_pdf = (PDF *) p_fdata;
    const auto kt2 = p_x[0], x1 = p_x[1], x2 = p_x[2];
    const auto s = s_hat(&x1, &x2);
    const auto t = t_hat(&x1, &x2, &kt2);
    const auto u = u_hat(&x1, &x2, &kt2);
    const auto subprocess_cs = sigma_gg_gg(&s, &t, &u);

    const complex<double> dummy = x1*x2*(x1*x2 -4*(kt2/g_s));
    const auto jacobian = - abs(pow(sqrt(dummy),-1));

    p_fval[0] = 0.5 * x1 * f_ses(&x1, &kt2, p_pdf)
                * x2 * f_ses(&x2, &kt2, p_pdf)
                * real(subprocess_cs * jacobian);

    //cout<<"kt2="<<kt2<<", x1="<<x1<<", x2="<<x2<<", s="<<s<<", t="<<t<<", u="<<u<<", cs="<<subprocess_cs<<", J="<<jacobian<<", I="<<p_fval[0]<<endl;
    return 0; // success
}


///s_hat
///Calculates and returns the subprocess mandelstam variable s. Return value can be complex.
///
///\param p_x1 = pointer to momentum fraction of the first particle
///\param p_x2 = pointer to momentum fraction of the other particle
///
complex<double> s_hat(const double * const p_x1, const double * const p_x2)
{
    return *p_x1 * *p_x2 * g_s;
}


///t_hat
///Calculates and returns the subprocess mandelstam variable t. Return value can be complex.
///
///\param p_x1 = pointer to momentum fraction of the first particle
///\param p_x2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///
complex<double> t_hat(const double * const p_x1, const double * const p_x2, const double * const p_kt2)
{
    const double x1 = *p_x1, x2 = *p_x2, kt2 = *p_kt2;
    const complex<double> dummy = x1*x2*(x1*x2 -4*(kt2/g_s));
    return 0.5*g_s*(sqrt(dummy) -x1*x2);
}


///u_hat
///Calculates and returns the subprocess mandelstam variable u. Return value can be complex.
///
///\param p_x1 = pointer to momentum fraction of the first particle
///\param p_x2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///
complex<double> u_hat(const double * const p_x1, const double * const p_x2, const double * const p_kt2)
{
    const double x1 = *p_x1, x2 = *p_x2, kt2 = *p_kt2;
    const complex<double> dummy = x1*x2*(x1*x2 -4*(kt2/g_s));
    return -0.5*g_s*(sqrt(dummy) +x1*x2);
}


///f_ses
///Calculates and returns the parton distribution function F_ses according to single effective
///subprocess approximation.
///
///\param p_x = pointer to momentum fraction of the particle in question
///\param p_q2 = pointer to momentum exchange squared
///\param p_pdf = pointer to a member of the PDF class from LHAPDF
///
double f_ses(const double * const p_x, const double * const p_q2, const PDF* p_pdf)
{
    double sum = 0;
    for (int flavour=-6; flavour<=6; ++flavour)
    {
        sum += (flavour == 0) ? p_pdf->xfxQ2(flavour, *p_x, *p_q2) : 4 * p_pdf->xfxQ2(flavour, *p_x, *p_q2) / 9;
    }
    //cout<<"x="<<*p_x<<", q2="<<*p_q2<<" SES="<<sum<<endl;
    return sum;
}


///sigma_gg_gg
///Calculates and returns the cross section of the subprocess gg->gg. Return value should be real.
///
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///
complex<double> sigma_gg_gg(const complex<double> * const p_s_hat, const complex<double> * const p_t_hat, const complex<double> * const p_u_hat)
{
    const complex<double> s = *p_s_hat, t = *p_t_hat, u = *p_u_hat, constant = 3;
    return 4.5*(constant-(u*t)/(s*s)-(u*s)/(t*t)-(s*t)/(u*u));
}
