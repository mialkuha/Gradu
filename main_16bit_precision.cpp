#include <iostream>
#include <fstream>
#include <stdio.h>

#include <complex>
#include <gsl/gsl_math.h>

#include "cubature.h"
#include "LHAPDF/LHAPDF.h"

using namespace std;
using namespace LHAPDF;



///GLOBALS///



const unsigned g_dim = 3;                        //Dimension fo the integrals
const string g_pdfsetname = "CT14lo";  //Name of the used pdf-set from LHAPDF
const int g_pdfsetmember = 0;                   //Member ID of the pdf in the set
const double g_kt2_lower_cutoff = 4;              //Lower cutoff kt used in integration
const double g_error_tolerance = 1e-4;           //Global error tolerance



///PROTOTYPES///


int integrand_function_ys(unsigned, const double*, void*, unsigned, double*);
int integrand_function_xs(unsigned, const double*, void*, unsigned, double*);
int phasespace_integral(const double[], const double[], double * const, double * const, int, const long double * const);
int calculate_sigma_jet(double * const, double * const, int, const long double * const);
long double f_ses(const long double * const, const long double * const, const PDF*);
long double diff_sigma_jet(const long double * const, const long double * const, const long double * const, const PDF*, const long double * const, const long double * const, const long double * const);
long double s_hat_from_ys(const long double * const, const long double * const, const long double * const);
long double t_hat_from_ys(const long double * const, const long double * const, const long double * const);
long double u_hat_from_ys(const long double * const, const long double * const, const long double * const);
long double s_hat_from_xs(const long double * const, const long double * const, const long double * const);
long double t_hat_from_xs(const long double * const, const long double * const, const long double * const, const long double * const);
long double u_hat_from_xs(const long double * const, const long double * const, const long double * const, const long double * const);
long double sigma_qiqj_qiqj(const long double * const, const long double * const, const long double * const, const long double);
long double sigma_qiqi_qiqi(const long double * const, const long double * const, const long double * const, const long double);
long double sigma_qiaqi_qjaqj(const long double * const, const long double * const, const long double * const, const long double);
long double sigma_qiaqi_qiaqi(const long double * const, const long double * const, const long double * const, const long double);
long double sigma_qiaqi_gg(const long double * const, const long double * const, const long double * const, const long double);
long double sigma_gg_qaq(const long double * const, const long double * const, const long double * const, const long double);
long double sigma_gq_gq(const long double * const, const long double * const, const long double * const, const long double);
long double sigma_gg_gg(const long double * const, const long double * const, const long double * const, const long double);



///MAIN///



int main()
{
    int not_success_xs, not_success_ys;
    double sigma_jet_xs, error_xs;
    double sigma_jet_ys, error_ys;
    long double mand_s;

    ofstream xs_data_sigma_jet, ys_data_sigma_jet;

    xs_data_sigma_jet.open ("xs_sigma_jet.dat");
    ys_data_sigma_jet.open ("ys_sigma_jet.dat");
    xs_data_sigma_jet << "# sigma_jet calculated in x-space, p_0=" << sqrt(g_kt2_lower_cutoff) << '\n';
    xs_data_sigma_jet << "# sqrt(s) sigma_jet error" << '\n';
    ys_data_sigma_jet << "# sigma_jet calculated in y-space, p_0=" << sqrt(g_kt2_lower_cutoff) << '\n';
    ys_data_sigma_jet << "# sqrt(s) sigma_jet error" << '\n';

    for(long double sqrt_mand_s = 10; sqrt_mand_s <= 5000; sqrt_mand_s *=1.5){
        mand_s = sqrt_mand_s*sqrt_mand_s;

        not_success_xs = calculate_sigma_jet(&sigma_jet_xs, &error_xs, 0, &mand_s);
        not_success_ys = calculate_sigma_jet(&sigma_jet_ys, &error_ys, 1, &mand_s);

        xs_data_sigma_jet << sqrt_mand_s << ' ' << sigma_jet_xs << ' ' << error_xs << '\n';
        ys_data_sigma_jet << sqrt_mand_s << ' ' << sigma_jet_ys << ' ' << error_ys << '\n';
    }

    xs_data_sigma_jet.close();
    ys_data_sigma_jet.close();

    if(not_success_xs||not_success_ys){

        cout<<"Bad stuff happened"<<endl;

        return 1;
    }
    else{
        return 0;
    }
}



///FUNCTIONS///



///calculate_sigma_jet
///Function that does the calculation of the sigma_jet. Calls phasespace_integral.
///
///\param p_value = pointer to the destination of the value of the sigma_jet
///\param p_error = pointer to the destination of the error output from cubature
///\param y_space = flag, 0 = integrate in x-space, 1 = integrate in y-space
///\param p_mand_s = pointer to the mandelstam variable s
///
///\return 0 on successful run
///
int calculate_sigma_jet(double * const p_value, double * const p_error, int y_space, const long double * const p_mand_s)
{
    int not_success;
    const double upper_limits [3] = {1, 1, 1};
    const double lower_limits [3] = {0, 0, 0};

    not_success = phasespace_integral(upper_limits, lower_limits, p_value, p_error, y_space, p_mand_s);
    return not_success;
}


///phasespace_integral
///Function that calculates the phasespace integral with cubature algorithm. Uses integrand_function.
///
///\param upper_limits = list of the upper limits of the integral. Order: kt2, x1, x2
///\param lower_limits = list of the lower limits of the integral. Order: kt2, x1, x2
///\param p_value = pointer to the destination of the value of the integral
///\param p_error = pointer to the destination of the error output from cubature
///\param y_space = flag, 0 = integrate in x-space, 1 = integrate in y-space
///\param p_mand_s = pointer to the mandelstam variable s
///
///\return 0 on successful run
///
int phasespace_integral(const double upper_limits [3], const double lower_limits [3], double * const p_value, double * const p_error, int y_space, const long double * const p_mand_s)
{

    const unsigned fdim=1;
    PDF* p_pdf = mkPDF(g_pdfsetname, g_pdfsetmember);
    pair<PDF*,const long double * const> fdata = { p_pdf , p_mand_s };

    int not_success;

    if(y_space){
    not_success = hcubature(fdim,               //Integrand dimension
                                integrand_function_ys, //Integrand function
                                &fdata,              //Pointer to additional arguments
                                g_dim,              //Variable dimension
                                lower_limits,       //Variables minimum
                                upper_limits,       //Variables maximum
                                0,                  //Max n:o of function evaluations
                                0,                  //Required absolute error
                                g_error_tolerance,  //Required relative error
                                ERROR_INDIVIDUAL,   //Enumerate of which norm is used on errors
                                p_value,            //Pointer to output
                                p_error);           //Pointer to error output
    }
    else{
    not_success = hcubature(fdim,               //Integrand dimension
                                integrand_function_xs, //Integrand function
                                &fdata,              //Pointer to additional arguments
                                g_dim,              //Variable dimension
                                lower_limits,       //Variables minimum
                                upper_limits,       //Variables maximum
                                0,                  //Max n:o of function evaluations
                                0,                  //Required absolute error
                                g_error_tolerance,  //Required relative error
                                ERROR_INDIVIDUAL,   //Enumerate of which norm is used on errors
                                p_value,            //Pointer to output
                                p_error);           //Pointer to error output
    }

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

///integrand_function_ys
///Function to be integrated by cubature algorithm. Calculates the value of the integrand in a
///phasespace point according to single effective subprocess approximation. TODO
///
///\param ndim = dimension of the integral variable
///\param p_x = pointer to the integral variable
///\param p_fdata = pointer to the additional arguments (constants in integration)
///\param fdim = dimension of the output
///\param p_fval = pointer to the output destination
///
///\return 0 on successful run
///
int integrand_function_ys(unsigned ndim, const double *p_x, void *p_fdata, unsigned fdim, double *p_fval)
{
    pair<PDF*,const long double * const> fdata = *(pair<PDF*,const long double * const>*) p_fdata;
    const PDF* p_pdf = fdata.first;
    const long double * const p_mand_s = fdata.second;
    const auto z1 = p_x[0], z2 = p_x[1], z3 = p_x[2];

    const auto kt2           = (long double) g_kt2_lower_cutoff + z1 * ((*p_mand_s/4.0L) - (long double) g_kt2_lower_cutoff);
    const auto sqrt_s_per_kt = sqrtl(*p_mand_s/ (long double) kt2);

    const auto y1_upper = acoshl(sqrt_s_per_kt/2.0L);
    //const auto y1_lower = -y1_upper;
    const auto y1       = ( -1.0L + (2.0L * z2) ) * y1_upper; //y1_lower + z2 * (y1_upper - y1_lower)

    const auto y2_upper = logl(sqrt_s_per_kt - expl(y1));
    const auto y2_lower = -logl(sqrt_s_per_kt - expl(-y1));
    const auto y2       = y2_lower + z3 * (y2_upper - y2_lower);

    const auto x1 = (expl(y1) + expl(y2)) / sqrt_s_per_kt;
    const auto x2 = (expl(-y1) + expl(-y2)) / sqrt_s_per_kt;

    const auto s_hat         = s_hat_from_ys(&y1, &y2, &kt2);
    const auto t_hat         = t_hat_from_ys(&y1, &y2, &kt2);
    const auto u_hat         = u_hat_from_ys(&y1, &y2, &kt2);
    const auto subprocess_cs = sigma_gg_gg(&s_hat, &t_hat, &u_hat, p_pdf->alphasQ2(kt2));

    const auto jacobian = ((*p_mand_s/4.0L) - g_kt2_lower_cutoff) * (2.0L*y1_upper) * (y2_upper - y2_lower);

    if (y2_upper < y2_lower || y1_upper < -y1_upper) p_fval[0]=0;
    else{
    //SES
/*
    p_fval[0] = 0.5 * f_ses(&x1, &kt2, p_pdf)
                * f_ses(&x2, &kt2, p_pdf)
                * subprocess_cs * jacobian;
*/
    //FULL SUMMATION
    p_fval[0] = diff_sigma_jet(&x1,&x2,&kt2,p_pdf,&s_hat,&t_hat,&u_hat)
                * jacobian;
    }


    //cout <<subprocess_cs<<' ';
/*
    ofstream data;
    data.open ("ys_mandelstams.dat",ios::app);
    data << s_hat << ' ' << t_hat << ' ' << u_hat << '\n';
    data.close();

    data.open ("ys_kt2x1x2.dat",ios::app);
    data << kt2 << ' ' << x1 << ' ' << x2 << '\n';
    data.close();

    data.open ("ys_kt2y1y2.dat",ios::app);
    data << kt2 << ' ' << y1 << ' ' << y2 << '\n';
    data.close();

    data.open ("ys_integrand.dat",ios::app);
    data << kt2 << ' ' << x1 << ' ' << x2 << ' ' << p_fval[0] << '\n';
    data.close();
*/
    //cout<<"kt2="<<kt2<<", x1="<<x1<<", x2="<<x2<<", s="<<s<<", t="<<t<<", u="<<u<<", cs="<<subprocess_cs<<", J="<<jacobian<<", I="<<p_fval[0]<<endl;

    return 0; // success
}



///integrand_function_xs
///Function to be integrated by cubature algorithm. Calculates the value of the integrand in a
///phasespace point according to single effective subprocess approximation. TODO
///
///\param ndim = dimension of the integral variable
///\param p_x = pointer to the integral variable
///\param p_fdata = pointer to the additional arguments (constants in integration)
///\param fdim = dimension of the output
///\param p_fval = pointer to the output destination
///
///\return 0 on successful run
///
int integrand_function_xs(unsigned ndim, const double *p_x, void *p_fdata, unsigned fdim, double *p_fval)
{

    pair<PDF*,const long double * const> fdata = *(pair<PDF*,const long double * const>*) p_fdata;
    const PDF* p_pdf = fdata.first;
    const long double * const p_mand_s = fdata.second;
    const auto z1 = p_x[0], z2 = p_x[1], z3 = p_x[2];

    const auto kt2 = (long double) g_kt2_lower_cutoff + z1 * ((*p_mand_s/4.0L) - (long double) g_kt2_lower_cutoff);

    const auto x1_upper = 1.0L;
    const auto x1_lower = 4.0L * kt2 / *p_mand_s;
    const auto x1       = x1_lower + z2 * (x1_upper - x1_lower);

    const auto x2_upper = 1.0L;
    const auto x2_lower = (4.0L * kt2 / *p_mand_s) / x1 +1e-16L; //shift by epsilon to avoid singularity in jacobian_from_ys_to_xs
    const auto x2       = x2_lower + z3 * (x2_upper - x2_lower);

    const auto s_hat         = s_hat_from_xs(&x1, &x2, p_mand_s);
    const auto t_hat         = t_hat_from_xs(&x1, &x2, &kt2, p_mand_s);
    const auto u_hat         = u_hat_from_xs(&x1, &x2, &kt2, p_mand_s);
    const auto subprocess_cs = sigma_gg_gg(&s_hat, &t_hat, &u_hat, p_pdf->alphasQ2(kt2));

    const auto jacobian = ((*p_mand_s/4) - g_kt2_lower_cutoff) * (x1_upper - x1_lower) * (x2_upper - x2_lower);

    auto dummy = x1*x2-4.0L*kt2/ *p_mand_s;
    const auto dummy2 = x1*x2*dummy;
    const auto jacobian_from_ys_to_xs = 1.0L/sqrt(dummy2);

    if (x2_upper < x2_lower || x1_upper < x1_lower) p_fval[0]=0;
    else{
    //SES
/*
    p_fval[0] = 2  //MAGIC NUMBER TODO
                * 0.5 * f_ses(&x1, &kt2, p_pdf)
                * f_ses(&x2, &kt2, p_pdf)
                * subprocess_cs * jacobian
                * jacobian_from_ys_to_xs;
*/
    //FULL SUMMATION
    p_fval[0] = 2  //MAGIC NUMBER TODO
                * diff_sigma_jet(&x1,&x2,&kt2,p_pdf,&s_hat,&t_hat,&u_hat)
                * jacobian * jacobian_from_ys_to_xs;

    }



    //cout <<subprocess_cs<<' ';

/*
    ofstream data;
    data.open ("xs_mandelstams.dat",ios::app);
    data << s_hat << ' ' << t_hat << ' ' << u_hat << '\n';
    data.close();

    auto y1 = log((x1+sqrt(x1*x2*(-4*sqrt(kt2/ *p_mand_s)*sqrt(kt2/ *p_mand_s)+x1*x2))/x2)/(2*sqrt(kt2/ *p_mand_s)));
    auto y2 = log((x1-sqrt(x1*x2*(-4*sqrt(kt2/ *p_mand_s)*sqrt(kt2/ *p_mand_s)+x1*x2))/x2)/(2*sqrt(kt2/ *p_mand_s)));
    data.open ("xs_kt2y1y2.dat",ios::app);
    data << kt2 << ' ' << y1 << ' ' << y2 << '\n';
    data.close();

    data.open ("xs_all.dat",ios::app);
    data<<dummy<<dummy2<<"kt2="<<kt2<<", x1="<<x1<<", x2="<<x2<<", s="<<s_hat<<", t="<<t_hat<<", u="<<u_hat<<", cs="<<subprocess_cs<<", J="<<jacobian<<", J2="<<jacobian_from_ys_to_xs<<", pdf1="<<f_ses(&x1, &kt2, p_pdf)<<", pdf2="<<f_ses(&x2, &kt2, p_pdf)<<", I="<<p_fval[0]<<'\n';
    data.close();

    data.open ("xs_kt2x1x2.dat",ios::app);
    data << kt2 << ' ' << x1 << ' ' << x2 << '\n';
    data.close();

    data.open ("xs_integrand.dat",ios::app);
    data << kt2 << ' ' << x1 << ' ' << x2 << ' ' << p_fval[0] << '\n';
    data.close();
*/
    //cout<<"kt2="<<kt2<<", x1="<<x1<<", x2="<<x2<<", s="<<s_hat<<", t="<<t_hat<<", u="<<u_hat<<", cs="<<subprocess_cs<<", J="<<jacobian<<", J2="<<jacobian_from_ys_to_xs<<", pdf1="<<f_ses(&x1, &kt2, p_pdf)<<", pdf2="<<f_ses(&x2, &kt2, p_pdf)<<", I="<<p_fval[0]<<endl;

    return 0; // success
}


///s_hat_from_ys
///Calculates and returns the subprocess mandelstam variable s. Return value can be complex.
///
///\param p_y1 = pointer to momentum fraction of the first particle
///\param p_y2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///
long double s_hat_from_ys(const long double * const p_y1, const long double * const p_y2, const long double * const p_kt2)
{
    const long double y1 = *p_y1, y2 = *p_y2, kt2 = *p_kt2;
    return 2*kt2 * (1 + cosh(y1 - y2));
}


///t_hat_from_ys
///Calculates and returns the subprocess mandelstam variable t. Return value can be complex.
///
///\param p_y1 = pointer to momentum fraction of the first particle
///\param p_y2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///
long double t_hat_from_ys(const long double * const p_y1, const long double * const p_y2, const long double * const p_kt2)
{
    const long double y1 = *p_y1, y2 = *p_y2, kt2 = *p_kt2;
    return -kt2 * (1 + exp(-y1 + y2));
}


///u_hat_from_ys
///Calculates and returns the subprocess mandelstam variable u. Return value can be complex.
///
///\param p_y1 = pointer to momentum fraction of the first particle
///\param p_y2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///
long double u_hat_from_ys(const long double * const p_y1, const long double * const p_y2, const long double * const p_kt2)
{
    const long double y1 = *p_y1, y2 = *p_y2, kt2 = *p_kt2;
    return -kt2 * (1 + exp(y1 - y2));
}


///s_hat_from_xs
///Calculates and returns the subprocess mandelstam variable s. Return value can be complex.
///
///\param p_x1 = pointer to momentum fraction of the first particle
///\param p_x2 = pointer to momentum fraction of the other particle
///\param p_mand_s = pointer to mandelstam variable s
///
long double s_hat_from_xs(const long double * const p_x1, const long double * const p_x2, const long double * const p_mand_s)
{
    return *p_x1 * *p_x2 * *p_mand_s;
}


///t_hat_from_xs
///Calculates and returns the subprocess mandelstam variable t. Return value can be complex.
///
///\param p_x1 = pointer to momentum fraction of the first particle
///\param p_x2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///\param p_mand_s = pointer to mandelstam variable s
///
long double t_hat_from_xs(const long double * const p_x1, const long double * const p_x2, const long double * const p_kt2, const long double * const p_mand_s)
{
    const long double x1 = *p_x1, x2 = *p_x2, kt2 = *p_kt2;
    const long double dummy = x1*x2*(x1*x2 -4*(kt2/ *p_mand_s));
    return 0.5* *p_mand_s *(sqrt(dummy) -x1*x2);
}


///u_hat_from_xs
///Calculates and returns the subprocess mandelstam variable u. Return value can be complex.
///
///\param p_x1 = pointer to momentum fraction of the first particle
///\param p_x2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///\param p_mand_s = pointer to mandelstam variable s
///
long double u_hat_from_xs(const long double * const p_x1, const long double * const p_x2, const long double * const p_kt2, const long double * const p_mand_s)
{
    const long double x1 = *p_x1, x2 = *p_x2, kt2 = *p_kt2;
    const long double dummy = x1*x2*(x1*x2 -4*(kt2/ *p_mand_s));
    return -0.5* *p_mand_s *(sqrt(dummy) +x1*x2);
}


///f_ses
///Calculates and returns the parton distribution function F_ses according to single effective
///subprocess approximation.
///
///\param p_x = pointer to momentum fraction of the particle in question
///\param p_q2 = pointer to momentum exchange squared
///\param p_pdf = pointer to a member of the PDF class from LHAPDF
///
long double f_ses(const long double * const p_x, const long double * const p_q2, const PDF* p_pdf)
{
    long double sum = 0;
    for (int flavour=-6; flavour<=6; ++flavour)
    {
        sum += (flavour == 0) ? p_pdf->xfxQ2(flavour, *p_x, *p_q2) : 4 * p_pdf->xfxQ2(flavour, *p_x, *p_q2) / 9;
    }
    //cout<<"x="<<*p_x<<", q2="<<*p_q2<<" SES="<<sum<<endl;
    return sum;
}


///diff_sigma_jet
///Calculates and returns the differential cross section.
///
///\param p_x1 = pointer to momentum fraction of the first particle in question
///\param p_x2 = pointer to momentum fraction of the first particle in question
///\param p_q2 = pointer to momentum exchange squared
///\param p_pdf = pointer to a member of the PDF class from LHAPDF
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///
long double diff_sigma_jet(const long double * const p_x1, const long double * const p_x2, const long double * const p_q2, const PDF* p_pdf, const long double * const p_s_hat, const long double * const p_t_hat, const long double * const p_u_hat)
{
    const long double x1 = *p_x1, x2 = *p_x2, q2 = *p_q2;
    long double sum = 0;

    sum += 0.5 * p_pdf->xfxQ2(0, x1, q2) * p_pdf->xfxQ2(0, x2, q2) * sigma_gg_gg(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2))
           + 6 * p_pdf->xfxQ2(0, x1, q2) * p_pdf->xfxQ2(0, x2, q2) * sigma_gg_qaq(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));

    for (int flavour=1; flavour<=6; ++flavour)
    {
        sum += p_pdf->xfxQ2(0, x1, q2) * p_pdf->xfxQ2(flavour, x2, q2) * sigma_gq_gq(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
        sum += p_pdf->xfxQ2(0, x1, q2) * p_pdf->xfxQ2(-flavour, x2, q2) * sigma_gq_gq(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
        sum += p_pdf->xfxQ2(flavour, x1, q2) * p_pdf->xfxQ2(0, x2, q2) * sigma_gq_gq(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
        sum += p_pdf->xfxQ2(-flavour, x1, q2) * p_pdf->xfxQ2(0, x2, q2) * sigma_gq_gq(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
    }

    for (int flavour=1; flavour<=6; ++flavour)
    {
        sum += 0.5 * p_pdf->xfxQ2(flavour, x1, q2) * p_pdf->xfxQ2(flavour, x2, q2) * sigma_qiqi_qiqi(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
        sum += 0.5 * p_pdf->xfxQ2(-flavour, x1, q2) * p_pdf->xfxQ2(-flavour, x2, q2) * sigma_qiqi_qiqi(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
    }

    for (int flavour1=1; flavour1<=6; ++flavour1)
    {
    for (int flavour2=1; flavour2<=6; ++flavour2)
    {
    if (flavour1 != flavour2){
        sum += p_pdf->xfxQ2(flavour1, x1, q2) * p_pdf->xfxQ2(flavour2, x2, q2) * sigma_qiqj_qiqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
        sum += p_pdf->xfxQ2(-flavour1, x1, q2) * p_pdf->xfxQ2(-flavour2, x2, q2) * sigma_qiqj_qiqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
        sum += p_pdf->xfxQ2(flavour1, x1, q2) * p_pdf->xfxQ2(-flavour2, x2, q2) * sigma_qiqj_qiqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
        sum += p_pdf->xfxQ2(-flavour1, x1, q2) * p_pdf->xfxQ2(flavour2, x2, q2) * sigma_qiqj_qiqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
    }}}

    for (int flavour=1; flavour<=6; ++flavour)
    {
        sum += p_pdf->xfxQ2(flavour, x1, q2) * p_pdf->xfxQ2(-flavour, x2, q2)
             *(sigma_qiaqi_qiaqi(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2))
             + 0.5*sigma_qiaqi_gg(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2))
             + 5*sigma_qiaqi_qjaqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2)));
        sum += p_pdf->xfxQ2(-flavour, x1, q2) * p_pdf->xfxQ2(flavour, x2, q2)
             *(sigma_qiaqi_qiaqi(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2))
             + 0.5*sigma_qiaqi_gg(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2))
             + 5*sigma_qiaqi_qjaqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2)));
    }

    return sum;
}


///sigma_qiqj_qiqj
///Calculates and returns the cross section of the subprocess qiqj->qiqj. Return value should be real.
///
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///\param p_alpha_s = pointer to strong interaction constant
///
long double sigma_qiqj_qiqj(const long double * const p_s_hat, const long double * const p_t_hat, const long double * const p_u_hat, const long double alpha_s)
{
    const long double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s))*4*((s*s+u*u)/(t*t))/9;
}


///sigma_qiqi_qiqi
///Calculates and returns the cross section of the subprocess qiqi->qiqi. Return value should be real.
///
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///\param p_alpha_s = pointer to strong interaction constant
///
long double sigma_qiqi_qiqi(const long double * const p_s_hat, const long double * const p_t_hat, const long double * const p_u_hat, const long double alpha_s)
{
    const long double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s))*4*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u)-(2*s*s)/(3*t*u))/9;
}


///sigma_qiaqi_qjaqj
///Calculates and returns the cross section of the subprocess qiaqi->qjaqj. Return value should be real.
///
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///\param p_alpha_s = pointer to strong interaction constant
///
long double sigma_qiaqi_qjaqj(const long double * const p_s_hat, const long double * const p_t_hat, const long double * const p_u_hat, const long double alpha_s)
{
    const long double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s))*4*((t*t+u*u)/(s*s))/9;
}


///sigma_qiaqi_qiaqi
///Calculates and returns the cross section of the subprocess qiaqi->qiaqi. Return value should be real.
///
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///\param p_alpha_s = pointer to strong interaction constant
///
long double sigma_qiaqi_qiaqi(const long double * const p_s_hat, const long double * const p_t_hat, const long double * const p_u_hat, const long double alpha_s)
{
    const long double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s))*4*((s*s+u*u)/(t*t)+(t*t+u*u)/(s*s)-(2*u*u)/(3*s*t))/9;
}


///sigma_qiaqi_gg
///Calculates and returns the cross section of the subprocess qiaqi->gg. Return value should be real.
///
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///\param p_alpha_s = pointer to strong interaction constant
///
long double sigma_qiaqi_gg(const long double * const p_s_hat, const long double * const p_t_hat, const long double * const p_u_hat, const long double alpha_s)
{
    const long double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s))*8*(t*t+u*u)*(4/(9*t*u)-1/(s*s))/3;
}


///sigma_gg_qaq
///Calculates and returns the cross section of the subprocess gg->qaq. Return value should be real.
///
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///\param p_alpha_s = pointer to strong interaction constant
///
long double sigma_gg_qaq(const long double * const p_s_hat, const long double * const p_t_hat, const long double * const p_u_hat, const long double alpha_s)
{
    const long double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s))*3*(t*t+u*u)*(4/(9*t*u)-1/(s*s))/8;
}


///sigma_gq_gq
///Calculates and returns the cross section of the subprocess gg->gg. Return value should be real.
///
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///\param p_alpha_s = pointer to strong interaction constant
///
long double sigma_gq_gq(const long double * const p_s_hat, const long double * const p_t_hat, const long double * const p_u_hat, const long double alpha_s)
{
    const long double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s))*(s*s+u*u)*(1/(t*t)-4/(9*s*u));
}


///sigma_gg_gg
///Calculates and returns the cross section of the subprocess gg->gg. Return value should be real.
///
///\param p_s_hat = pointer to subprocess mandelstam variable s
///\param p_t_hat = pointer to subprocess mandelstam variable t
///\param p_u_hat = pointer to subprocess mandelstam variable u
///\param p_alpha_s = pointer to strong interaction constant
///
long double sigma_gg_gg(const long double * const p_s_hat, const long double * const p_t_hat, const long double * const p_u_hat, const long double alpha_s)
{
    const long double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s))*4.5*(3.0-(u*t)/(s*s)-(u*s)/(t*t)-(s*t)/(u*u));
}
