#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "cubature.h"
#include "LHAPDF/LHAPDF.h"

using namespace std;
using namespace LHAPDF;



///GLOBALS///



const unsigned g_dim = 3;                        //Dimension fo the integrals
const string g_pdfsetname = "CT14lo";  //Name of the used pdf-set from LHAPDF
const int g_pdfsetmember = 0;                   //Member ID of the pdf in the set
const double g_error_tolerance = 1e-4;           //Global error tolerance
const int g_eikonal_sum_term_count = 3;         //N:o of eikonal sum terms to study
const int g_data_point_count = 13;                //N:o of data point constants
const double g_sqrt_s_list [g_data_point_count] = {100, 200, 400, 550, 1000, 1800, 4000, 7000, 8000, 13000, 25000, 50000, 100000};


///PROTOTYPES///


int eikonal_integrand(unsigned, const double*, void*, unsigned, double*);
int integrand_function_ys(unsigned, const double*, void*, unsigned, double*);
double find_kt2_lower_cutoff(const double * const, const double * const, double *, double *, PDF*);
int phasespace_integral(const double[], const double[], double * const, double * const, const double * const, const double * const, PDF*);
int calculate_sigma_inel_from_sigma_jet(double * const, double * const);
int calculate_sigma_tot_from_sigma_jet(double * const, double * const);
int calculate_sigma_inel_from_integral(double * const, double[]);
int calculate_sigma_tot_from_integral(double * const, double[]);
int calculate_sigma_jet(double * const, double * const, const double * const, const double * const, PDF*);
void find_eikonal_sum_contributions(const double * const, double [], const double, double[], PDF*);
double f_ses(const double * const, const double * const, const PDF*);
double diff_sigma_jet(const double * const, const double * const, const double * const, const PDF*, const double * const, const double * const, const double * const);
double s_hat_from_ys(const double * const, const double * const, const double * const);
double t_hat_from_ys(const double * const, const double * const, const double * const);
double u_hat_from_ys(const double * const, const double * const, const double * const);
double sigma_qiqj_qiqj(const double * const, const double * const, const double * const, const double);
double sigma_qiqi_qiqi(const double * const, const double * const, const double * const, const double);
double sigma_qiaqi_qjaqj(const double * const, const double * const, const double * const, const double);
double sigma_qiaqi_qiaqi(const double * const, const double * const, const double * const, const double);
double sigma_qiaqi_gg(const double * const, const double * const, const double * const, const double);
double sigma_gg_qaq(const double * const, const double * const, const double * const, const double);
double sigma_gq_gq(const double * const, const double * const, const double * const, const double);
double sigma_gg_gg(const double * const, const double * const, const double * const, const double);
double sigma_jet_product_integrand(double *, size_t, void *);
int calculate_sigma_jet_product(double *, double *, const double * const, const double * const, PDF*, int);


///MAIN///



int main()
{
    double kt2_lower_cutoff [g_data_point_count];
    double eikonal_sum_contributions [g_data_point_count][g_eikonal_sum_term_count];
    double data_sigma_inel;
    double sigma_jet[g_eikonal_sum_term_count+1], sigma_inel, sigma_el, sigma_tot;
    ofstream data, histogram;

    PDF* p_pdf = mkPDF(g_pdfsetname, g_pdfsetmember);

    histogram.open("histogram_data.dat");

    for(int i=0; i<g_data_point_count; i++)
    {
        data_sigma_inel = 42.6*pow(g_sqrt_s_list[i],-0.92)
                         -33.4*pow(g_sqrt_s_list[i],-1.09)
                         +0.307*pow(log(pow(g_sqrt_s_list[i],2)/29.1),2)
                         +35.5 -(//FIT BY COMPETE FOR SIGMA_TOT
                         11.7 -1.59*log(pow(g_sqrt_s_list[i],2))
                         +0.134*pow(log(pow(g_sqrt_s_list[i],2)),2));//FIT BY TOTEM FOR SIGMA_EL

        ostringstream fn;
        fn << "data_sqrt(s)=" << g_sqrt_s_list[i] << ".dat";

        data.open(fn.str().c_str());

        kt2_lower_cutoff[i] = find_kt2_lower_cutoff(&g_sqrt_s_list[i], &data_sigma_inel, sigma_jet, &sigma_inel, p_pdf);

        data << "#kt² lower cutoff: " << kt2_lower_cutoff[i] << '\n';
        data << "#data_sqrt(s): " << g_sqrt_s_list[i] << '\n';


        calculate_sigma_tot_from_integral(&sigma_tot, sigma_jet);

        find_eikonal_sum_contributions(kt2_lower_cutoff+i, eikonal_sum_contributions[i], sigma_inel, sigma_jet, p_pdf);

        sigma_el = sigma_tot - sigma_inel;

        double dummy = 0;
        histogram << g_sqrt_s_list[i] << ' ';

        data << "#sigma_tot: " << sigma_tot << '\n';
        data << "#sigma_el: " << sigma_el << '\n';
        data << "#sigma_inel: " << sigma_inel << '\n';
        data << "#sigma_jet: " << sigma_jet[0] << '\n';
        data << "#eikonal sum term contributions:" << '\n';
        for(int j=0; j<g_eikonal_sum_term_count; ++j)
        {
            data << j+1 << ". " << eikonal_sum_contributions[i][j] << '%' << '\n';
            dummy += eikonal_sum_contributions[i][j];
            histogram << dummy << ' ';
        }
        histogram << '\n';
        data.close();
    }

    histogram.close();

    return 0;
}



///FUNCTIONS///



///sigma_jet_product_integrand
///Integrand function for calculating products of sigma_jets with momentum conservation taken into account.
///
///\param k = pointer to the integration variable
///\param fulldim = dimension of the integration variable
///\param p_fdata = void pointer to the parameters of the function
///
double sigma_jet_product_integrand(double *p_x, size_t fulldim, void *p_fdata){

    pair<PDF*,pair<const double * const,const double * const> >  fdata = *(pair<PDF*,pair<const double * const,const double * const> > *) p_fdata;

    const PDF* p_pdf = fdata.first;
    const double * const p_mand_s = fdata.second.first;
    const double * const p_kt2_lower_cutoff = fdata.second.second;

    const int dim = fulldim/3;
    double z[3][dim];

    double kt2[dim];
    double sqrt_s_per_kt[dim];

    double y1_upper[dim];
    double y1[dim];

    double y2_upper[dim];
    double y2_lower[dim];
    double y2[dim];

    double x1[dim];
    double x2[dim];

    for(int i=0; i<dim; ++i){
        z[0][i] = p_x[3*i];
        z[1][i] = p_x[3*i +1];
        z[2][i] = p_x[3*i +2];

        kt2[i] = *p_kt2_lower_cutoff + z[0][i] * ((*p_mand_s/4) - *p_kt2_lower_cutoff);

        sqrt_s_per_kt[i] = sqrt(*p_mand_s/kt2[i]);

        y1_upper[i] = acosh(sqrt_s_per_kt[i]/2);
        //y1_lower = -y1_upper;
        y1[i] = ( -1 + (2 * z[1][i]) ) * y1_upper[i]; //y1_lower + z2 * (y1_upper - y1_lower)

        y2_upper[i] = log(sqrt_s_per_kt[i] - exp(y1[i]));
        y2_lower[i] = -log(sqrt_s_per_kt[i] - exp(-y1[i]));
        y2[i] = y2_lower[i] + z[2][i] * (y2_upper[i] - y2_lower[i]);

        x1[i] = (exp(y1[i]) + exp(y2[i])) / sqrt_s_per_kt[i];
        x2[i] = (exp(-y1[i]) + exp(-y2[i])) / sqrt_s_per_kt[i];
    }

    auto sum_x1s=0;
    auto sum_x2s=0;

    for (auto x:x1) sum_x1s+=x;
    for (auto x:x2) sum_x2s+=x;

    if (sum_x1s>1 || sum_x2s>1) return 0; //MOMENTUM CONSERVATION



    double s_hat[dim];
    double t_hat[dim];
    double u_hat[dim];
    //double subprocess_cs[dim];
    double jacobian[dim];
    double sigma[dim];

    for(int i=0; i<dim; ++i){
        s_hat[i] = s_hat_from_ys(&y1[i], &y2[i], &kt2[i]);
        t_hat[i] = t_hat_from_ys(&y1[i], &y2[i], &kt2[i]);
        u_hat[i] = u_hat_from_ys(&y1[i], &y2[i], &kt2[i]);

        //subprocess_cs[i] = sigma_gg_gg(&s_hat[i], &t_hat[i], &u_hat[i], p_pdf->alphasQ2(kt2[i]));

        jacobian[i] = ((*p_mand_s/4) - *p_kt2_lower_cutoff) * (2*y1_upper[i]) * (y2_upper[i] - y2_lower[i]);

        //SES
        /*
        sigma[i] = 0.5 * f_ses(&x1[i], &kt2[i], p_pdf)
                    * f_ses(&x2[i], &kt2[i], p_pdf)
                    * subprocess_cs[i] * jacobian[i]; */
        //FULL SUMMATION
        sigma[i] = 0.5 * diff_sigma_jet(&x1[i],&x2[i],&kt2[i],p_pdf,&s_hat[i],&t_hat[i],&u_hat[i])
                    * jacobian[i];
    }

    auto result = 1;
    for (auto s:sigma) result*=s;

    return result;
}


///find_eikonal_sum_contributions
///Function that finds the lower cutoff for kt2 that fits the data with secant method.
///
///\param p_kt2_lower_cutoff = pointer to the p_kt2_lower_cutoff
///\param eikonal_sum_contributions = array to be filled wtih eikonal_sum_contributions
///\param sigma_inel = calculated sigma_inel
///\param sigma_jet = calculated sigma_jet
///\param p_pdf = pointer to the PDF object
///
void find_eikonal_sum_contributions(const double * const p_kt2_lower_cutoff, double eikonal_sum_contributions [g_eikonal_sum_term_count], const double sigma_inel, double sigma_jet[g_eikonal_sum_term_count+1], PDF* p_pdf)
{
    double eikonal_sum_errors [g_eikonal_sum_term_count];
//    ofstream ys_data_sigma_jet;

    const double upper_limits = 1;
    const double lower_limits = 0;

    for(int i=0; i<g_eikonal_sum_term_count; i++)
    {

        pair< double, pair< double , int > > fdata = {sigma_jet[0],{ sigma_jet[i+1], i+1 }};

        hcubature(1,               //Integrand dimension
                  eikonal_integrand, //Integrand function
                  &fdata,              //Pointer to additional arguments
                  1,              //Variable dimension
                  &lower_limits,       //Variables minimum
                  &upper_limits,       //Variables maximum
                  0,                  //Max n:o of function evaluations
                  0,                  //Required absolute error
                  g_error_tolerance,  //Required relative error
                  ERROR_INDIVIDUAL,   //Enumerate of which norm is used on errors
                  eikonal_sum_contributions +i,            //Pointer to output
                  eikonal_sum_errors +i);           //Pointer to error output

        eikonal_sum_contributions[i] = 100 * eikonal_sum_contributions[i]/ sigma_inel;
    }
}


///find_kt2_lower_cutoff
///Function that finds the lower cutoff for kt2 that fits the data with secant method.
///
///\param p_sqrt_s = pointer to the sqrt_s
///\param p_data_sigma_inel = pointer to the sigma_inel from data
///\param p_sigma_jet_res = pointer to the destination of sigma jet
///\param p_sigma_inel_res = pointer to the destination of sigma inel
///\param p_pdf = pointer to the PDF object
///
///\return kt2 lower cutoff that fits the data
///
double find_kt2_lower_cutoff(const double * const p_sqrt_s, const double * const p_data_sigma_inel, double * p_sigma_jet_res, double * p_sigma_inel_res, PDF* p_pdf)
{
    int not_success_ys=0;
    double sigma_jet_ys [20][g_eikonal_sum_term_count+1], sigma_jet_error_ys [20][g_eikonal_sum_term_count+1];
    double sigma_inel_ys [20];
    double error [20];
    double mand_s , kt2_lower_cutoff [20];
//    ofstream ys_data_sigma_jet;

    mand_s = *p_sqrt_s **p_sqrt_s;

    kt2_lower_cutoff[0] = 2;
    kt2_lower_cutoff[1] = 4;

    for(int i=0; i<20; ++i)
    {
        if (i>1)
        {
            kt2_lower_cutoff[i] = (kt2_lower_cutoff[i-2] * error[i-1] - kt2_lower_cutoff[i-1] * error[i-2])
                                  / (error[i-1]-error[i-2]);
        }
        if (kt2_lower_cutoff[i] <= 0) kt2_lower_cutoff[i] = 0.5;

        /*
            ys_data_sigma_jet.open ("ys_sigma_jet.dat");
            ys_data_sigma_jet << "# sigma_jet calculated in y-space, p_0=" << sqrt(kt2_lower_cutoff) << '\n';
            ys_data_sigma_jet << "# sqrt(s) sigma_jet error" << '\n';*/

        //for(double sqrt_mand_s = 10; sqrt_mand_s <= 5000; sqrt_mand_s *=1.2){
        //    mand_s = sqrt_mand_s*sqrt_mand_s;

        not_success_ys = calculate_sigma_jet(&sigma_jet_ys[i][0], &sigma_jet_error_ys[i][0], &mand_s, kt2_lower_cutoff+i, p_pdf);
            cout<<0<<' '<<sigma_jet_ys[i][0]<<endl;

        for (int j=1; j<g_eikonal_sum_term_count+1; ++j){
            calculate_sigma_jet_product(&sigma_jet_ys[i][j], &sigma_jet_error_ys[i][j], &mand_s, kt2_lower_cutoff+i, p_pdf, j);
            cout<<j<<' '<<sigma_jet_ys[i][j]<<endl;
        }
        //    ys_data_sigma_jet << mand_s << ' ' << sigma_jet_ys << ' ' << sigma_jet_error_ys << '\n';
        //}

        //ys_data_sigma_jet.close();

        not_success_ys = calculate_sigma_inel_from_integral(sigma_inel_ys+i, sigma_jet_ys[i]);

        cout << kt2_lower_cutoff[i] << ' ' << sigma_inel_ys[i] << ' ' << *p_data_sigma_inel << endl;

        error[i] = (sigma_inel_ys[i] - *p_data_sigma_inel);

        if (not_success_ys) cout << "error at i="<<i<<" kt2_lower_cutoff=" <<kt2_lower_cutoff[i]<<endl;

        if ( abs(error[i] / *p_data_sigma_inel) <= 1e-3){
            *p_sigma_inel_res = sigma_inel_ys[i];
            *p_sigma_jet_res = sigma_jet_ys[i][0];
            return kt2_lower_cutoff[i];
        }
    }

    *p_sigma_inel_res = sigma_inel_ys[19];
    *p_sigma_jet_res = sigma_jet_ys[19][0];
    return kt2_lower_cutoff[19];
}


///calculate_sigma_inel_from_sigma_jet
///Function that does the calculation of the sigma_inel based on sigma_jet.
///
///\param p_value = pointer to the destination of the value of the sigma_inel
///\param p_sigma_jet = pointer to the value of sigma_jet
///
///\return 0 on successful run
///
int calculate_sigma_inel_from_sigma_jet(double * const p_value, double * const p_sigma_jet)
{
    int not_success = 0;

    double variation = 4.72;
    double dummy = *p_sigma_jet / (4 * M_PI * variation);

    *p_value = (4 * M_PI * variation) * (M_EULER + log(dummy) + gsl_sf_expint_E1(dummy)) ;

    return not_success;
}


///calculate_sigma_inel_from_integral
///Function that does the calculation of the sigma_inel by integrating the series expansion.
///
///\param p_value = pointer to the destination of the value of the sigma_inel
///\param sigma_jet = array of sigma_jet products
///
///\return 0 on successful run
///
double f_inel (double t, void * params){
    double x = t/(1-t);
    double jacobian = pow(1-t,-2);
    double sigma_jets[g_eikonal_sum_term_count+1];
    for(int i=0; i<g_eikonal_sum_term_count+1; i++){
        sigma_jets[i] = *((double *) params+i);
    }
    double variation = 4.72;
    double A = exp(-x/(4*variation))/ (4 * M_PI * variation);
    double dummy = exp(-A*sigma_jets[0]);

    double res=0;

    for(int i=1; i<=g_eikonal_sum_term_count; i++){
        int factorial = 1;
        for(int j=1; j<=i; j++){
            factorial *=j;
        }
        res += sigma_jets[i]*pow(A,i)/factorial;
    }
    return M_PI*res*dummy*jacobian;
}
///
int calculate_sigma_inel_from_integral(double * const p_value, double sigma_jet[g_eikonal_sum_term_count+1])
{
    int not_success = 0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    F.function = &f_inel;
    F.params = sigma_jet;

    gsl_integration_qags(&F, 0, 1, 0, g_error_tolerance, 1000, w, p_value, &error);

    gsl_integration_workspace_free(w);

    return not_success;
}


///calculate_sigma_tot_from_integral
///Function that does the calculation of the sigma_inel by integrating the series expansion.
///
///\param p_value = pointer to the destination of the value of the sigma_inel
///\param sigma_jet = array of sigma_jet products
///
///\return 0 on successful run
///
double f_tot (double t, void * params){
    double x = t/(1-t);
    double jacobian = pow(1-t,-2);
    double sigma_jets[g_eikonal_sum_term_count+1];
    for(int i=0; i<g_eikonal_sum_term_count+1; i++){
        sigma_jets[i] = *((double *) params+i);
    }
    double variation = 4.72;
    double A = exp(-x/(4*variation))/ (8 * M_PI * variation);
    double dummy = exp(-A*sigma_jets[0]);

    double res=0;

    for(int i=1; i<=g_eikonal_sum_term_count; i++){
        int factorial = 1;
        for(int j=1; j<=i; j++){
            factorial *=j;
        }
        res += sigma_jets[i]*pow(A,i)/factorial;
    }
    return 2*M_PI*res*dummy*jacobian;
}
///
int calculate_sigma_tot_from_integral(double * const p_value, double sigma_jet[g_eikonal_sum_term_count+1])
{
    int not_success = 0;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    F.function = &f_tot;
    F.params = sigma_jet;

    gsl_integration_qags(&F, 0, 1, 0, g_error_tolerance, 1000, w, p_value, &error);

    gsl_integration_workspace_free(w);

    return not_success;
}


///calculate_sigma_tot_from_sigma_jet
///Function that does the calculation of the sigma_tot based on sigma_jet.
///
///\param p_value = pointer to the destination of the value of the sigma_tot
///\param p_sigma_jet = pointer to the value of sigma_jet
///
///\return 0 on successful run
///
int calculate_sigma_tot_from_sigma_jet(double * const p_value, double * const p_sigma_jet)
{
    int not_success = 0;

    double variation = 4.72;
    double dummy = *p_sigma_jet / (8 * M_PI * variation);

    *p_value = (8 * M_PI * variation) * (M_EULER + log(dummy) + gsl_sf_expint_E1(dummy)) ;

    return not_success;
}


///calculate_sigma_jet
///Function that does the calculation of the sigma_jet. Calls phasespace_integral.
///
///\param p_value = pointer to the destination of the value of the sigma_jet
///\param p_error = pointer to the destination of the error output from cubature
///\param p_mand_s = pointer to the mandelstam variable s
///\param p_kt2_lower_cutoff = pointer to the lower cutoff of kt²
///\param p_pdf = pointer to the PDF object
///
///\return 0 on successful run
///
int calculate_sigma_jet(double * const p_value, double * const p_error, const double * const p_mand_s, const double * const p_kt2_lower_cutoff, PDF* p_pdf)
{
    int not_success;
    const double upper_limits [3] = {1, 1, 1};
    const double lower_limits [3] = {0, 0, 0};

    not_success = phasespace_integral(upper_limits, lower_limits, p_value, p_error, p_mand_s, p_kt2_lower_cutoff, p_pdf);
    return not_success;
}


///calculate_sigma_jet_product
///Function that does the calculation of the sigma_jets with the momentum conservation limiting.
///
///\param p_value = pointer to the destination of the value of the sigma_jet
///\param p_error = pointer to the destination of the error output from cubature
///\param p_mand_s = pointer to the mandelstam variable s
///\param p_kt2_lower_cutoff = pointer to the lower cutoff of kt²
///\param p_pdf = pointer to the PDF object
///\param nof_dijets = number of similar reactions, affects momentum conservation.
///
///\return 0 on successful run
///
int calculate_sigma_jet_product(double * p_value, double * p_error, const double * const p_mand_s, const double * const p_kt2_lower_cutoff, PDF* p_pdf, int nof_dijets)
{
    size_t dim = 3*nof_dijets;
    pair<PDF*,pair<const double * const,const double * const> > fdata = { p_pdf , { p_mand_s , p_kt2_lower_cutoff } };

    double xl[dim];
    double xu[dim];

    for (auto& x:xl) x=0;
    for (auto& x:xu) x=1;

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G = { &sigma_jet_product_integrand, dim, &fdata};

    size_t calls = 100000;

    gsl_rng_env_setup();

    T=gsl_rng_default;
    r=gsl_rng_alloc(T);

    gsl_monte_miser_state *s =gsl_monte_miser_alloc(dim);

    gsl_monte_miser_integrate(&G, xl, xu, dim, calls, r, s, p_value, p_error);

    cout<<*p_value<<endl;

    gsl_monte_miser_free(s);

    /*

    gsl_monte_vegas_state *s =gsl_monte_vegas_alloc(dim);

    gsl_monte_vegas_integrate(&G, xl, xu, dim, 10000, r, s, p_value, p_error);

    cout<<*p_value<<' '<<gsl_monte_vegas_chisq(s)<<endl;

    do{
        gsl_monte_vegas_integrate(&G, xl, xu, dim, calls, r, s, p_value, p_error);

    cout<<*p_value<<' '<<gsl_monte_vegas_chisq(s)<<endl;

    }while(fabs(gsl_monte_vegas_chisq(s) -1.0) > 0.1);

    gsl_monte_vegas_free(s);*/

//    not_success = phasespace_integral(upper_limits, lower_limits, p_value, p_error, p_mand_s, p_kt2_lower_cutoff, p_pdf);
    return 0;
}



///phasespace_integral
///Function that calculates the phasespace integral with cubature algorithm. Uses integrand_function.
///
///\param upper_limits = list of the upper limits of the integral. Order: kt2, x1, x2
///\param lower_limits = list of the lower limits of the integral. Order: kt2, x1, x2
///\param p_value = pointer to the destination of the value of the integral
///\param p_error = pointer to the destination of the error output from cubature
///\param p_mand_s = pointer to the mandelstam variable s
///\param p_kt2_lower_cutoff = pointer to the lower cutoff of kt²
///\param p_pdf = pointer to the PDF object
///
///\return 0 on successful run
///
int phasespace_integral(const double upper_limits [3], const double lower_limits [3], double * const p_value, double * const p_error, const double * const p_mand_s, const double * const p_kt2_lower_cutoff, PDF* p_pdf)
{

    const unsigned fdim=1;
    pair<PDF*,pair<const double * const,const double * const> > fdata = { p_pdf , { p_mand_s , p_kt2_lower_cutoff } };

    int not_success;


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


    if(not_success)
    {
        cout<<"Problem with integration"<<endl;
        return 1;
    }
    else
    {
        //printf("Computed integral = %0.10g +/- %g\n", *p_value, *p_error);
        return 0;
    }
}

///eikonal_integrand
///Function to be integrated by cubature algorithm.
///
///\param ndim = dimension of the integral variable
///\param p_x = pointer to the integral variable
///\param p_fdata = pointer to the additional arguments (constants in integration)
///\param fdim = dimension of the output
///\param p_fval = pointer to the output destination
///
///\return 0 on successful run
///
int eikonal_integrand(unsigned ndim, const double *p_x, void *p_fdata, unsigned fdim, double *p_fval)
{
    pair< double, pair< double , int > > fdata = *(pair< double, pair< double , int > > *) p_fdata;
    double sigma_jet = fdata.first;
    double sigma_jet_product = fdata.second.first;
    int k = fdata.second.second;

    double variation = 4.72;
    double dummy = sigma_jet / (4 * M_PI * variation);

    int factorial = 1;

    double x = *p_x/(1- *p_x);
    double jacobian = 1/pow(1- *p_x,2);

    for(int i=1; i<=k; i++)
    {
        factorial *=i;
    }
    p_fval[0] = (4 * M_PI * variation)*exp(-dummy*exp(-x)) *(pow(exp(-x)/ (4 * M_PI * variation),k)/factorial) * sigma_jet_product * jacobian;

    return 0; // success
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
    pair<PDF*,pair<const double * const,const double * const> >  fdata = *(pair<PDF*,pair<const double * const,const double * const> > *) p_fdata;
    const PDF* p_pdf = fdata.first;
    const double * const p_mand_s = fdata.second.first;
    const double * const p_kt2_lower_cutoff = fdata.second.second;
    const auto z1 = p_x[0], z2 = p_x[1], z3 = p_x[2];

    const auto kt2           = *p_kt2_lower_cutoff + z1 * ((*p_mand_s/4) - *p_kt2_lower_cutoff);
    const auto sqrt_s_per_kt = sqrt(*p_mand_s/kt2);

    const auto y1_upper = acosh(sqrt_s_per_kt/2);
    //const auto y1_lower = -y1_upper;
    const auto y1       = ( -1 + (2 * z2) ) * y1_upper; //y1_lower + z2 * (y1_upper - y1_lower)

    const auto y2_upper = log(sqrt_s_per_kt - exp(y1));
    const auto y2_lower = -log(sqrt_s_per_kt - exp(-y1));
    const auto y2       = y2_lower + z3 * (y2_upper - y2_lower);

    const auto x1 = (exp(y1) + exp(y2)) / sqrt_s_per_kt;
    const auto x2 = (exp(-y1) + exp(-y2)) / sqrt_s_per_kt;

    const auto s_hat         = s_hat_from_ys(&y1, &y2, &kt2);
    const auto t_hat         = t_hat_from_ys(&y1, &y2, &kt2);
    const auto u_hat         = u_hat_from_ys(&y1, &y2, &kt2);
    //const auto subprocess_cs = sigma_gg_gg(&s_hat, &t_hat, &u_hat, p_pdf->alphasQ2(kt2));

    const auto jacobian = ((*p_mand_s/4) - *p_kt2_lower_cutoff) * (2*y1_upper) * (y2_upper - y2_lower);

    if (y2_upper < y2_lower || y1_upper < -y1_upper) p_fval[0]=0;
    else
    {
        //SES
        /*
        p_fval[0] = 0.5 * f_ses(&x1, &kt2, p_pdf)
                    * f_ses(&x2, &kt2, p_pdf)
                    * subprocess_cs * jacobian; */
        //FULL SUMMATION
        p_fval[0] = 0.5 * diff_sigma_jet(&x1,&x2,&kt2,p_pdf,&s_hat,&t_hat,&u_hat)
                    * jacobian;
    }

    /*
        //cout <<subprocess_cs<<' ';

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

        //cout<<"kt2="<<kt2<<", x1="<<x1<<", x2="<<x2<<", s="<<s<<", t="<<t<<", u="<<u<<", cs="<<subprocess_cs<<", J="<<jacobian<<", I="<<p_fval[0]<<endl;
    */
    return 0; // success
}

///s_hat_from_ys
///Calculates and returns the subprocess mandelstam variable s. Return value can be complex.
///
///\param p_y1 = pointer to momentum fraction of the first particle
///\param p_y2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///
double s_hat_from_ys(const double * const p_y1, const double * const p_y2, const double * const p_kt2)
{
    const double y1 = *p_y1, y2 = *p_y2, kt2 = *p_kt2;
    return 2*kt2 * (1 + cosh(y1 - y2));
}


///t_hat_from_ys
///Calculates and returns the subprocess mandelstam variable t. Return value can be complex.
///
///\param p_y1 = pointer to momentum fraction of the first particle
///\param p_y2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///
double t_hat_from_ys(const double * const p_y1, const double * const p_y2, const double * const p_kt2)
{
    const double y1 = *p_y1, y2 = *p_y2, kt2 = *p_kt2;
    return -kt2 * (1 + exp(-y1 + y2));
}


///u_hat_from_ys
///Calculates and returns the subprocess mandelstam variable u. Return value can be complex.
///
///\param p_y1 = pointer to momentum fraction of the first particle
///\param p_y2 = pointer to momentum fraction of the other particle
///\param p_kt2 = pointer to the transverse momentum squared
///
double u_hat_from_ys(const double * const p_y1, const double * const p_y2, const double * const p_kt2)
{
    const double y1 = *p_y1, y2 = *p_y2, kt2 = *p_kt2;
    return -kt2 * (1 + exp(y1 - y2));
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
double diff_sigma_jet(const double * const p_x1, const double * const p_x2, const double * const p_q2, const PDF* p_pdf, const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat)
{
    const double x1 = *p_x1, x2 = *p_x2, q2 = *p_q2;
    double sum = 0;

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
            if (flavour1 != flavour2)
            {
                sum += p_pdf->xfxQ2(flavour1, x1, q2) * p_pdf->xfxQ2(flavour2, x2, q2) * sigma_qiqj_qiqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
                sum += p_pdf->xfxQ2(-flavour1, x1, q2) * p_pdf->xfxQ2(-flavour2, x2, q2) * sigma_qiqj_qiqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
                sum += p_pdf->xfxQ2(flavour1, x1, q2) * p_pdf->xfxQ2(-flavour2, x2, q2) * sigma_qiqj_qiqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
                sum += p_pdf->xfxQ2(-flavour1, x1, q2) * p_pdf->xfxQ2(flavour2, x2, q2) * sigma_qiqj_qiqj(p_s_hat,p_t_hat,p_u_hat, p_pdf->alphasQ2(q2));
            }
        }
    }

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
double sigma_qiqj_qiqj(const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat, const double alpha_s)
{
    const double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
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
double sigma_qiqi_qiqi(const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat, const double alpha_s)
{
    const double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
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
double sigma_qiaqi_qjaqj(const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat, const double alpha_s)
{
    const double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
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
double sigma_qiaqi_qiaqi(const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat, const double alpha_s)
{
    const double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
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
double sigma_qiaqi_gg(const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat, const double alpha_s)
{
    const double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
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
double sigma_gg_qaq(const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat, const double alpha_s)
{
    const double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
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
double sigma_gq_gq(const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat, const double alpha_s)
{
    const double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
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
double sigma_gg_gg(const double * const p_s_hat, const double * const p_t_hat, const double * const p_u_hat, const double alpha_s)
{
    const double s = *p_s_hat, t = *p_t_hat, u = *p_u_hat;
    return (M_PI * alpha_s * alpha_s / (s * s))*4.5*(3.0-(u*t)/(s*s)-(u*s)/(t*t)-(s*t)/(u*u));
}
