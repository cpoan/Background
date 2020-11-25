#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/LogLikelihoodFCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "HFitInterface.h"
using std::cout;
using std::setprecision;
using std::fixed;
using std::endl;
using std::vector;
using std::ofstream;
/*
    `   #N_range:  
            Number of muon energy ranges.
        
        #N_source: 
            Number of isotopes.
            ex: 1. Muon uncorrelated; 2. 9Li&8He
        
        #N_source_max: 
            We can add N12 and B12 to be the potential bkg contribution.
        
        #N_sliceTypes:
            Number of quantity we apply slicing. 
            If no special requirement, it would be Ep, Ed, dT, dist.
        
        #N_slices:
            The number of slices we apply to quantity.
            ex: N_slices = 4, Ep will be divided to 4 categories.
            So as the other quantities Ed, dT, dist.
        
        #N_combinedFit_pars:
            Number of free parameters in this combinedFit. Not yet consider we will fix some or not.
        
        #N_TF1_pars:
            Each single histogram in this joint fit will be given a TF1, the number of parameters that each TF1 needs is N_TF1_pars.
        
        #combinedFit_parNames[]:
            Array that contain the names of parameters in jointFit.
        
        #TF1_parNames[];
            Array that contain the names of parameters in each TF1.
        
        #scaleFactor:
            Scale the normalize factor due to effect of bin size.
        
        #fitMin, fitMax:
            Define the range of fitting. unit is "ms".
        
        #EpsBegin:
            A index that tell program where is the parameter epsilon begin.
        
        #livetime[]:
            Array that contain livetime information (Consider the muon veto efficiency already.),
            which will be used to calculate daily rates.
*/
const int N_range = 3;
const int N_sources = 2;
const int N_sources_max = 4;
const int N_sliceTypes = 4;
const int N_slices = 5;
const int N_combinedFit_pars = N_range+1+N_range*(N_sources-1)+N_sources-1+N_sources*N_sliceTypes*N_slices;
const int N_TF1_pars = 1+1+N_sources+(N_sources-1);
const int N_sites = 3;
const char* combinedFit_parNames[N_combinedFit_pars];
const char* TF1_parNames[N_TF1_pars];
const double scalefactor = 10.;
const double fitMin = 15.;
const double fitMax = 5000;
const int epsBegin = N_range+1+N_range*(N_sources-1)+N_sources-1-1; 
//double livetime[8] ={1580.34,1833.13,1907.23,2190.06,2189.62,2189.29,2005.62,1744.43};
double livetime[8] ={1277.49,1431.6,1489.37,1710.97,1710.58,1710.3,1526.58,1327.46};
//double livetime[8] ={0,1277.49,1833.13,1907.23,2190.06,2189.62,2189.29,2005.62,1744.43};
double livetime_site[N_sites] = {
        livetime[0]+livetime[1],
        livetime[2]+livetime[7],
        livetime[3]+livetime[4]+livetime[5]+livetime[6]
    };
/*
        double scale            = par[0];
        double Ru               = par[1];
        double N_dc             = par[2];
        double N_li9he8         = par[3];
        double Lifetime_li9he8  = par[4];
        -------Ususlly the parameters we needs only up to here-------
        double N_b12            = par[5];
        double Lifetime_b12     = par[6];
        double N_n12            = par[7];
        double Lifetime_n12     = par[8];
*/
const char *pdf_source[4] = {
	"[2] * [1] * exp(-[1] * x )",
	"[3] * ([1] + 1 / [4]) * exp(-([1] + 1 / [4]) * x )",
	"[5] * ([1] + 1 / [6]) * exp(-([1] + 1 / [6]) * x )",
	"[7] * ([1] + 1 / [8]) * exp(-([1] + 1 / [8]) * x )",
};

/*
        -------Global varaiable-------
        #par[]:
            Array of parameters that really been fitted,
            each site will use this global variable in order.

        #TF1_pars[N_range][N_sliceTypes][N_slices][N_TF1_pars]:
            Array of TF1 parameters,
            it's N_range x N_sliceTypes x N_slices' N_TF1_pars dimension array,
            each site will use this global variable in order.

        #array_TF1[N_range][N_sliceTypes][N_slices]:
            Array of TF1, 
            each site will use this global variable in order.

        #array_h[N_range][N_sliceTypes][N_slices]:
            Array of histograms,
            we will access all histograms from one rootfile.

        During the fitting,
        par[] is the fitting target,
        but we will parametrize them to each each TF1_pars,
        then it's more convinent to do the calculation of FCN in unit of histogram,
        and it's more convinent to do some visualization also.
        For the details of parametrization please access the Docdb (...).

        #hist_content and hist_center:
            Array of vector,
            just transfer the information of histograms to vector.

        #array_gradient[N_range][N_sliceTypes][N_slices][N_TF1_pars]:
            Save the likelihood partial derivative of parameters in TF1_pars,
            the relation of these parametetrs and true fitted parameters can
            further be used to calculate likelihood partial derivative of true fitted parameters. 
        #texfile:
            Output some important results in form of latex table.

        The rest of global parameters are less important and not essential rules in this program.
*/
double par[N_combinedFit_pars];
double TF1_pars[N_range][N_sliceTypes][N_slices][N_TF1_pars];
TF1 *array_TF1[N_range][N_sliceTypes][N_slices];
TH1D *array_h[N_range][N_sliceTypes][N_slices];
vector<double> hist_content[N_range][N_sliceTypes][N_slices];
vector<double> hist_center[N_range][N_sliceTypes][N_slices];
ROOT::Math::IMultiGenFunction *fChi2[N_range][N_sliceTypes][N_slices];
ROOT::Math::WrappedMultiTF1 *wf[N_range][N_sliceTypes][N_slices];
double array_gradient[N_range][N_sliceTypes][N_slices][N_TF1_pars];
ofstream texfile;

/*
        -------Function illustration-------
        #void fill_combinedFit_parNames(const char** array_names);
            Fill the fit parameter names to array_names.

        #void showCombinedFitParNames(const char** array_names);
            Show the fit parameter names.

        #void fill_TF1_parNames(const char** array_names);
            Fill the par names of TF1 (All same).

        #void showTF1ParNames(const char** array_names);
            Show the par names of TF1.

        #void parTran(const double* par, double tf1_pars[N_range][N_sliceTypes][N_slices][N_TF1_pars]);
            Parametrize the true fitted parameters to TF1_pars.

        #double function_single(double x, double* par);
            double* par points to the array of certain set of TF1 pars,
            using these parameters and specified the "x",
            then we will have the expected value in histograms.
        
        #void function_gradient_single(double x, double *par, double *grad);
            double* par points to the array of certain set of TF1 pars,
            using these parameters and specified the "x",
            then they will change the (double *) grad,
            "ith" index match partial derivative of "ith" parameter of single TF1 fit function.

        #double likelihood_single(vector<double> n, vector<double> x, double *par);
            'single' means it's the likelihood of one histogram,
            therefore the argument is the vector contain bin center, bin content,
            and of course the (double *) par, true fitted parameters.

        #double likelihood(const double *par);
            Real likelihood function in this adjoint fit,
            equivelant to summation of all likelihood_single.
        
        #double likelihood_derivative(const double *par, int coord);
            Funciton be given true fitted parameters and a index of parameter "coord",
            then function will return the partial derivative of likelihood with respect to this parameter[coord].
            Joint fit results obtain by using this function as the core function of searching minimum it way more stable than using likelihood(const double* par).
            
        #void do_fit(int site, double results[N_sites][N_combinedFit_pars][2], double daily_rates[N_sites][2], const char* preDir);
            Tell do fit function the site you want to perform combined fit,
            it will fill the fit result (parameters/uncertainties) to "double results[N_sites][N_combinedFit_pars][2]" the[2] in array is for fitted mean/uncertainty for parameter.
            Also the daily rates results will be contained by double daily_rates[N_sites][2].
            Somes plots will be save in ./Slices/"preDir"/ and ./FitPlot/"preDir"/.

        #void initialize_minimizer(int site, ROOT::Math::Minimizer *mini);
            Set initial value of those fitted parameters and some configuration to minimizer.

        #double chi2(const double *par);
            Calculate the chi2 to judge goodness of fit.

        #void plotSlices(int site, ROOT::Math::Minimizer* mini, const char* preDir);
        #void plotFit(int site, ROOT::Math::Minimizer* mini, const char* preDir);
            These function are called by do_fit(),
            plots will be save in ./Slices/"preDir"/ and ./FitPlot/"preDir"/.
        
*/
void fill_combinedFit_parNames(const char** array_names);
void showCombinedFitParNames(const char** array_names);
void fill_TF1_parNames(const char** array_names);
void showTF1ParNames(const char** array_names);
void parTran(const double* par, double tf1_pars[N_range][N_sliceTypes][N_slices][N_TF1_pars]);
double function_single(double x, double* par);
void function_gradient_single(double x, double *par, double *grad);
double likelihood_single(vector<double> n, vector<double> x, double *par);
double likelihood(const double *par);
double likelihood_derivative(const double *par, int coord);
void do_fit(int site, double results[N_sites][N_combinedFit_pars][2], double daily_rates[N_sites][2], const char* preDir);
void initialize_minimizer(int site, ROOT::Math::Minimizer *mini);
double chi2(const double *par);
void plotSlices(int site, ROOT::Math::Minimizer* mini, const char* preDir);
void plotFit(int site, ROOT::Math::Minimizer* mini, const char* preDir);

int main(int argc, char** argv){
    //cout<<fixed<<setprecision(2);
    //for(int site = 0;site<3;site++) cout << livetime_site[site] << " ";
    //cout<<"\n";
    if(argc!=3) std::cout<<"Wrong number of input, please retry.\n";
    texfile.open("table.tex",std::ofstream::app);
    fill_combinedFit_parNames(combinedFit_parNames);
    //showCombinedFitParNames(combinedFit_parNames);
    fill_TF1_parNames(TF1_parNames);
    //showTF1ParNames(TF1_parNames);

    for(int idx = 0;idx<N_combinedFit_pars;idx++) par[idx] = idx;
    parTran(par,TF1_pars);

    double results[N_sites][N_combinedFit_pars][2];
    double daily_rates[N_sites][2];

    char preDir[100];
    sprintf(preDir,"%s/s_%s",argv[1],argv[2]);
    printf("%s",preDir);

    do_fit(1,results,daily_rates,preDir);
    do_fit(2,results,daily_rates,preDir);
    do_fit(3,results,daily_rates,preDir);
    texfile.close();
}

void fill_combinedFit_parNames(const char** array_names){
    static char buf[N_combinedFit_pars+20][100];
    int idx = 0;
    const char rangeNames[N_range][100] = {
        "low",
        "mid",
        "high"
    };
    const char sourceNames[N_sources_max][100] = {
        "dc",
        "li9he8",
        "b12",
        "n12"
    };
    const char sliceTypeNames[N_sliceTypes][100] = {
        "prompt",
        "delay",
        "time",
        "distance"
    };
    array_names[idx]="N_dc";
    idx++;
    for(int range = 0;range<N_range;range++){
        sprintf(buf[idx],"R_u_%s",rangeNames[range]);
        array_names[idx] = buf[idx]; idx++;
        for(int source = 1;source<N_sources;source++){
            sprintf(buf[idx],"N_%s_%s",sourceNames[source],rangeNames[range]);
            array_names[idx] = buf[idx]; idx++;
        }
    }
    for(int source = 1;source<N_sources;source++){
        sprintf(buf[idx],"Lifetime_%s",sourceNames[source]);
        array_names[idx] = buf[idx]; idx++;
    }
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++){
                sprintf(buf[idx],"Epsilon_%s_%s_%i",sourceNames[source],sliceTypeNames[sliceType],slice);
                array_names[idx] = buf[idx]; idx++;
            }
    //cout<<"N_combinedFit_pars from calculation :"<< N_combinedFit_pars << "\n";
    //cout<<"N_combinedFit_pars from program : "<<idx<<"\n";
}

void showCombinedFitParNames(const char** array_names){
    cout<<"Here's combined fit parameter names.---------------------------------------\n";
    for(int idx = 0;idx<N_combinedFit_pars;idx++){
        if(idx>epsBegin-1&&(idx-epsBegin-1)%5==0) cout<<"\n";
        cout<<"No.["<<idx<<"] "<<array_names[idx]<<"\n";
    } 
    cout<<"\nEnd---------------------------------------\n";
}

void fill_TF1_parNames(const char** array_names){
    array_names[0] = "scale";
    array_names[1] = "Ru";
    array_names[2] = "N_dc";
    array_names[3] = "N_li9he8";
    array_names[4] = "Lifetime_li9he8";
    if(N_sources<3) return;
    array_names[5] = "N_b12";
    array_names[6] = "Lifetime_b12";
    array_names[7] = "N_n12";
    array_names[8] = "Lifetime_n12";
}

void showTF1ParNames(const char** array_names){
    cout<<"Here's TF1 parameter names.---------------------------------------\n";
    for(int idx = 0;idx<N_TF1_pars;idx++){
        cout<<"No.["<<idx<<"] "<<array_names[idx]<<"\n";
    } 
    cout<<"\nEnd---------------------------------------\n";
}

void parTran(const double* par, double tf1_pars[N_range][N_sliceTypes][N_slices][N_TF1_pars]){
    double N_dc;
    double R_u[N_range];
    double N_source_range[N_sources-1][N_range];
    double Lifetime_source[N_sources-1];
    double Epsilon_source_sliceType_slice[N_sources][N_sliceTypes][N_slices];
    int idx = 0;
    N_dc = par[idx];
    idx++;
    for(int range = 0;range<N_range;range++){
        R_u[range] = par[idx]; idx++;
        for(int source = 0;source<N_sources-1;source++){
            N_source_range[source][range] = par[idx]; idx++;
        }
    }
    for(int source = 0;source<N_sources-1;source++){
        Lifetime_source[source] = par[idx]; idx++; 
    }
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++){
                Epsilon_source_sliceType_slice[source][sliceType][slice] = par[idx]; idx++;
            }
    //=========================epsilon scaling==========================
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++){
            double tmp_epsilon_sum = 0;
            for(int slice = 0;slice<N_slices;slice++){
                tmp_epsilon_sum += Epsilon_source_sliceType_slice[source][sliceType][slice];
            }
            for(int slice = 0;slice<N_slices;slice++){
                Epsilon_source_sliceType_slice[source][sliceType][slice] /= tmp_epsilon_sum;
                //cout<<Epsilon_source_sliceType_slice[source][sliceType][slice]<<"\t";
            }
            //cout<<"\n";
        }
    //=========================Calculating N_dc_fake & Filling====================
    double N_dc_fake_source[N_sources-1];
    for(int range = 0;range<N_range;range++){
        for(int source = 0;source<N_sources-1;source++){
            N_dc_fake_source[source] = 0;
            for(int range2 = 0;range2<N_range;range2++){
                if(range2!=range)
                    N_dc_fake_source[source]+=N_source_range[source][range2];
            }
        }
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++){
            for(int slice = 0;slice<N_slices;slice++){
                tf1_pars[range][sliceType][slice][0] = scalefactor; 
                tf1_pars[range][sliceType][slice][1] = R_u[range];
                tf1_pars[range][sliceType][slice][2] = N_dc*Epsilon_source_sliceType_slice[0][sliceType][slice];
                for(int source = 0;source<N_sources-1;source++)
                    tf1_pars[range][sliceType][slice][2]+=N_dc_fake_source[source]*Epsilon_source_sliceType_slice[source+1][sliceType][slice];
                tf1_pars[range][sliceType][slice][3] = N_source_range[0][range]*Epsilon_source_sliceType_slice[1][sliceType][slice];
                tf1_pars[range][sliceType][slice][4] = Lifetime_source[0];
                if(N_sources<3) continue;
                tf1_pars[range][sliceType][slice][5] = N_source_range[1][range]*Epsilon_source_sliceType_slice[2][sliceType][slice];
                tf1_pars[range][sliceType][slice][6] = Lifetime_source[1];
                tf1_pars[range][sliceType][slice][7] = N_source_range[2][range]*Epsilon_source_sliceType_slice[3][sliceType][slice];
                tf1_pars[range][sliceType][slice][8] = Lifetime_source[2];
            }
        }
    }
}

double function_single(double x, double* par){
    double scale = par[0];
    double Ru = par[1];
    double N_dc = par[2];
    double N_li9he8 = par[3];
    double Lifetime_li9he8 = par[4];
    double N_b12 = 0;
    double Lifetime_b12 = 0;
    double N_n12 = 0;
    double Lifetime_n12 = 0;
    if(N_sources>3){
        N_b12 = par[5];
        Lifetime_b12 = par[6];
        N_n12 = par[7];
        Lifetime_n12 = par[8];
    }

    double term1 = N_dc * Ru * exp(-Ru * x);
    double term2 = N_li9he8 * (Ru + 1.0 / Lifetime_li9he8) * exp( -(Ru + 1.0 / Lifetime_li9he8) * x);
    double term3 = 0;
    double term4 = 0;
    if(N_sources>3){
        term3 = N_b12 * (Ru + 1.0 / Lifetime_b12) * exp( -(Ru + 1.0 / Lifetime_b12) * x);
        term4 = N_n12 * (Ru + 1.0 / Lifetime_n12) * exp( -(Ru + 1.0 / Lifetime_n12) * x);
    }

    return scale * ( term1 + term2 + term3 + term4 );
}

void function_gradient_single(double x, double *par, double *grad){

    double scale = par[0];
    double Ru = par[1];
    double N_dc = par[2];
    double N_li9he8 = par[3];
    double Lifetime_li9he8 = par[4];
    double N_b12 = 0;
    double Lifetime_b12 = 0;
    double N_n12 = 0;
    double Lifetime_n12 = 0;
    if(N_sources>3){
        N_b12 = par[5];
        Lifetime_b12 = par[6];
        N_n12 = par[7];
        Lifetime_n12 = par[8];
    }

    double term1 = N_dc * Ru * exp(-Ru * x);
    double lambda1 = Ru + 1.0 / Lifetime_li9he8;
    double term2 = N_li9he8 * lambda1 * exp( -lambda1 * x);
    double lambda2 = 0;
    double term3 = 0;
    double lambda3 = 0;
    double term4 = 0;

    if(N_sources>3){
        lambda2 = Ru + 1.0 / Lifetime_b12;
        term3 = N_b12 * lambda2 * exp( -lambda2 * x);
        lambda3 = Ru + 1.0 / Lifetime_n12;
        term4 = N_n12 * lambda3 * exp( -lambda3 * x);
    }

    grad[0] = 0;
    if(N_sources>3)
        grad[1] = scale * ( N_dc * (1 - Ru * x) * exp(-Ru *x) + 
                    N_li9he8 * (1 - lambda1 * x) * exp(-lambda1 * x) +
                    N_b12 * (1 - lambda2 * x) * exp(-lambda2 * x) +
                    N_n12 * (1 - lambda3 * x) * exp(-lambda3 * x) );
    else
        grad[1] = scale * ( N_dc * (1 - Ru * x) * exp(-Ru *x) + 
                    N_li9he8 * (1 - lambda1 * x) * exp(-lambda1 * x)); 
    grad[2] = scale * Ru * exp(-Ru * x);
    grad[3] = scale * lambda1 * exp(-lambda1 * x);
    grad[4] = 0;
    if(N_sources<3) return;
    grad[5] = scale * lambda2 * exp(-lambda2 * x);
    grad[6] = 0;
    grad[7] = scale * lambda3 * exp(-lambda3 * x);
    grad[8] = 0;
    return;
}

double likelihood_single(vector<double> n, vector<double> x, double *par){
    double r = 0;
    for(size_t b=0; b<n.size(); ++b){
        double _f = function_single(x[b], par);
        if(n[b] > 0){
            r += _f - n[b] * (1 + log(_f / n[b]));
        }else{
            r += _f;
        }
    }
    return r;
}

void likelihood_gradient_single(vector<double> n, vector<double> x, double *par, double *g){
    for(int i=0;i<N_TF1_pars;++i)
        g[i] = 0;
    double _g[N_TF1_pars];
    for(size_t b=0; b<n.size();++b){
        double _f = function_single(x[b], par);
        function_gradient_single(x[b], par, _g);
        for(int i=0;i<N_TF1_pars;++i)
            g[i] += (1 - n[b] / _f) * _g[i];
    } 
}

double likelihood(const double *par){
	parTran(par, TF1_pars);
	double result = 0;
	for(int range = 0;range<N_range;range++){
		for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++){
		    for(int slice = 0;slice<N_slices;slice++){
                double _tmp = likelihood_single(hist_content[range][sliceType][slice], hist_center[range][sliceType][slice], TF1_pars[range][sliceType][slice]);
				result += _tmp;
			}
	    }
	}
	return result;
}

double likelihood_derivative(const double *par, int coord){
	parTran(par, TF1_pars);
    double N_dc;
    double R_u[N_range];
    double N_source_range[N_sources-1][N_range];
    double Lifetime_source[N_sources-1];
    double Epsilon_source_sliceType_slice[N_sources][N_sliceTypes][N_slices];
    double Epsilon_sum_source_sliceType[N_sources][N_sliceTypes];
    int idx = 0;
    N_dc = par[idx];
    idx++;
    for(int range = 0;range<N_range;range++){
        R_u[range] = par[idx]; idx++;
        for(int source = 0;source<N_sources-1;source++){
            N_source_range[source][range] = par[idx]; idx++;
        }
    }
    for(int source = 0;source<N_sources-1;source++){
        Lifetime_source[source] = par[idx]; idx++; 
    }
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++){
                Epsilon_source_sliceType_slice[source][sliceType][slice] = par[idx]; idx++;
            }
    //=========================epsilon scaling==========================
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++){
            double tmp_epsilon_sum = 0;
            Epsilon_sum_source_sliceType[source][sliceType] = 0;
            for(int slice = 0;slice<N_slices;slice++){
                tmp_epsilon_sum += Epsilon_source_sliceType_slice[source][sliceType][slice];
            }
            Epsilon_sum_source_sliceType[source][sliceType] = tmp_epsilon_sum;
            for(int slice = 0;slice<N_slices;slice++){
                Epsilon_source_sliceType_slice[source][sliceType][slice] /= tmp_epsilon_sum;
            }
        }
    //=====================Epsilon preparing over=======================
    for(int range = 0;range<N_range;range++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice  = 0;slice<N_slices;slice++) 
                likelihood_gradient_single(
                        hist_content[range][sliceType][slice],
                        hist_center[range][sliceType][slice],
                        TF1_pars[range][sliceType][slice],
                        array_gradient[range][sliceType][slice]
                        );
    double _sum = 0;
    //===========Here calculation for gradient from Ndc==================
    if(coord == 0){
        for(int range = 0;range<N_range;range++)
            for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
                for(int slice = 0;slice<N_slices;slice++)
                    _sum += Epsilon_source_sliceType_slice[0][sliceType][slice]*array_gradient[range][sliceType][slice][2];
        return _sum;
    }else
    //===========Here calculation for gradient from Ru==================
    if(coord == 1 || coord == 1+N_sources || coord == 1+2*N_sources ){
        int range;
        range = (coord-1)/N_sources;
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++)
                _sum += array_gradient[range][sliceType][slice][1];
        return _sum;
    }else if( coord<=N_sources*N_range){ 
        int range = (coord-1) / N_sources;
        int source = (coord-1) % N_sources;
        for(int range2 = 0;range2<N_range;range2++)
            for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
                for(int slice = 0;slice<N_slices;slice++){
                    if(range2==range)
                        _sum+=Epsilon_source_sliceType_slice[source][sliceType][slice]*array_gradient[range2][sliceType][slice][3+2*(source-1)];
                    else
                        _sum+=Epsilon_source_sliceType_slice[source][sliceType][slice]*array_gradient[range2][sliceType][slice][2];
                }
        return _sum;
    }else if(coord>N_sources*N_range&&coord<N_sources*N_range+N_sources){
        return 0;
    }else if(coord>=N_sources*N_range+N_sources&&coord<(N_sources*N_range+N_sources+N_slices*N_sliceTypes)){
        int sliceType = (coord - (N_sources*N_range+N_sources)) / N_slices;
        int slice = (coord - (N_sources*N_range+N_sources)) % N_slices;
        //printf("dc_coord_%i_sliceType_slice_%i_%i\n",coord,sliceType,slice);
        for(int range = 0;range<N_range;range++)
            for(int slice2 = 0;slice2<N_slices;slice2++){
                if(slice==slice2)
                    _sum+=array_gradient[range][sliceType][slice2][2]*(Epsilon_sum_source_sliceType[0][sliceType]-Epsilon_source_sliceType_slice[0][sliceType][slice]*Epsilon_sum_source_sliceType[0][sliceType]*(1.))/pow(Epsilon_sum_source_sliceType[0][sliceType],2)*N_dc;
                else
                    _sum+=array_gradient[range][sliceType][slice2][2]*(-1.*Epsilon_source_sliceType_slice[0][sliceType][slice2]*Epsilon_sum_source_sliceType[0][sliceType]*(1.))/pow(Epsilon_sum_source_sliceType[0][sliceType],2)*N_dc;
            }
        return _sum;
    }else if(coord>=N_sources*N_range+N_sources+N_slices*N_sliceTypes){
        int begin = N_sources*N_range+N_sources+N_slices*N_sliceTypes;
        int source = (coord - begin) / (N_sliceTypes*N_slices);
        int sliceType = (coord - begin) % (N_sliceTypes*N_slices) / N_slices;
        int slice = (((coord - begin) % (N_sliceTypes*N_slices) )% N_slices);
        //printf("coord_%i_source_sliceType_slice_%i_%i_%i\n",coord,source,sliceType,slice);
        for(int range = 0;range<N_range;range++){
            for(int slice2 = 0;slice2<N_slices;slice2++){
                for(int range2=0;range2<N_range;range2++){
                    if(range2!=range){
                        if(slice2!=slice)
                            _sum+=array_gradient[range][sliceType][slice2][2]*N_source_range[source][range2]*(-1.*Epsilon_sum_source_sliceType[source+1][sliceType]*Epsilon_source_sliceType_slice[source+1][sliceType][slice2]*(1.))/pow(Epsilon_sum_source_sliceType[source+1][sliceType],2);
                        else
                            _sum+=array_gradient[range][sliceType][slice2][2]*N_source_range[source][range2]*(Epsilon_sum_source_sliceType[source+1][sliceType]-Epsilon_sum_source_sliceType[source+1][sliceType]*Epsilon_source_sliceType_slice[source+1][sliceType][slice]*(1.))/pow(Epsilon_sum_source_sliceType[source+1][sliceType],2);
                    }else{
                        if(slice2!=slice)
                            _sum+=array_gradient[range][sliceType][slice2][3+source*2]*N_source_range[source][range2]*(-1.*Epsilon_sum_source_sliceType[source+1][sliceType]*Epsilon_source_sliceType_slice[source+1][sliceType][slice2]*(1.))/pow(Epsilon_sum_source_sliceType[source+1][sliceType],2);
                        else
                            _sum+=array_gradient[range][sliceType][slice2][3+source*2]*N_source_range[source][range2]*(Epsilon_sum_source_sliceType[source+1][sliceType]-Epsilon_sum_source_sliceType[source+1][sliceType]*Epsilon_source_sliceType_slice[source+1][sliceType][slice]*(1.))/pow(Epsilon_sum_source_sliceType[source+1][sliceType],2);
                    }
                }
            }
        }
        return _sum;
    }
    return _sum;
}

void do_fit(int site, double results[N_sites][N_combinedFit_pars][2], double daily_rates[N_sites][2], const char* preDir){
    int ndf = 0;

	ROOT::Fit::DataRange rangeD(fitMin, fitMax);
	ROOT::Fit::BinData *data_chi2[N_range][N_sliceTypes][N_slices];
	ROOT::Fit::DataOptions opt_chi2;

	char _f_sum[512] = "";
    char _f_tmp[512] = "";
	for(int source=0;source<N_sources;source++){
        if(source!=0) 
            strcat(_f_tmp, " + ");
        strcat(_f_tmp, pdf_source[source]);
    }
    sprintf(_f_sum,"[0] * ( %s )", _f_tmp);
    TFile* histogramRootFile = new TFile("./p17b.root","READ");
    for(int range = 0;range<N_range;range++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++){
                array_h[range][sliceType][slice] = (TH1D*)histogramRootFile->Get(TString::Format("dtlSH_EH%i_Nslice%i_%i_%i_%i",site,N_slices,range,sliceType,slice));
                array_TF1[range][sliceType][slice] = new TF1(TString::Format("f_%i_%i_%i",range,sliceType,slice),_f_sum,fitMin,fitMax);
				wf[range][sliceType][slice] = new ROOT::Math::WrappedMultiTF1(*array_TF1[range][sliceType][slice], 1);
                wf[range][sliceType][slice]->SetDerivPrecision(1e-7);

				data_chi2[range][sliceType][slice] = new ROOT::Fit::BinData(opt_chi2, rangeD);
				ROOT::Fit::FillData(*data_chi2[range][sliceType][slice], array_h[range][sliceType][slice]);
				fChi2[range][sliceType][slice] = new ROOT::Fit::Chi2FCN<ROOT::Math::IBaseFunctionMultiDim>(*data_chi2[range][sliceType][slice], *wf[range][sliceType][slice]);	
                ndf += data_chi2[range][sliceType][slice]->NPoints();
                //cout<<data_chi2[range][sliceType][slice]->NPoints()<<"here ndf\n";

                size_t N_bins = array_h[range][sliceType][slice]->GetNbinsX();
                hist_content[range][sliceType][slice].clear();
                hist_content[range][sliceType][slice].reserve(N_bins);
                hist_center[range][sliceType][slice].clear();
                hist_center[range][sliceType][slice].reserve(N_bins);
                for(size_t bin=1;bin<=N_bins;bin++){
                   double _center = array_h[range][sliceType][slice]->GetBinCenter(bin);
                    if(_center < fitMin || _center > fitMax) continue;
                    hist_content[range][sliceType][slice].push_back(array_h[range][sliceType][slice]->GetBinContent(bin));
                    hist_center[range][sliceType][slice].push_back(array_h[range][sliceType][slice]->GetBinCenter(bin));
                }
            }

	//ROOT::Math::Functor f_tmp(&likelihood, N_combinedFit_pars);
	//ROOT::Math::Minimizer *mini = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
    ROOT::Math::GradFunctor f_tmp(&likelihood, &likelihood_derivative, N_combinedFit_pars);
    ROOT::Math::Minimizer* mini = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    mini->SetFunction(f_tmp);
    initialize_minimizer(site, mini);

    mini->Minimize();
    mini->SetStrategy(3);
    mini->Hesse();
    mini->Minimize();
    mini->Hesse();
    mini->Minimize();
    mini->Hesse();
    cout<<"COV STATUS"<<mini->CovMatrixStatus()<<endl;
    const double *x = mini->X();
    const double *x_err = mini->Errors();
    for(int idx = 0;idx<N_combinedFit_pars;idx++){
        results[site-1][idx][0] = x[idx];
        results[site-1][idx][1] = x_err[idx];
    }
    double cov = 0;
    daily_rates[site-1][0] = 0;
    for(int range = 0; range<N_range;range++){
        if(range!=2){
            daily_rates[site-1][0] += x[2+N_sources*range];
            cov+=pow(x_err[2+N_sources*range],2);
        }

        else{
            daily_rates[site-1][0] += 0.21*x[2+N_sources*range];
            cov+=pow(0.21*x_err[2+N_sources*range],2);
        }
    }
    daily_rates[site-1][0] = (daily_rates[site-1][0])/(livetime_site[site-1]);
    daily_rates[site-1][1] = sqrt(cov)/(livetime_site[site-1]);

    //====================Here is derivative checking========
    double _x[1000];
    for(int i = 0;i<N_combinedFit_pars;i++)
        _x[i] = mini->X()[i];
    for(int i=0;i<N_combinedFit_pars;++i){
        double _eps = _x[i] * 1e-6;
        double l0 = likelihood(_x);
        _x[i] += _eps;
        double l1 = likelihood(_x);
        _x[i] -= 2*_eps;
        //cout << combinedFit_parNames[i] << " " << (l1-l0) / _eps << " " << likelihood_derivative(_x, i) << endl;
    }
    
    //cout<<"Resutls\n";
    //for(int idx = 0;idx<N_combinedFit_pars;idx++)
    //    cout<<"No.["<<idx<<"]:"<<combinedFit_parNames[idx]<<"="<<x[idx]<<"+-"<<x_err[idx]<<"\n";
    parTran(x,TF1_pars);

    ndf -= mini->NFree();
	cout << "Chi2/NDF: " << chi2(mini->X()) << "/" << ndf << endl;
    plotSlices(site,mini,preDir);
    plotFit(site,mini,preDir);
    if(site==1)
        texfile<<"\\hline\n";
    if(site==2)
        texfile<<"\t"<<N_slices<<"\t";
    else
        texfile<<"\t\t";
    texfile<<"&\tEH"<<site<<"\t&\t"<<fixed<<setprecision(2)<<daily_rates[site-1][0]<<"+-"<<setprecision(3)<<daily_rates[site-1][1]<<"\t&\t"<<setprecision(0)<<chi2(mini->X())<<"/"<<ndf<<"\\\\\n";
}

void initialize_minimizer(int site, ROOT::Math::Minimizer *mini){

	const char *init_str = "Initializing parameter %2d %s %10.2e %10.2e\n";

    cout<<"INITIALIZING MINIMIZER\n"; 
    mini->SetErrorDef(0.5);
    mini->SetTolerance(0.01);
    mini->SetPrintLevel(2);
    mini->SetStrategy(1);
    mini->SetMaxFunctionCalls(1000000);
    mini->SetMaxIterations(1000000);
    //=================Create initial value===================
    const double R_u_init_site_range[N_sites][N_range] = {
        { 1.178e-2, 8.613e-3, 1.441e-4 },
        { 8.118e-3, 7.128e-3, 1.412e-4 },
        { 5.825e-4, 4.615e-4, 1.373e-5 }
    };
    const double N_dc_init_site[N_sites] = {
        3500000.,3500000.,1450000.
    };
    const double N_init_site_source_range[N_sites][N_sources_max-1][N_range] = {
        {  //site range source
            { 0.01*livetime_site[0], 1.65*livetime_site[0], 4.60*livetime_site[0]},
            { 1.58*livetime_site[0], 1.66*livetime_site[0], 1.36*livetime_site[0]},
            { 0.84*livetime_site[0], 0.11*livetime_site[0], 0.92*livetime_site[0]}
        },
        {
            { 0.01*livetime_site[1], 1.14*livetime_site[1], 3.71*livetime_site[1]},
            { 1.14*livetime_site[1], 1.30*livetime_site[1], 1.29*livetime_site[1]},
            { 0.05*livetime_site[1], 0.05*livetime_site[1], 0.47*livetime_site[1]}
        },
        {
            { 0.01*livetime_site[2], 0.15*livetime_site[2], 0.45*livetime_site[2]},
            { 0.09*livetime_site[2], 0.13*livetime_site[2], 0.17*livetime_site[2]},
            { 0.01*livetime_site[2], 0.01*livetime_site[2], 0.09*livetime_site[2]}
        }
    };
    //=================Start Filling==========================
    double step_ratio = 1.0e-4;
    double N_dc;
    double R_u[N_range];
    double N_source_range[N_sources-1][N_range];
    double Lifetime_source[N_sources_max-1] = {
        256.366, 29.142, 15.9
    };
    double Epsilon_source_sliceType_slice[N_sources][N_sliceTypes][N_slices];
    N_dc = N_dc_init_site[site-1];
    for(int range = 0;range<N_range;range++){
        R_u[range] = R_u_init_site_range[site-1][range];
        for(int source = 0;source<N_sources-1;source++)
            N_source_range[source][range] = N_init_site_source_range[site-1][source][range];
    }
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++){
                //if(sliceType==1&&source==1){
                //    if(slice==0)
                //        Epsilon_source_sliceType_slice[source][sliceType][slice] = 10.0;
                //    else if(slice==1)
                //        Epsilon_source_sliceType_slice[source][sliceType][slice] = 1.0;
                //    else if(slice==2)
                //        Epsilon_source_sliceType_slice[source][sliceType][slice] = 4.0;
                //    else
                //        Epsilon_source_sliceType_slice[source][sliceType][slice] = 0.1;

                //}else
                    Epsilon_source_sliceType_slice[source][sliceType][slice] = 5.0;
            }

    int idx = 0;
    //N_dc = par[idx];
    mini->SetLowerLimitedVariable(idx,combinedFit_parNames[idx],N_dc,N_dc*step_ratio,0);
    idx++;
    for(int range = 0;range<N_range;range++){
        //R_u[range] = par[idx]; idx++;
        mini->SetLowerLimitedVariable(idx,combinedFit_parNames[idx],R_u[range],R_u[range]*step_ratio, 0);
        idx++;
        for(int source = 0;source<N_sources-1;source++){
            //N_source_range[source][range] = par[idx]; idx++;
            mini->SetLowerLimitedVariable(idx,combinedFit_parNames[idx],N_source_range[source][range],N_source_range[source][range]*step_ratio,0);
            idx++;
        }
    }
    for(int source = 0;source<N_sources-1;source++){
        //Lifetime_source[source] = par[idx]; idx++; 
        mini->SetFixedVariable(idx,combinedFit_parNames[idx],Lifetime_source[source]);
        idx++;
    }
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++){
                //Epsilon_source_sliceType_slice[source][sliceType][slice] = par[idx]; idx++;
                mini->SetLowerLimitedVariable(
                        idx,
                        combinedFit_parNames[idx],
                        Epsilon_source_sliceType_slice[source][sliceType][slice],
                        Epsilon_source_sliceType_slice[source][sliceType][slice]*step_ratio,
                        0);
                idx++;
            }
    //for(int i = 0;i<N_combinedFit_pars;i++)
    //    printf(init_str,i,combinedFit_parNames[i],mini->X()[i], mini->Errors()[i]);
    //cout<<"Number of free parameter in combined fit :"<< mini->NFree()<<"\n";
}

double chi2(const double *par){

	parTran(par, TF1_pars);
	double result = 0;
    for(int range = 0; range<N_range;range++)
        for(int sliceType = 0; sliceType<N_sliceTypes; sliceType++)
            for(int slice = 0; slice<N_slices;slice++)
                result += (*fChi2[range][sliceType][slice])(TF1_pars[range][sliceType][slice]);

	return result;

}

void plotSlices(int site, ROOT::Math::Minimizer* mini, const char* preDir){
    char location[100];
    sprintf(location,"./Slices/%s/",preDir);

    const double* par = mini->X();
    const double* par_err = mini->Errors();
    double Epsilon_source_sliceType_slice[N_sources][N_sliceTypes][N_slices];
    double Epsilon_source_sliceType_slice_error[N_sources][N_sliceTypes][N_slices];
    double Epsilon_sum_source_sliceType[N_sources][N_sliceTypes];
    int idx = epsBegin+1;
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++){
                Epsilon_source_sliceType_slice[source][sliceType][slice] = par[idx];
                Epsilon_source_sliceType_slice_error[source][sliceType][slice] = par_err[idx]; idx++;
                //printf("source_sliceType_slice_%i_%i_%i: %f+-%f\n",source,sliceType,slice,Epsilon_source_sliceType_slice[source][sliceType][slice],Epsilon_source_sliceType_slice_error[source][sliceType][slice]);
            }
    //=========================epsilon scaling==========================
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++){
            double tmp_epsilon_sum = 0;
            for(int slice = 0;slice<N_slices;slice++){
                tmp_epsilon_sum += Epsilon_source_sliceType_slice[source][sliceType][slice];
            }
            Epsilon_sum_source_sliceType[source][sliceType] = tmp_epsilon_sum;
            for(int slice = 0;slice<N_slices;slice++){
                double term1 = pow(Epsilon_source_sliceType_slice_error[source][sliceType][slice]/Epsilon_source_sliceType_slice[source][sliceType][slice],2);
                int index = (epsBegin+1)+source*N_sliceTypes*N_slices+sliceType*N_slices+slice;
                int start_idx = (epsBegin+1)+source*N_sliceTypes*N_slices+sliceType*N_slices;
                int end_idx = start_idx+N_slices; 
                //printf("%i_%i_%i\n",index,start_idx,end_idx);
                double term2 = 0;
                for(int i = start_idx;i<end_idx;i++)
                    for(int j = start_idx;j<end_idx;j++)
                        term2 += mini->CovMatrix(i,j);
                term2 /= pow(Epsilon_sum_source_sliceType[source][sliceType],2);
                double term3 = 0;
                for(int i = start_idx;i<end_idx;i++)
                    term3 += mini->CovMatrix(index,i);
                term3 /= (Epsilon_sum_source_sliceType[source][sliceType]*Epsilon_source_sliceType_slice[source][sliceType][slice]);
                term3 = 2.*term3;
                //term2 = 0;
                //term3 = 0;
                Epsilon_source_sliceType_slice[source][sliceType][slice] /= Epsilon_sum_source_sliceType[source][sliceType];
                Epsilon_source_sliceType_slice_error[source][sliceType][slice] = Epsilon_source_sliceType_slice[source][sliceType][slice]*sqrt(term1+term2-term3);
            }
        }

    char sourcesName[N_sources_max][100] = {
        "Uncorrelated",
        "Li9He8",
        "b12",
        "n12"
    };

    char plot_Title[N_sliceTypes][255] = {
        "Prompt energy distribution",
        "Delay energy distribution",
        "Capture time distribution",
        "p-d distance distribution"
    };
    char xaxis_title[N_sliceTypes][255] = {
        "MeV",
        "MeV",
        "{#mu}s",
        "mm"
    };
    double start[N_sliceTypes] = {
        1.5,
        1.9,
        1000,
        0
    };
    double ratio[N_sliceTypes] = {
        (12-0.5)/N_slices,
        (12-1.9)/N_slices,
        (400000-1000)/N_slices,
        (500)/N_slices
    };
    TGraphErrors* gr_source_sliceType[N_sources][N_sliceTypes];
    for(int source = 0;source<N_sources;source++){
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++){
            gr_source_sliceType[source][sliceType] = new TGraphErrors();
            for(int slice = 0;slice<N_slices;slice++){
                gr_source_sliceType[source][sliceType]->SetPoint(slice,(slice+0.5)*ratio[sliceType]+start[sliceType],Epsilon_source_sliceType_slice[source][sliceType][slice]);
                gr_source_sliceType[source][sliceType]->SetPointError(slice,0,Epsilon_source_sliceType_slice_error[source][sliceType][slice]);
            }
        }
    }
    TLegend* lg;
    TCanvas* c = new TCanvas("c","c",1200,800);
    c->cd();
    gPad->SetGrid(1);
    for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++){
            lg = new TLegend(0.8,0.8,1.0,1.0);
        for(int source = 0;source<N_sources;source++){
            gr_source_sliceType[source][sliceType]->SetTitle(plot_Title[sliceType]);
            gr_source_sliceType[source][sliceType]->SetMarkerStyle(8);
            gr_source_sliceType[source][sliceType]->SetMarkerColor(source+1);
            if(source==0)
                gr_source_sliceType[source][sliceType]->SetLineColor(kBlue);
            else
                gr_source_sliceType[source][sliceType]->SetLineColor(kGreen);
            gr_source_sliceType[source][sliceType]->SetLineWidth(3);
                gr_source_sliceType[source][sliceType]->SetMinimum(0.0);
            if(sliceType==2)
                gr_source_sliceType[source][sliceType]->SetMaximum(0.8);
            else if(sliceType==3)
                gr_source_sliceType[source][sliceType]->SetMaximum(0.55);
            if(source==0){
                gr_source_sliceType[source][sliceType]->Draw("apc");
                gr_source_sliceType[source][sliceType]->GetXaxis()->SetTitle(xaxis_title[sliceType]);
            }
            else
                gr_source_sliceType[source][sliceType]->Draw("pcsame");
            lg->AddEntry(gr_source_sliceType[source][sliceType],sourcesName[source],"pl");
        }
        lg->Draw("same");
        c->SaveAs(TString::Format("%s_EH%i_sliceType_%i.png",location,site,sliceType));
        delete lg;
    }
    TH2D *h_cov = new TH2D("h_cov","h_cov",N_combinedFit_pars,0,N_combinedFit_pars,N_combinedFit_pars,0,N_combinedFit_pars);
    for(int i = 0;i<N_combinedFit_pars;i++)
        for(int j = 0;j<N_combinedFit_pars;j++)
            h_cov->SetBinContent(i,j,mini->CovMatrix(i,j)/(mini->Errors()[i]*mini->Errors()[j]));
    gStyle->SetOptStat(0);
    h_cov->Draw("colz");
    c->SaveAs(TString::Format("CovMatrix_EH%i.png",site));
    delete c;
}

void plotFit(int site, ROOT::Math::Minimizer* mini, const char* preDir){
    char location[100];
    sprintf(location,"./FitPlot/%s/",preDir);
	char _f_sum[512] = "";
    char _f_tmp[512] = "";
	for(int source=0;source<N_sources;source++){
        if(source!=0) 
            strcat(_f_tmp, " + ");
        strcat(_f_tmp, pdf_source[source]);
    }
    sprintf(_f_sum,"[0] * ( %s )", _f_tmp);

    char fdc_name[512]="";
    sprintf(_f_tmp,"%s",pdf_source[0]);
    sprintf(fdc_name,"[0] * ( %s )", _f_tmp);
    char flihe_name[512]="";
    sprintf(_f_tmp,"%s",pdf_source[1]);
    sprintf(flihe_name,"[0] * ( %s )", _f_tmp);

    parTran(mini->X(),TF1_pars);
    TF1* fdc[N_range][N_sliceTypes][N_slices];
    TF1* flihe[N_range][N_sliceTypes][N_slices];

    TCanvas* c1 = new TCanvas("c1","c1",1200,800);
    TCanvas* c2 = new TCanvas("c2","c2",N_slices*600,N_sliceTypes*450);
    gStyle->SetOptStat(1);
    for(int range = 0;range<N_range;range++){
        c2->Clear();
        c2->Divide(N_slices,N_sliceTypes);
        int c2count = 1;
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++){
            for(int slice = 0;slice<N_slices;slice++){
                fdc[range][sliceType][slice] = new TF1(TString::Format("%i_%i_%i",range,sliceType,slice)/*,TString::Format("dc_%i_%i_%i",range,sliceType,slice)*/,fdc_name,fitMin,fitMax);
                fdc[range][sliceType][slice]->SetLineColor(kBlue);
                flihe[range][sliceType][slice] = new TF1(TString::Format("%i_%i_%i",range,sliceType,slice)/*,TString::Format("lihe_%i_%i_%i",range,sliceType,slice)*/,flihe_name,fitMin,fitMax);
                flihe[range][sliceType][slice]->SetLineColor(kGreen); 
                array_TF1[range][sliceType][slice]->SetLineColor(kRed);
                for(int idx = 0;idx<N_TF1_pars;idx++){
                    array_TF1[range][sliceType][slice]->SetParameter(idx,TF1_pars[range][sliceType][slice][idx]);
                    fdc[range][sliceType][slice]->SetParameter(idx,TF1_pars[range][sliceType][slice][idx]);
                    flihe[range][sliceType][slice]->SetParameter(idx,TF1_pars[range][sliceType][slice][idx]);
                }
                c1->cd();
                array_h[range][sliceType][slice]->SetAxisRange(10,fitMax);
                array_h[range][sliceType][slice]->SetMaximum(array_h[range][sliceType][slice]->GetBinContent(1)*1./0.7);
                array_TF1[range][sliceType][slice]->SetMinimum(0.);
                fdc[range][sliceType][slice]->SetMinimum(0.);
                flihe[range][sliceType][slice]->SetMinimum(0.);
                array_h[range][sliceType][slice]->GetXaxis()->SetMoreLogLabels(1);
                array_h[range][sliceType][slice]->GetXaxis()->SetTitle("ms");
                array_h[range][sliceType][slice]->GetYaxis()->SetTitle("Counts");
                array_h[range][sliceType][slice]->GetXaxis()->CenterTitle(1);
                array_h[range][sliceType][slice]->GetYaxis()->CenterTitle(1);
                array_h[range][sliceType][slice]->GetXaxis()->SetTitleOffset(1.30);
                array_h[range][sliceType][slice]->GetYaxis()->SetTitleOffset(1.45);
                array_h[range][sliceType][slice]->SetMinimum(0.);
                array_h[range][sliceType][slice]->Draw("E");
                array_TF1[range][sliceType][slice]->Draw("same");
                fdc[range][sliceType][slice]->Draw("same");
                flihe[range][sliceType][slice]->Draw("same");
                gPad->SetLogx();
                c1->SaveAs(TString::Format("%shist_EH%i_%i_%i_%i.png",location,site,range,sliceType,slice));
                c1->Clear();
                c2->cd(c2count++);
                gPad->SetLogx();
                array_h[range][sliceType][slice]->Draw("E");
                array_TF1[range][sliceType][slice]->Draw("same");
                fdc[range][sliceType][slice]->Draw("same");
                flihe[range][sliceType][slice]->Draw("same");
            }
        }
        c2->SaveAs(TString::Format("./thesisPlot/s%i/EH%i_range%i.png",N_slices,site,range));
    }
}
