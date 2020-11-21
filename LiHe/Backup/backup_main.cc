#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
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
using std::endl;
using std::vector;

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
    //double scale = par[0];
    //double Ru = par[1];
    //double N_dc = par[2];
    //double N_li9he8 = par[3];
    //double Lifetime_li9he8 = par[4];
    //double N_b12 = par[5];
    //double Lifetime_b12 = par[6];
    //double N_n12 = par[7];
    //double Lifetime_n12 = par[8];
const char *pdf_source[4] = {
	"[2] * [1] * exp(-[1] * x )",
	"[3] * ([1] + 1 / [4]) * exp(-([1] + 1 / [4]) * x )",
	"[5] * ([1] + 1 / [6]) * exp(-([1] + 1 / [6]) * x )",
	"[7] * ([1] + 1 / [8]) * exp(-([1] + 1 / [8]) * x )",
};
double par[N_combinedFit_pars];
double TF1_pars[N_range][N_sliceTypes][N_slices][N_TF1_pars];
TF1 *array_TF1[N_range][N_sliceTypes][N_slices];
TH1D *array_h[N_range][N_sliceTypes][N_slices];
vector<double> hist_content[N_range][N_sliceTypes][N_slices];
vector<double> hist_center[N_range][N_sliceTypes][N_slices];
ROOT::Math::IMultiGenFunction *fChi2[N_range][N_sliceTypes][N_slices];
ROOT::Math::WrappedMultiTF1 *wf[N_range][N_sliceTypes][N_slices];
double array_gradient[N_range][N_sliceTypes][N_slices][N_TF1_pars];

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
void do_fit(int site, double results[N_sites][N_combinedFit_pars][2], double daily_rates[N_sites][2]);
void initialize_minimizer(int site, ROOT::Math::Minimizer *mini);
double chi2(const double *par);
void plotSlices(int site, const double fitresults[N_sites][N_combinedFit_pars][2]);

int main(int argc, char** argv){
    fill_combinedFit_parNames(combinedFit_parNames);
    showCombinedFitParNames(combinedFit_parNames);
    fill_TF1_parNames(TF1_parNames);
    showTF1ParNames(TF1_parNames);

    for(int idx = 0;idx<N_combinedFit_pars;idx++) par[idx] = idx;
    parTran(par,TF1_pars);

    double results[N_sites][N_combinedFit_pars][2];
    double daily_rates[N_sites][2];

    do_fit(1,results,daily_rates);
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
    cout<<"N_combinedFit_pars from calculation :"<< N_combinedFit_pars << "\n";
    cout<<"N_combinedFit_pars from program : "<<idx<<"\n";
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
        int sliceType = (coord - (N_sources*N_range+N_sources)) / 5;
        int slice = (coord - (N_sources*N_range+N_sources)) % 5;
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

void do_fit(int site, double results[N_sites][N_combinedFit_pars][2], double daily_rates[N_sites][2]){
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
    TFile* histogramRootFile = new TFile("lihe.data3.root");
    for(int range = 0;range<N_range;range++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++){
                array_h[range][sliceType][slice] = (TH1D*)histogramRootFile->Get(TString::Format("dtlSH_EH%i_%i_%i_%i",site,range,sliceType,slice));
                array_TF1[range][sliceType][slice] = new TF1(TString::Format("f_%i_%i_%i",range,sliceType,slice),_f_sum,N_TF1_pars);
				wf[range][sliceType][slice] = new ROOT::Math::WrappedMultiTF1(*array_TF1[range][sliceType][slice], 1);
                wf[range][sliceType][slice]->SetDerivPrecision(1e-7);

				data_chi2[range][sliceType][slice] = new ROOT::Fit::BinData(opt_chi2, rangeD);
				ROOT::Fit::FillData(*data_chi2[range][sliceType][slice], array_h[range][sliceType][slice]);
				fChi2[range][sliceType][slice] = new ROOT::Fit::Chi2FCN<ROOT::Math::IBaseFunctionMultiDim>(*data_chi2[range][sliceType][slice], *wf[range][sliceType][slice]);	
                ndf += data_chi2[range][sliceType][slice]->NPoints();
                cout<<data_chi2[range][sliceType][slice]->NPoints()<<"here ndf\n";

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
    mini->SetStrategy(2);
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
    //daily_rates[site][0] = x[2];
    //daily_rates[site][0] += x[6];
    //daily_rates[site][0] += x[10];
    double cov = 0;
    //cov = mini->CovMatrix(2,2);
    //cov += mini->CovMatrix(2,6);
    //cov += mini->CovMatrix(6,2);
    //cov += mini->CovMatrix(6,6);
    //cov += mini->CovMatrix(2,10);
    //cov += mini->CovMatrix(10,2);
    //cov += mini->CovMatrix(10,10);
    //cov += mini->CovMatrix(6,10);
    //cov += mini->CovMatrix(10,6);
    //cov = sqrt(cov);
    //daily_rates[site][1] = cov;
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
        cout << combinedFit_parNames[i] << " " << (l1-l0) / _eps << " " << likelihood_derivative(_x, i) << endl;
    }
    
    cout<<"Resutls\n";
    for(int idx = 0;idx<N_combinedFit_pars;idx++)
        cout<<"No.["<<idx<<"]:"<<combinedFit_parNames[idx]<<"="<<x[idx]<<"+-"<<x_err[idx]<<"\n";
    parTran(x,TF1_pars);

    ndf -= mini->NFree();
	cout << "Chi2/NDF: " << chi2(mini->X()) << "/" << ndf << endl;
    plotSlices(site,results);
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
        4000000.,4100000.,4100000.
    };
    double livetime[9] ={0,1544*0.82,1744*0.82,1744*0.85,1744*0.98,1744*0.98,1744*0.98,1544*0.98,1544*0.85}; 
    double livetime_site[N_sites] = {
        livetime[1]+livetime[2],
        livetime[3]+livetime[8],
        livetime[4]+livetime[5]+livetime[6]+livetime[7]
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
    double step_ratio = 1.0e-2;
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
            for(int slice = 0;slice<N_slices;slice++)
                Epsilon_source_sliceType_slice[source][sliceType][slice] = 0.2;

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
    for(int i = 0;i<N_combinedFit_pars;i++)
        printf(init_str,i,combinedFit_parNames[i],mini->X()[i], mini->Errors()[i]);
    cout<<"Number of free parameter in combined fit :"<< mini->NFree()<<"\n";
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

void plotSlices(int site, const double par_result[N_sites][N_combinedFit_pars][2]){
    double Epsilon_source_sliceType_slice[N_sources][N_sliceTypes][N_slices];
    double Epsilon_source_sliceType_slice_error[N_sources][N_sliceTypes][N_slices];
    int idx = epsBegin+1;
    for(int source = 0;source<N_sources;source++)
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++)
            for(int slice = 0;slice<N_slices;slice++){
                Epsilon_source_sliceType_slice[source][sliceType][slice] = par_result[site-1][idx][0];
                Epsilon_source_sliceType_slice_error[source][sliceType][slice] = par_result[site-1][idx][1]; idx++;
                printf("source_sliceType_slice_%i_%i_%i: %f+-%f\n",source,sliceType,slice,Epsilon_source_sliceType_slice[source][sliceType][slice],Epsilon_source_sliceType_slice_error[source][sliceType][slice]);
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
                Epsilon_source_sliceType_slice_error[source][sliceType][slice] /= tmp_epsilon_sum;
                //cout<<Epsilon_source_sliceType_slice[source][sliceType][slice]<<"\t";
            }
            //cout<<"\n";
        }

    char sourcesName[N_sources_max][100] = {
        "dc",
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
    TGraphErrors* gr_source_sliceType[N_sources][N_sliceTypes];
    for(int source = 0;source<N_sources;source++){
        for(int sliceType = 0;sliceType<N_sliceTypes;sliceType++){
            gr_source_sliceType[source][sliceType] = new TGraphErrors();
            for(int slice = 0;slice<N_slices;slice++){
                gr_source_sliceType[source][sliceType]->SetPoint(slice,slice+0.5,Epsilon_source_sliceType_slice[source][sliceType][slice]);
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
            gr_source_sliceType[source][sliceType]->SetMarkerColor(source+1);
            if(sliceType==2)
                gr_source_sliceType[source][sliceType]->SetMaximum(0.8);
            if(source==0)
                gr_source_sliceType[source][sliceType]->Draw("ap");
            else
                gr_source_sliceType[source][sliceType]->Draw("psame");
            lg->AddEntry(gr_source_sliceType[source][sliceType],sourcesName[source],"pl");
        }
        lg->Draw("same");
        c->SaveAs(TString::Format("Slice_%i_sliceType_%i.png",N_slices,sliceType));
        delete lg;
    }
}
