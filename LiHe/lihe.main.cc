// this is made in 2018/05/31 
// I (PoAn,chenpoan@outlook.com), cheange the time cut window to 1000us,
// so as to fit the time distribution,
// then extract the time cut efficiency at time cut set at 400us.
//
// you can use bash command "make timemain" to generate timemain, which is operating file.
#include <iostream>
#include <fstream>
#include "TChain.h"
#include "TFile.h"
#include "AdSimple.C"
#include "CalibStats.C"
#include "PhyEvent.h"
#include <cmath>
#include "TH1D.h"
#include "TString.h"

bool EH1=0;
bool EH2=0;
bool EH3=1;
int EH=0;

using namespace std;
const long long preVetoOfWP = 2000;
const long long pVetoOfWP   = 600000;
const long long pVetoOfAD   = 1400000;
const long long pVetoOfSH   = 400000000;
const long long timeCut     = 400000;
const long long LongTimeCut = 400000;
const double    distCut     = 500;

TH1D* h_high = new TH1D("h_high","h_high",500,0,5000); 
TH1D* h_mid = new TH1D("h_mid","h_mid",500,0,5000); 
TH1D* h_low = new TH1D("h_low","h_low",500,0,5000); 

void Fill_dtlSH(TH1D* h[5][3][3][4][8], double dt_high, double dt_mid, double dt_low, double prompt, double delay, double dt, double distance);

bool epLock[5]= {false,false,false,false,false};
bool edLock[5]= {false,false,false,false,false};
bool IBDReady[5]= {false,false,false,false,false};
bool SingleReady[5]= {false,false,false,false,false};
bool FastNepLock[5]= {false,false,false,false,false};

int main(int argc, char** argv){
    if(EH1)EH=0;
    if(EH3)EH=3;

    /////////////////parameter space//////////////////
    double a_FastN_ep[5]={0,0,0,0,0};    
    double a_FastN_xp[5]={0,0,0,0,0};    
    double a_FastN_yp[5]={0,0,0,0,0};    
    double a_FastN_zp[5]={0,0,0,0,0};    
    double a_FastN_tp[5]={0,0,0,0,0};    

    double a_ep[5],a_xp[5],a_yp[5],a_zp[5],a_ed[5],a_xd[5],a_yd[5],a_zd[5];
    long long a_tlSH_high[5],a_tlSH_mid[5],a_tlSH_low[5];
    long long a_tlSH_high_tmp[5],a_tlSH_mid_tmp[5],a_tlSH_low_tmp[5];
    double a_dist[5];
    long long a_dT[5],a_tp[5],a_td[5];
    int detNo;
    double a_ibd_ep[5],a_ibd_xp[5],a_ibd_yp[5],a_ibd_zp[5],a_ibd_ed[5],a_ibd_xd[5],a_ibd_yd[5],a_ibd_zd[5];
    long long a_ibd_dtlSH_high[5],a_ibd_dtlSH_mid[5],a_ibd_dtlSH_low[5];
    long long a_ibd_dT[5],a_ibd_tp[5],a_ibd_td[5];

    double ep,xp,yp,zp,ed,xd,yd,zd;
    double dtlSH_high,dtlSH_mid,dtlSH_low;
    double dist;
    long long dT,tp,td;

    double a_E[5],a_X[5],a_Y[5],a_Z[5];
    double E,X,Y,Z;
    bool isVeto;
    long long a_T[5];
    long long T;
    int sdetNo;

    int theDetNo;
    long long theAbsRunTime,theLiveTime,totalTime;

    long long Date,AmcCount,liveTimeCount,N_up,N_down,N_07,N_15;
    long long NIBD              = 0;

    long long NmuonArray[5]={0,0,0,0,0};
    long long Nmuon = 0;

    long long NaccSingleArr[5]={0,0,0,0,0};
    long long NaccSingle=0;
    long long a_NIBD[5]={0,0,0,0,0};

    double multiEff = 0;
    double multiEff_Err = 0;
    long long Nsample[5]={0,0,0,0,0};
    long long Nsuccess[5]={0,0,0,0,0};

    ////////////////parameter space end////////////////

    ifstream list(argv[1]);
    string fname;

    TChain *adch = new TChain("/Event/Rec/AdSimple");
    TChain *csch = new TChain("/Event/Data/CalibStats");

    while(getline(list,fname)){
        cout<<fname<<endl;
        adch->Add(fname.c_str());
        csch->Add(fname.c_str());
    }

    AdSimple adrec(adch);
    CalibStats calst(csch);

    TFile* output = new TFile(argv[2],"RECREATE");
    output->cd();
    /////Li9He8 histogram//////
    for(int nslice = 4;nslice<9;nslice++)
        for(int site = 0;site<3;site++)
            for(int range = 0;range<3;range++)
                for(int sliceType = 0;sliceType<4;sliceType++)
                    for(int slice = 0;slice<nslice;slice++)
                            array_h[nslice-4][site][range][sliceType][slice] = new TH1D(TString::Format("dtlSH_EH%i_Nslice%i_%i_%i_%i",site+1,nslice,range,sliceType,slice),TString::Format("dtlSH_EH%i_%i_%i_%i",site+1,range,sliceType,slice),500,0,5000);
    /////////////////Fast Neutron Relevant/////////////
    TH1D* prompttest = new TH1D("prompttest","prompttest",50,0,100);
    TH1D* a_FastN_Prompt[8] = {0};
    a_FastN_Prompt[0] = new TH1D("FastN_AD1","FastN_AD1",50,0,100);
    a_FastN_Prompt[1] = new TH1D("FastN_AD2","FastN_AD2",50,0,100);
    a_FastN_Prompt[2] = new TH1D("FastN_AD3","FastN_AD3",50,0,100);
    a_FastN_Prompt[3] = new TH1D("FastN_AD4","FastN_AD4",50,0,100);
    a_FastN_Prompt[4] = new TH1D("FastN_AD5","FastN_AD5",50,0,100);
    a_FastN_Prompt[5] = new TH1D("FastN_AD6","FastN_AD6",50,0,100);
    a_FastN_Prompt[6] = new TH1D("FastN_AD7","FastN_AD7",50,0,100);
    a_FastN_Prompt[7] = new TH1D("FastN_AD8","FastN_AD8",50,0,100);
    /////////////////Fast Neutron End/////////////
    
    TTree* tree = new TTree("IBD","ibd event");
    TTree* tree2 = new TTree("Singles","singles event");
    TTree* tree3 = new TTree("amc","Information of Amc Singles rate per every day");
    TTree* tree4 = new TTree("live","liveTime for a run");
    tree->Branch("ep",&ep);
    tree->Branch("ed",&ed);
    tree->Branch("xp",&xp);
    tree->Branch("yp",&yp);
    tree->Branch("zp",&zp);
    tree->Branch("xd",&xd);
    tree->Branch("yd",&yd);
    tree->Branch("zd",&zd);
    tree->Branch("dist",&dist);
    tree->Branch("dT",&dT);
    tree->Branch("detNo",&detNo);

    tree2->Branch("E",&E);
    tree2->Branch("X",&X);
    tree2->Branch("Y",&Y);
    tree2->Branch("Z",&Z);
    tree2->Branch("T",&T);
    tree2->Branch("detNo",&sdetNo);
    tree2->Branch("isVeto",&isVeto);

    tree3->Branch("N",&AmcCount);
    tree3->Branch("N07",&N_07);
    tree3->Branch("N15",&N_15);
    tree3->Branch("N_up",&N_up);
    tree3->Branch("N_down",&N_down);
    tree3->Branch("Date",&Date);
    tree3->Branch("liveTime",&liveTimeCount);
    tree3->Branch("detNo",&detNo);
    tree3->Branch("Nmuon",&Nmuon);
    tree3->Branch("Nacc",&NaccSingle);

    tree4->Branch("liveTime",&theLiveTime);
    tree4->Branch("T",&theAbsRunTime);
    tree4->Branch("Date",&Date);
    tree4->Branch("detNo",&theDetNo);
    tree4->Branch("totalTime",&totalTime);
    tree4->Branch("NIBD",&NIBD);
    tree4->Branch("multiEff",&multiEff);
    tree4->Branch("multiEff_Err",&multiEff_Err);
    
    long long totalEventNumber  = adrec.fChain->GetEntries();
    long long thisEventTime     = 0;
    long long startTime         = 0;
    long long lastTime          = 0;
    long long MuonVetoEnd[5]    ={0,0,0,0,0};
    long long FastNVetoEnd[5]    ={0,0,0,0,0};
    long long liveTime[5]       ={0,0,0,0,0};
    long long preLiveTime[5]    ={0,0,0,0,0};
    long long AmcSingles[5]     ={0,0,0,0,0};
    long long Nup[5]            ={0,0,0,0,0};
    long long Ndown[5]          ={0,0,0,0,0};
    long long N07[5]            ={0,0,0,0,0};
    long long N15[5]            ={0,0,0,0,0};
    int day                     =1;
    long long closestPromptTime[5]={0,0,0,0,0};
    long long closestFastNPromptTime[5]={0,0,0,0,0};
    long long closestMuonTime[5]={0,0,0,0,0};
    long long closestFastNMuonTime[5]={0,0,0,0,0};
    long long closestWPMuonTime=0;
    long long closestOWSMuonTime=0;
    long long closestIWSMuonTime=0;
    long long preSingleTime[5]={0,0,0,0,0};

    bool sampleLock[5]={0,0,0,0,0};
    int stage[5] = {0,0,0,0,0};
    long long samp_pre_prompt_time[5]={0,0,0,0,0};
    long long samp_delayed_time[5]={0,0,0,0,0};

    for(long long i=0;i<totalEventNumber;i++){
        adrec.GetEntry(i);
        calst.GetEntry(i);
        PhyEvent phyEvent(calst,adrec);

        int AdNo = phyEvent.detector;
        if(EH2){
            if(AdNo==1){detNo=3;sdetNo=3;}
            if(AdNo==2){detNo=8;sdetNo=8;}
        }else{
            detNo=AdNo+EH;
            sdetNo=AdNo+EH;
        }

        if(i==0){
            thisEventTime   = phyEvent.t;
            startTime       = phyEvent.t;
            for(int j=1;j<5;j++){
                MuonVetoEnd[j]  = startTime;
                preSingleTime[j] = startTime;
                closestPromptTime[j] = 1000;
            }
        }
        else if(i==totalEventNumber-1){         //This else if is in purpose to deal with the end
            lastTime        = phyEvent.t;       //of entries, to fill in the live time
            cout<<lastTime-startTime<<endl;     //(correlated to efficient) and also the total
            totalTime       = lastTime-startTime;
            theAbsRunTime=phyEvent.trigTime_s;  // running time.
            for(int j=1;j<5;j++){
                if(EH1){
                    if(j>2) continue;
                    theDetNo=j;
                    theLiveTime=liveTime[j];
                    NIBD=a_NIBD[j];
                    multiEff=100.*double(Nsuccess[j])/double(Nsample[j]);
                    multiEff_Err=100.*sqrt(double(Nsuccess[j]))/double(Nsample[j]);
                    Date        =phyEvent.trigTime_s;
                    tree4->Fill();
                }
                if(EH2){
                    if(j==1) theDetNo=3;
                    if(j==2) theDetNo=8;
                    if(j>2) continue;
                    theLiveTime=liveTime[j];
                    NIBD=a_NIBD[j];
                    multiEff=100.*double(Nsuccess[j])/double(Nsample[j]);
                    multiEff_Err=100.*sqrt(double(Nsuccess[j]))/double(Nsample[j]);
                    Date        =phyEvent.trigTime_s;
                    tree4->Fill();
                }
                if(EH3){
                    theDetNo=j+EH;
                    theLiveTime=liveTime[j];
                    NIBD=a_NIBD[j];
                    multiEff=100.*double(Nsuccess[j])/double(Nsample[j]);
                    multiEff_Err=100.*sqrt(double(Nsuccess[j]))/double(Nsample[j]);
                    Date        =phyEvent.trigTime_s;
                    tree4->Fill();
                }
            }
        }else thisEventTime      = phyEvent.t;
        
        bool isMuon = false;
        if(!phyEvent.isGood) continue;
        if(phyEvent.isWPMu || phyEvent.isSHMu || phyEvent.isADMu) isMuon = true;
        if(isMuon){
            if(phyEvent.isADMu||phyEvent.isSHMu){
                if(phyEvent.isSHMu){
                    //if(thisEventTime-closestMuonTime[AdNo]>1000000)
                        a_tlSH_high_tmp[AdNo] = thisEventTime;
                }else if(thisEventTime-closestMuonTime[AdNo]>1000000){
                    if(phyEvent.isMulow)
                        a_tlSH_low_tmp[AdNo] = thisEventTime;
                    else if(phyEvent.isMumid)
                        a_tlSH_mid_tmp[AdNo] = thisEventTime;
                }
                closestMuonTime[AdNo]=thisEventTime;
            }
            if(phyEvent.isWPMu){
                closestWPMuonTime = thisEventTime;
                for(int j=1;j<5;j++){
                    if(IBDReady[j]==true){
                        if(thisEventTime-preVetoOfWP<a_ibd_td[j]){
                            IBDReady[j]=false;
                        }
                    }
                    if(thisEventTime-preVetoOfWP>MuonVetoEnd[j]) liveTime[j]+=thisEventTime-preVetoOfWP-MuonVetoEnd[j];
                    MuonVetoEnd[j]      = max(MuonVetoEnd[j],thisEventTime+pVetoOfWP);
                    NmuonArray[j]++;
                }
                //if(phyEvent.isIWSMu){
                //    if(thisEventTime-closestOWSMuonTime>2000){
                //        for(int j =1; j<5;j++){
                //            FastNVetoEnd[j]     = max(FastNVetoEnd[j],thisEventTime+pVetoOfWP);
                //        }
                //        closestIWSMuonTime = thisEventTime;
                //    }
                //}else/* if(phyEvent.isOWSMu)*/{
                //    closestOWSMuonTime = thisEventTime;
                //}
                continue;
            }else if(phyEvent.isADMu){
                if(thisEventTime>MuonVetoEnd[AdNo]) liveTime[AdNo]+=thisEventTime-MuonVetoEnd[AdNo];
                MuonVetoEnd[AdNo]   = max(MuonVetoEnd[AdNo],thisEventTime+pVetoOfAD);
                if(thisEventTime-closestOWSMuonTime>2000){
                    FastNVetoEnd[AdNo]  = max(FastNVetoEnd[AdNo],thisEventTime+pVetoOfAD);
                }
                NmuonArray[AdNo]++;
                continue;
            }else/* if(phyEvent.isSHMu)*/{
                epLock[AdNo]=false;
                if(thisEventTime>MuonVetoEnd[AdNo]) liveTime[AdNo]+=thisEventTime-MuonVetoEnd[AdNo];
                //MuonVetoEnd[AdNo]   = max(MuonVetoEnd[AdNo],thisEventTime+pVetoOfSH);
                FastNVetoEnd[AdNo]  = max(FastNVetoEnd[AdNo],thisEventTime+pVetoOfSH);
                NmuonArray[AdNo]++;
                continue;
            }
        }
        a_tlSH_high[AdNo] = a_tlSH_high_tmp[AdNo];
        if(thisEventTime-closestMuonTime[AdNo]>1000000){
            a_tlSH_mid[AdNo] = a_tlSH_mid_tmp[AdNo];
            a_tlSH_low[AdNo] = a_tlSH_low_tmp[AdNo];
        }
        ///////////AmcSingles rate counting.//////////
        if(thisEventTime-startTime>86400000000000*day){
            for(int j=1;j<5;j++){
                N_07        =N07[j];
                N_15        =N15[j];
                N_up        =Nup[j];
                N_down      =Ndown[j];
                AmcCount    =AmcSingles[j];
                Date        =phyEvent.trigTime_s;
                Nmuon       =NmuonArray[j];
                NaccSingle  =NaccSingleArr[j];
                if(MuonVetoEnd[j]>thisEventTime){
                    liveTimeCount=liveTime[j]-preLiveTime[j];
                    preLiveTime[j]=liveTime[j];
                }else{
                    liveTimeCount=liveTime[j]-preLiveTime[j]+thisEventTime-MuonVetoEnd[j];
                    preLiveTime[j]=liveTime[j]+thisEventTime-MuonVetoEnd[j];
                }
                if(EH2){
                    if(j==1){detNo=3;}
                    if(j==2){detNo=8;}
                }else{
                    detNo=j+EH;
                }
                if(N_up!=0){
                    tree3->Fill();
                }
                AmcSingles[j]=0;
                N07[j]=0;
                N15[j]=0;
                Nup[j]=0;
                Ndown[j]=0;
                NmuonArray[j]=0;
                NaccSingleArr[j]=0;
            }
                day+=1;
                detNo=sdetNo;
        }
        //////////////IBD candidates seletion start here////////////
        //if(phyEvent.isPrompt){
        //    if(thisEventTime-preSingleTime[AdNo]>1000&&thisEventTime>MuonVetoEnd[AdNo]){
        //        E=phyEvent.e;
        //        X=phyEvent.x;
        //        Y=phyEvent.y;
        //        Z=phyEvent.z;
        //        T=phyEvent.t;
        //        tree2->Fill();
        //        NaccSingleArr[AdNo]++;
        //        if(E>0.7){N07[AdNo]++;}
        //        if(E>1.5){N15[AdNo]++;}
        //        if(E>6.&&Z>0.){Nup[AdNo]++;AmcSingles[AdNo]++;}     //here two if are for
        //        if(E>6.&&Z<0.){Ndown[AdNo]++;AmcSingles[AdNo]--;}   // amc singles.
        //    }
        //    preSingleTime[AdNo]=thisEventTime;
        //}
        if(phyEvent.isPrompt){
            if(IBDReady[AdNo]){
                if(phyEvent.isDelay){
                    IBDReady[AdNo]=false;
                    if(thisEventTime-a_ibd_td[AdNo]>timeCut){
                        ep=a_ibd_ep[AdNo]; 
                        xp=a_ibd_xp[AdNo];
                        yp=a_ibd_yp[AdNo];
                        zp=a_ibd_zp[AdNo];
                        tp=a_ibd_tp[AdNo];
                        ed=a_ibd_ed[AdNo];
                        xd=a_ibd_xd[AdNo];
                        yd=a_ibd_yd[AdNo];
                        zd=a_ibd_zd[AdNo];
                        td=a_ibd_td[AdNo];
                        dtlSH_high = a_ibd_dtlSH_high[AdNo];
                        dtlSH_mid = a_ibd_dtlSH_mid[AdNo];
                        dtlSH_low = a_ibd_dtlSH_low[AdNo];
                        dT=td-tp;
                        dist=sqrt(pow(xd-xp,2)+pow(yd-yp,2)+pow(zd-zp,2));                    
                        Fill_dtlSH(array_h,dtlSH_high,dtlSH_mid,dtlSH_low,ep,ed,dT,dist);
                        //tree->Fill();
                        a_NIBD[AdNo]++;
                    }
                }
            }
            if(epLock[AdNo]==false){
                epLock[AdNo]=true;
                a_ep[AdNo]=phyEvent.e;
                a_xp[AdNo]=phyEvent.x;
                a_yp[AdNo]=phyEvent.y;
                a_zp[AdNo]=phyEvent.z;
                a_tp[AdNo]=phyEvent.t;
            }else{
                if(thisEventTime<MuonVetoEnd[AdNo]
                        ||thisEventTime-a_tp[AdNo]>timeCut
                        ||thisEventTime-a_tp[AdNo]<1000
                        ||phyEvent.isDelay!=true
                        ||thisEventTime-closestPromptTime[AdNo]<2*timeCut
                        ||pow(a_xp[AdNo]-phyEvent.x,2)+pow(a_yp[AdNo]-phyEvent.y,2)+pow(a_zp[AdNo]-phyEvent.z,2)>25.e4){
                    closestPromptTime[AdNo]=a_tp[AdNo];
                    a_ep[AdNo]=phyEvent.e;
                    a_xp[AdNo]=phyEvent.x;
                    a_yp[AdNo]=phyEvent.y;
                    a_zp[AdNo]=phyEvent.z;
                    a_tp[AdNo]=phyEvent.t;
                }else{
                    a_ibd_ep[AdNo]=a_ep[AdNo];
                    a_ibd_xp[AdNo]=a_xp[AdNo];
                    a_ibd_yp[AdNo]=a_yp[AdNo];
                    a_ibd_zp[AdNo]=a_zp[AdNo];
                    a_ibd_tp[AdNo]=a_tp[AdNo];
                    a_ibd_ed[AdNo]=phyEvent.e;
                    a_ibd_xd[AdNo]=phyEvent.x;
                    a_ibd_yd[AdNo]=phyEvent.y;
                    a_ibd_zd[AdNo]=phyEvent.z;
                    a_ibd_td[AdNo]=phyEvent.t;
                    a_ibd_dtlSH_high[AdNo]=double(a_ibd_tp[AdNo]-a_tlSH_high[AdNo]);
                    a_ibd_dtlSH_mid[AdNo]=double(a_ibd_tp[AdNo]-a_tlSH_mid[AdNo]);
                    a_ibd_dtlSH_low[AdNo]=double(a_ibd_tp[AdNo]-a_tlSH_low[AdNo]);
                    IBDReady[AdNo]=true;
                    epLock[AdNo]=false;
                    closestPromptTime[AdNo]=thisEventTime;
                }
            }
        }
    }
    h_high->Write();
    h_mid->Write();
    h_low->Write();
    output->Write();
}

void Fill_dtlSH(TH1D* h[5][3][3][4][8], double dt_high, double dt_mid, double dt_low, double prompt, double delay, double dt, double distance){
    int site = EH1*1+EH2*2+EH3*3-1;
    if(dt_high/1.e6<5.e3)
        h_high->Fill(dt_high/1.e6);
    h_mid->Fill(dt_mid/1.e6);
    h_low->Fill(dt_low/1.e6);
    for(int nslice = 4;nslice<9;nslice++){
        int slice_prompt = int((prompt-0.7)/(12.001-0.7)*1.*nslice);
        int slice_delay = int((delay-1.5)/(12.001-1.5)*1.*nslice);
        int slice_time = int((dt-1000.)/(400001.-1000.)*1.*nslice);
        int slice_distance = int((distance*1.)/(500)*1.*nslice);
        h[nslice-4][site][0][0][slice_prompt]->    Fill(dt_low/1.e6);
        h[nslice-4][site][0][1][slice_delay]->     Fill(dt_low/1.e6);
        h[nslice-4][site][0][2][slice_time]->      Fill(dt_low/1.e6);
        h[nslice-4][site][0][3][slice_distance]->  Fill(dt_low/1.e6);
        h[nslice-4][site][1][0][slice_prompt]->    Fill(dt_mid/1.e6);
        h[nslice-4][site][1][1][slice_delay]->     Fill(dt_mid/1.e6);
        h[nslice-4][site][1][2][slice_time]->      Fill(dt_mid/1.e6);
        h[nslice-4][site][1][3][slice_distance]->  Fill(dt_mid/1.e6);
        h[nslice-4][site][2][0][slice_prompt]->    Fill(dt_high/1.e6);
        h[nslice-4][site][2][1][slice_delay]->     Fill(dt_high/1.e6);
        h[nslice-4][site][2][2][slice_time]->      Fill(dt_high/1.e6);
        h[nslice-4][site][2][3][slice_distance]->  Fill(dt_high/1.e6);
    }
}
