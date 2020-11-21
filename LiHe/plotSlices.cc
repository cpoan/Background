void plotSlices(){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    double dailyRates[7][3][2] = {
        {//4
            {5.58,1.240},
            {4.19,1.174},
            {0.39,0.043}
        },
        {//5
            {4.99,1.083},
            {4.69,0.925},
            {0.38,0.046}
        },
        {//6
            {6.00,1.466},
            {4.11,0.942},
            {0.40,0.047}
        },
        {//7
            {5.35,1.109},
            {4.38,0.901},
            {0.37,0.045}
        },
        {//8
            {5.62,1.365},
            {4.28,0.796},
            {0.39,0.045}
        },
        {//9
            {5.84,1.296},
            {4.75,0.916},
            {0.38,0.044}
        },
        {//10
            {4.82,1.083},
            {4.02,0.699},
            {0.39,0.037}
        }
    };
    TGraphErrors* gr[3];
    for(int eh = 0;eh<3;eh++){
        gr[eh] = new TGraphErrors();
        gr[eh]->SetMarkerStyle(8);
        gr[eh]->SetMarkerSize(1.8);
        gr[eh]->SetLineWidth(2);
        if(eh==0){
            gr[eh]->SetMarkerColor(kRed);
            gr[eh]->SetLineColor(kRed);
        }else if(eh==1){
            gr[eh]->SetMarkerColor(kBlue);
            gr[eh]->SetLineColor(kBlue);
        }else{
            gr[eh]->SetMarkerColor(kGreen);
            gr[eh]->SetLineColor(kGreen);
        }
        for(int slice = 0;slice<7;slice++){
            gr[eh]->SetPoint(slice,slice+4,dailyRates[slice][eh][0]);
            gr[eh]->SetPointError(slice,0,dailyRates[slice][eh][1]);
        }
    }
    TGraphErrors* ref_gr = new TGraphErrors();
    ref_gr->SetPoint(0,3.5,-1);
    ref_gr->SetPoint(1,10.5,8);
    ref_gr->SetMarkerStyle(8);
    ref_gr->SetMarkerSize(0.01);
    ref_gr->GetXaxis()->CenterTitle(kTRUE);
    ref_gr->GetYaxis()->CenterTitle(kTRUE);
    ref_gr->SetTitle("Daily Rate R_{LiHe};Slices;R[day^{-1}]");
    TCanvas* c = new TCanvas("c","c",1600,900);
    c->cd();
    gPad->SetGrid();
    ref_gr->Draw("AP");
    for(int i = 0 ;i<3;i++){
        gr[i]->Draw("pSAME");
    }
    TLegend* lg = new TLegend(0.75,0.75,0.95,0.95);
    lg->AddEntry(gr[0],"EH1","pl");
    lg->AddEntry(gr[1],"EH2","pl");
    lg->AddEntry(gr[2],"EH3","pl");
    lg->Draw("SAME");
    TCanvas* c2 = new TCanvas("c2","c2",3000,500);
    c2->Divide(3,1);
    for(int i = 0;i<3;i++){
        c2->cd(i+1);
        gPad->SetGrid();
        gr[i]->Fit("pol0","Q");
        gr[i]->Draw("ap");
        if(i==0)
            gr[i]->GetFunction("pol0")->SetLineColor(kRed);
        else if(i==1)
            gr[i]->GetFunction("pol0")->SetLineColor(kBlue);
        else
            gr[i]->GetFunction("pol0")->SetLineColor(kGreen);
    }
    for(int eh = 0;eh<3;eh++){
        double rateError = 0;
        for(int slice = 0;slice<7;slice++){
            rateError+=abs(dailyRates[slice][eh][1])/7.;
        }
        cout << "EH" << eh+1 << " " << rateError << "," << "SliceVariation = " << gr[eh]->GetRMS(2) << "\n";
    }
}
