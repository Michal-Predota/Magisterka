#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"


int temperature()
{
	TFile *file = new TFile("./pty_pip_fullstat.root","read");
	TFile *output_temp = new TFile("./files/temp_010_pip.root","recreate");
	
	
	TH2D *pty_c010 = (TH2D*)file->Get("ptvsy_pair_c010");
	TH2D *pty_c1020 = (TH2D*)file->Get("ptvsy_pair_c1020");
	TH2D *pty_c2030 = (TH2D*)file->Get("ptvsy_pair_c2030");
	TH2D *pty_c3040 = (TH2D*)file->Get("ptvsy_pair_c3040");
	
	
	
	
	
	TH1D *pt_0y02 = new TH1D("p_T 0<y<0.2","p_{T} 0<y<0.2",100,0,2000);
	TH1D *pt_02y04 = new TH1D("p_T 0.2<y<0.4","p_{T} 0.2<y<0.4",100,0,2000);
	TH1D *pt_04y06 = new TH1D("p_T 0.4<y<0.6","p_{T} 0.4<y<0.6",100,0,2000);
	TH1D *pt_06y08 = new TH1D("p_T 0.6<y<0.8","p_{T} 0.6<y<0.8",100,0,2000);
	TH1D *pt_08y10 = new TH1D("p_T 0.8<y<1.0","p_{T} 0.8<y<1.0",100,0,2000);
	TH1D *pt_10y12 = new TH1D("p_T 1.0<y<1.2","p_{T} 1.0<y<1.2",100,0,2000);
	TH1D *pt_12y14 = new TH1D("p_T 1.2<y<1.4","p_{T} 1.2<y<1.4",100,0,2000);
	TH1D *pt_14y16 = new TH1D("p_T 1.4<y<1.6","p_{T} 1.4<y<1.6",100,0,2000);
	TH1D *pt_16y18 = new TH1D("p_T 1.6<y<1.8","p_{T} 1.6<y<1.8",100,0,2000);
	TH1D *pt_18y20 = new TH1D("p_T 1.8<y<2.0","p_{T} 1.8<y<2.0",100,0,2000);
	
	TH2D *test = new TH2D("h","h",10,0,1,10,0,1);
	
	
	
	
	TH1D *pt = new TH1D("pt","pt",100,0,2000);
	
	vector<TH1D*> histos = {pt_0y02,pt_02y04,pt_04y06,pt_06y08,pt_08y10,pt_10y12,pt_12y14,pt_14y16,pt_16y18,pt_18y20};
	
	
	
	int values[10][100];
	
	for(int i=0; i<10; i++)
	{
		for(int j=0; j<100; j++)
		{
			values[i][j]=pty_c010->GetBinContent(pty_c010->GetBin(i+1,j+1));
		}
	}
	
	TF1* fit002 = new TF1("fit002","[0]*exp(-[1]*x)");
	TF1* fit0204 = new TF1("fit0204","[0]*exp(-[1]*x)");
	TF1* fit0406 = new TF1("fit0406","[0]*exp(-[1]*x)");
	TF1* fit0608 = new TF1("fit0608","[0]*exp(-[1]*x)");
	TF1* fit0810 = new TF1("fit0810","[0]*exp(-[1]*x)");
	TF1* fit1012 = new TF1("fit1012","[0]*exp(-[1]*x)");
	TF1* fit1214 = new TF1("fit1214","[0]*exp(-[1]*x)");
	TF1* fit1416 = new TF1("fit1416","[0]*exp(-[1]*x)");
	TF1* fit1618 = new TF1("fit1618","[0]*exp(-[1]*x)");
	TF1* fit1820 = new TF1("fit1820","[0]*exp(-[1]*x)");
	
	vector<TF1*> fits = {fit002,fit0204,fit0406,fit0608,fit0810,fit1012,fit1214,fit1416,fit1618,fit1820};
	
	
	
	
	TCanvas canvas("canvas");
	canvas.Print("./outputs/temperature fits.pdf["); 
	
	canvas.SetLogy();
	
	
	for(int i=0; i<8; i++)
	{
		for(int j=0; j<100; j++)
		{
			histos.at(i)->SetBinContent(j+1, values[i][j]);
		}
		
		fits.at(i)->SetRange(1300, 2000);
		fits.at(i)->SetParameter(1,1/150);
		fits.at(i)->SetParameter(0,1);
		
		
		
		double l,h;
		fits.at(i)->GetRange(l,h);
		//cout<<l<<" "<<h<<endl;
		histos.at(i)->Sumw2();
		histos.at(i)->SetMarkerStyle(8);
		histos.at(i)->SetMarkerSize(0.5);
		histos.at(i)->SetMarkerColor(histos.at(i)->GetLineColor());
		histos.at(i)->Fit(fits.at(i), "RE");
		
		cout<<"probability: "<<fits.at(i)->GetProb()<<endl;
		cout<<"temperature:  "<<1/fits.at(i)->GetParameter(1)<<endl;
		
		
		canvas.Clear();
	//	canvas.SetLogy();
		
		TPaveText text(1600,0.4*histos.at(i)->GetMaximum(),2000,1.9*histos.at(i)->GetMaximum());
		char temp[100];
		sprintf(temp, "T_{eff} = %f", 1/fits.at(i)->GetParameter(1));
		text.AddText(temp);
		text.SetBorderSize(1);
		
		fits.at(i)->SetLineWidth(2);
		
		histos.at(i)->GetXaxis()->SetTitle("p_{T} (MeV/c)");
		histos.at(i)->GetYaxis()->SetTitle("dN/dp_{T}");
		
		histos.at(i)->GetYaxis()->SetTitleSize(0.05);
		histos.at(i)->GetXaxis()->SetTitleSize(0.05);
		histos.at(i)->GetYaxis()->SetLabelSize(0.04);
		histos.at(i)->GetXaxis()->SetLabelSize(0.04);
		
		
		histos.at(i)->SetStats(0);
		histos.at(i)->Draw();
		text.Draw();
		
		
		
		gPad->SetBottomMargin(0.12);
		gPad->SetTopMargin(0.08);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.12);
		
		
		const char *name;
		name = histos.at(i)->GetName();
		
		char fileName[100];
		sprintf(fileName, "./images/pip/temperature/%s_010.png",name);
		
		canvas.SaveAs(fileName);
		
		canvas.Print("./outputs/temperature fits.pdf"); 
	}
	TCanvas *c = new TCanvas("c","c",800,800);
	c->cd();
//	histos.at(0)->Draw();
	
	
	TH1D *tmp_vs_y = new TH1D("#Delta^{0}","#Delta^{0}", 8,0,1.6);
	
	
	for(int i=0; i<8; i++)
	{
		cout<<1/fits.at(i)->GetParameter(1)<<endl;
		tmp_vs_y->SetBinContent(i+1, 1/fits.at(i)->GetParameter(1));
		tmp_vs_y->SetBinError(i+1,1/fits.at(i)->GetParameter(1)*fits.at(i)->GetParError(1));
		tmp_vs_y->SetMarkerStyle(8);
		tmp_vs_y->SetMarkerSize(0.5);
		
		tmp_vs_y->GetXaxis()->SetTitle("y");
		tmp_vs_y->GetYaxis()->SetTitle("T_{eff}");
		
	}
	
		output_temp->cd();
		tmp_vs_y->Write();
	
	canvas.Clear();
	tmp_vs_y->Draw();
	canvas.Print("./outputs/temperature fits.pdf"); 
	
	
	

	
	canvas.Print("./outputs/temperature fits.pdf]"); 
	
	output_temp->Save();
	

//	test->Draw("text");
	
	
	
	
	
	
	return 0;
}