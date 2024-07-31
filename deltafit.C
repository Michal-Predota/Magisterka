#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"


int deltafit()
{
	
	TFile *file1 = new TFile("./files/delta++_3040.root","read");
	TFile *output = new TFile("./files/fits_d++_3040.root","recreate");
	
	//TF1* breit_wigner = new TF1("breit_wigner", "[0]*(x*x*0.117)/(TMath::Power((1.232*1.232-x*x), 2)+(x*x*0.117*0.117))", 0, 2);//breit wigner in HAL
	
	
	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	
	TCanvas canvas("canvas");
	canvas.Print("./outputs/delta++fit_3040.pdf["); 
	
	vector<double> mass;
	vector<double> width;
		
	TH2D *masses = new TH2D("masses","masses",8,0,1.6,8,0,1.6);
	masses->GetXaxis()->SetTitle("y");
	masses->GetYaxis()->SetTitle("pT");
	TH2D *widths = new TH2D("widths","widths",8,0,1.6,8,0,1.6);
	widths->GetXaxis()->SetTitle("y");
	widths->GetYaxis()->SetTitle("pT");
	
	
		
	TH2D *masses_errors = new TH2D("masses_errors","masses_errors",8,0,1.6,8,0,1.6);
	masses_errors->GetXaxis()->SetTitle("y");
	masses_errors->GetYaxis()->SetTitle("pT");
	TH2D *widths_errors = new TH2D("widths_errors","widths_errors",8,0,1.6,8,0,1.6);
	widths_errors->GetXaxis()->SetTitle("y");
	widths_errors->GetYaxis()->SetTitle("pT");
	TH2D *probability = new TH2D("probability","probability",8,0,1.6,8,0,1.6);
	probability->GetXaxis()->SetTitle("y");
	probability->GetYaxis()->SetTitle("pT");
	TH2D *chi2ndf = new TH2D("chi2ndf","chi2ndf",8,0,1.6,8,0,1.6);
	chi2ndf->GetXaxis()->SetTitle("y");
	chi2ndf->GetYaxis()->SetTitle("pT");
		
	for(int ynum=0; ynum<8; ynum++)
	{		
		for(int ptnum=0; ptnum<8; ptnum++)
		{
				char histName[100];
				sprintf(histName, "delta %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				TH1D *delta = (TH1D*)file1->Get(histName);
				
				
				for(int i=0; i<delta->GetNbinsX(); i++)
				{
					if(delta->GetBinContent(i+1)<0)
						delta->SetBinContent(i+1,0);
				}
				
			//	delta->Rebin(10);
				
				cout<<histName<<endl;
				char fitName[100];
				sprintf(fitName, "breit_wigner %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				TF1* breit_wigner = new TF1(fitName, "[0]*(x*x*[2])/(TMath::Power(([1]*[1]-x*x), 2)+(x*x*[2]*[2]))", 1120, 1225);//[0]-scale, [1]-mass, [2]-width
				
		//		breit_wigner->SetRange(0.9*delta->GetMean(), 1.1*delta->GetMean());
				
				breit_wigner->SetParameter(1, 1232);
				breit_wigner->SetParameter(2, 117);
				
				
				delta->Fit(breit_wigner, "RE");
				
				mass.push_back(breit_wigner->GetParameter(1));
				width.push_back(breit_wigner->GetParameter(2));
				
				cout<<breit_wigner->GetProb()<<endl;
				
				
				
				double massVal = breit_wigner->GetParameter(1);
				double widthVal = breit_wigner->GetParameter(2);
				
				double massErr = breit_wigner->GetParError(1);
				double widthErr = breit_wigner->GetParError(2);
				
				if(massVal>0 && massVal!=1232 && widthVal!=117 && widthVal>0 && massErr<0.1*massVal && widthErr<0.1*widthVal)
				{
					masses->SetBinContent(masses->GetBin(ynum+1,ptnum+1),breit_wigner->GetParameter(1));
					widths->SetBinContent(widths->GetBin(ynum+1,ptnum+1),breit_wigner->GetParameter(2));
				
					masses_errors->SetBinContent(masses_errors->GetBin(ynum+1,ptnum+1),breit_wigner->GetParError(1));
					widths_errors->SetBinContent(widths_errors->GetBin(ynum+1,ptnum+1),breit_wigner->GetParError(2));
					probability->SetBinContent(probability->GetBin(ynum+1,ptnum+1),breit_wigner->GetProb());
					chi2ndf->SetBinContent(chi2ndf->GetBin(ynum+1,ptnum+1),breit_wigner->GetChisquare()/breit_wigner->GetNDF());
				}
				
				canvas.Clear();
	//			canvas.SetLogy();
				delta->GetXaxis()->SetRangeUser(1000,2000);
				delta->Draw();
				canvas.Print("./outputs/delta++fit_3040.pdf"); 
				output->cd();
				breit_wigner->Write();
				
				cout<<endl<<endl;
		}
	}
	
	
	
	canvas.SetLogy(0);
	
	TH1D* hMass = new TH1D("hMass","hMass",100,1000,1300);
	
	for(int i=0;i<mass.size();i++)
	{
		
		hMass->Fill(mass.at(i));
	}
	
	canvas.Clear();
	hMass->Draw("CONT");
	canvas.Print("./outputs/delta++fit_3040.pdf"); 
	
	
	TH1D* hWidth = new TH1D("hWidth","hWidth",100,0,500);
	
	for(int i=0;i<mass.size();i++)
	{
		hWidth->Fill(width.at(i));
	}
	
	canvas.Clear();
	hWidth->Draw();
	canvas.Print("./outputs/delta++fit_3040.pdf"); 
	
	canvas.Clear();
	masses->SetStats(0);
	masses->Draw("colzTEXT");
	canvas.Print("./outputs/delta++fit_3040.pdf"); 
	
	canvas.Clear();
	widths->SetStats(0);
	widths->Draw("colzTEXT");
	canvas.Print("./outputs/delta++fit_3040.pdf"); 
	
	canvas.Clear();
	masses_errors->SetStats(0);
	masses_errors->Draw("colzTEXT");
	canvas.Print("./outputs/delta++fit_3040.pdf"); 
	
	canvas.Clear();
	widths_errors->SetStats(0);
	widths_errors->Draw("colzTEXT");
	canvas.Print("./outputs/delta++fit_3040.pdf"); 
	
	canvas.Clear();
	probability->SetStats(0);
	probability->Draw("colzTEXT");
	canvas.Print("./outputs/delta++fit_3040.pdf"); 
	
	canvas.Clear();
	chi2ndf->SetStats(0);
	chi2ndf->Draw("colzTEXT");
	canvas.Print("./outputs/delta++fit_3040.pdf"); 
	
	
	canvas.Print("./outputs/delta++fit_3040.pdf]"); 
	
	
	output->cd();
	masses->Write();
	widths->Write();
	masses_errors->Write();
	widths_errors->Write();
	probability->Write();
	chi2ndf->Write();
	
	
	output->Save();
	
	return 0;
}