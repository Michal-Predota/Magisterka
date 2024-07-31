#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

using namespace std;

int subtract()
{
	TFile *pimp_file = new TFile("everything_fullstat_pim_c010.root", "read");
	TFile *pimp_file_1 = new TFile("everything_fullstat_pim_c010.root", "read");
	TFile *pimp_file_norp = new TFile("everything_fullstat_pim_norp_c010.root", "read");
	
	TFile *output = new TFile("./files/subtracted.root", "recreate");
	
	
	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	
	
	TCanvas canvas("canvas");
	canvas.Print("./outputs/subtracted.pdf["); 
	
		
	for(int ynum=0; ynum<8; ynum++)
	{		
		for(int ptnum=0; ptnum<8; ptnum++)
		{
			
			char histName[100];	
			
			sprintf(histName, "same-event %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D *sameevent_exp = (TH1D*)pimp_file->Get(histName);
			
			sprintf(histName, "bckg %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* bckg_exp = (TH1D*)pimp_file->Get(histName);
			TH1D* bckg_exp_norp = (TH1D*)pimp_file_norp->Get(histName);
			TH1D* bckg_exp_copy = (TH1D*)pimp_file_1->Get(histName);
			
			
			
			
			bckg_exp_norp->Rebin(4);
			bckg_exp_copy->Rebin(4);
			sameevent_exp->Rebin(4);
			bckg_exp->Rebin(4);
			
		
			bckg_exp->Scale(sameevent_exp->GetEntries()/bckg_exp->GetEntries(), "nosw2");
				
			
			TH1D *exp_subtracted = new TH1D(*sameevent_exp);
			
			sprintf(histName, "same-mixed %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			exp_subtracted->SetNameTitle(histName,histName);
			
			exp_subtracted->Add(bckg_exp, -1);
			bckg_exp->SetLineColor(kRed);
			exp_subtracted->SetStats(0);
			exp_subtracted->GetXaxis()->SetRangeUser(1000,2000);
			exp_subtracted->GetXaxis()->SetRangeUser(1000,2000);
			
			
			
			
			
			canvas.Clear();
			
			TLegend *legend = new TLegend(0.7,0.9,0.9,0.8);
			legend->AddEntry(exp_subtracted, "same - mixed rp", ",l");
			exp_subtracted->Draw("");
			
			legend->Draw();
			canvas.Print("./outputs/subtracted.pdf");   
			
			
			output->cd();
			exp_subtracted->Write();
			
		}
	}
	
	output->Save();
	canvas.Print("./outputs/subtracted.pdf]");   
	
	
	return 0;
}