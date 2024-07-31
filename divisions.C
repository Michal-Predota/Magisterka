#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

using namespace std;

int divisions()
{
	TFile *sameevent_HAL = new TFile("./files/file_nodelta_flow_010.root", "read");
	TFile *background_HAL = new TFile("./files/file_nodelta_noflow_010.root", "read");
	
	TFile *pimp_file = new TFile("everything_fullstat_pim_c010.root", "read");
	TFile *reco_file = new TFile("everything_pim_reco_c010.root", "read");
	
	
	TFile *pimp_file_norp = new TFile("everything_fullstat_pim_norp_c010.root", "read");
	TFile *reco_file_norp = new TFile("everything_pim_sim_reco_norp.root", "read");
	
	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	
	
	
	TCanvas canvas("canvas");
	canvas.Print("./outputs/bckg rp div no rp.pdf["); 
	
	
	for(int ynum=0; ynum<8; ynum++)
	{		
		for(int ptnum=0; ptnum<8; ptnum++)
		{
			char histName[100];	
			
			sprintf(histName, "same-event %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D *sameevent_exp = (TH1D*)pimp_file->Get(histName);
			TH1D *sameevent_exp_norp = (TH1D*)pimp_file_norp->Get(histName);
			
			
			sprintf(histName, "bckg %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* bckg_exp = (TH1D*)pimp_file->Get(histName);
			TH1D* bckg_exp_norp = (TH1D*)pimp_file_norp->Get(histName);
			
			
			sprintf(histName, "same-event %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D *sameevent_MC_hist = (TH1D*)reco_file->Get(histName);
			
			
			sprintf(histName, "bckg %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* bckg_MC_hist = (TH1D*)reco_file->Get(histName);
			
			
			sprintf(histName, "M_{inv} %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D *sameevent_HAL_hist = (TH1D*)sameevent_HAL->Get(histName);
			TH1D* bckg_HAL_hist = (TH1D*)background_HAL->Get(histName);
			
			
			
			
			sameevent_exp->Rebin(10);
			bckg_exp->Rebin(10);
			sameevent_MC_hist->Rebin(10);
			bckg_MC_hist->Rebin(10);
			sameevent_HAL_hist->Rebin(10);
			bckg_HAL_hist->Rebin(10);
			bckg_exp_norp->Rebin(10);
			sameevent_exp_norp->Rebin(10);
			
			
			sameevent_exp->Scale(1/sameevent_exp->Integral());
			bckg_exp->Scale(1/bckg_exp->Integral());
			bckg_exp_norp->Scale(1/bckg_exp_norp->Integral());
			sameevent_MC_hist->Scale(1/sameevent_MC_hist->Integral());
			bckg_MC_hist->Scale(1/bckg_MC_hist->Integral());
			sameevent_HAL_hist->Scale(1/sameevent_HAL_hist->Integral());
			bckg_HAL_hist->Scale(1/bckg_HAL_hist->Integral());
			sameevent_exp_norp->Scale(1/sameevent_exp_norp->Integral());
			
			
			
		//	sprintf(histName, "HAL %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D *exp_divided = new TH1D(*sameevent_exp);
			TH1D *HAL_divided = new TH1D(*sameevent_HAL_hist);
			TH1D *MC_divided = new TH1D(*sameevent_MC_hist);
			TH1D *sameevent_exp_divided = new TH1D(*sameevent_exp);
			
			
	
			exp_divided->Divide(bckg_exp);
			
			bckg_exp->Divide(bckg_exp_norp);
			
			HAL_divided->Divide(bckg_HAL_hist);
	
			MC_divided->Divide(bckg_MC_hist);
			
			sameevent_exp_divided->Divide(sameevent_exp_norp);
			
		//	exp_divided->Divide(MC_divided);
		//	HAL_divided->Divide(MC_divided);
		//	exp_divided->Divide(HAL_divided);
			
	
			
			sprintf(histName, "exp/HAL %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			exp_divided->GetYaxis()->SetRangeUser(0.8,1.2);
			exp_divided->GetXaxis()->SetRangeUser(1000,1600);
			exp_divided->SetNameTitle(histName,histName);
			exp_divided->SetStats(0);
			
			
			sprintf(histName, "HAL/MC %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			HAL_divided->SetNameTitle(histName,histName);
			HAL_divided->SetLineColor(kRed);
		//	HAL_divided->GetYaxis()->SetRangeUser(0.95,1.1);
			
			
			sprintf(histName, "bckg rp/no rp %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			bckg_exp->SetNameTitle(histName,histName);
			bckg_exp->GetYaxis()->SetRangeUser(0.95,1.05);
			
			sprintf(histName, "sameevent exp rp/no rp %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			sameevent_exp_divided->SetNameTitle(histName,histName);
			sameevent_exp_divided->GetYaxis()->SetRangeUser(0.95,1.05);
			
			canvas.Clear();
		//	bckg_exp->Draw();
		
			bckg_exp->Draw();
			
		//	exp_divided->Draw();
		//	HAL_divided->Draw("same");
		//	TLegend *legend = new TLegend(0.7,0.9,0.9,0.8);
		//	legend->AddEntry(exp_divided, "exp/MC", ",l");
		//	legend->AddEntry(HAL_divided, "HAL/MC", ",l");
		//	legend->AddEntry(MC_divided, "MC", ",l");
			
		//	legend->Draw();
		
			canvas.Print("./outputs/bckg rp div no rp.pdf");   
		}
	}
	
	
	canvas.Print("./outputs/bckg rp div no rp.pdf]");   
	
	return 0;
}