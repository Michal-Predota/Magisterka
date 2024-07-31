#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

using namespace std;

int samevsbckg()
{
	
	
	TFile *sameevent_file = new TFile("sameevent_pim_010_100k.root", "read");
	TFile *background_file = new TFile("background_pim_010_10k.root", "read");
	
	TFile *pimp_file = new TFile("everything_pim_010_norp.root", "read");
	TFile *smallMixing = new TFile("results_pim.root", "read");
	
	TFile *recoMixing = new TFile("everything_pim_sim_reco.root", "read");
	
	
	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	
	TCanvas canvas("canvas");
	canvas.Print("./outputs/same over bckg exp.pdf["); 
	for(int ynum=0; ynum<8; ynum++)
	{
		for(int ptnum=0; ptnum<8; ptnum++)
		{
			TLegend *legend = new TLegend(0.65,0.7,0.9,0.9);
			
			char histName[100];
			
			sprintf(histName, "bckg %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* hist_bckg = (TH1D*)pimp_file->Get(histName);
			hist_bckg->SetLineColor(kRed);
			legend->AddEntry(hist_bckg, "background", "l");
			
			sprintf(histName, "same-event %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* hist_sameevent = (TH1D*)pimp_file->Get(histName);
			legend->AddEntry(hist_sameevent, "same event", ",l");
			
			hist_sameevent->Rebin(10);
			hist_bckg->Rebin(10);
			
			hist_bckg->Scale(1/hist_bckg->Integral());
			hist_sameevent->Scale(1/hist_sameevent->Integral());
			
			sprintf(histName, "divided MC %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			
			hist_sameevent->Divide(hist_bckg);
			hist_sameevent->SetNameTitle(histName,histName);
			
			canvas.Clear();
			//gPad->SetLogy();
			hist_bckg->SetStats(0);
			//hist_bckg->Draw("hist");  
			hist_sameevent->GetYaxis()->SetRangeUser(0.95,1.05);
			hist_sameevent->Draw();
			legend->SetTextSize(0.05);
			//legend->Draw();
			canvas.Print("./outputs/same over bckg exp.pdf");    
			
			
		}
	}
	
	
//	TH1D* same_all = (TH1D*)sameevent_file->Get("#pi^{-} p");
//	TH1D* bckg_all = (TH1D*)background_file->Get("invmass");
	
	
//	same_all->Rebin(4);
//	bckg_all->Rebin(4);
	
	
//	same_all->Scale(1/same_all->Integral());
//	bckg_all->Scale(1/bckg_all->Integral());
	
//	same_all->Divide(bckg_all);
	
	
	
	
//	canvas.Clear();
//	gPad->SetLogy(0);
//	same_all->SetStats(0);
//	same_all->Sumw2(false);
//	same_all->GetYaxis()->SetRangeUser(0.5,1.5);
//	same_all->Draw("hist");
//	canvas.Print("./outputs/same event vs background pim.pdf");    
			
	
	
	
	canvas.Print("./outputs/same over bckg exp.pdf]");    
	
	
	return 0;
}