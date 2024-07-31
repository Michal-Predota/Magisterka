#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"


using namespace std;

int exphalcomp()
{
	
	
	TFile *nodelta_file = new TFile("./files/file_nodelta_flow_010.root", "read");
	TFile *nodelta0f_file = new TFile("./files/file_nodelta_noflow_010.root", "read");
	TFile *sameevent_file = new TFile("sameevent_pim_010.root", "read");
	TFile *background_file = new TFile("background_pim_010_10k.root", "read");
	
	
	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	
	TH1D *plots[8][8];
	TCanvas* canvases[8];
	TCanvas canvas("canvas");
	canvas.Print("./outputs/HAL and exp pim.pdf["); 
	int filenum=0;
	for(int ynum=0; ynum<8; ynum++)
	{
		for(int ptnum=0; ptnum<8; ptnum++)
		{
			char histName[100];
			sprintf(histName, "M_{inv} %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			cout<<histName<<endl;
			TH1D* hist_nodelta = (TH1D*)nodelta_file->Get(histName);
			TH1D* hist_nodelta0f = (TH1D*)nodelta0f_file->Get(histName);
			
			
			
			hist_nodelta->Rebin(10);
			hist_nodelta0f->Rebin(10);
			
			hist_nodelta->Scale(1/hist_nodelta->Integral(), "nosw2");
			hist_nodelta0f->Scale(1/hist_nodelta0f->Integral(), "nosw2");
			
			TH1D* hist_divided = new TH1D(*hist_nodelta);
			hist_divided->Divide(hist_nodelta0f);
		
		
			
			sprintf(histName, "flow %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			hist_divided->SetNameTitle(histName,  histName);
			hist_divided->GetYaxis()->SetRangeUser(0, 2);
			
			
			
			
			sprintf(histName, "bckg %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* hist_bckg = (TH1D*)background_file->Get(histName);
			sprintf(histName, "same-event %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* hist_sameevent = (TH1D*)sameevent_file->Get(histName);
			
			hist_sameevent->Rebin(10);
			hist_bckg->Rebin(10);
			
			hist_sameevent->Scale(1/hist_sameevent->Integral(), "nosw2");
			hist_bckg->Scale(1/hist_bckg->Integral(), "nosw2");
			
			TH1D* hist_divided_exp = new TH1D(*hist_sameevent);
			
			hist_divided_exp->Divide(hist_bckg);
			hist_divided_exp->SetLineColor(kRed);
			//hist_divided_exp->GetYaxis()->SetRangeUser(0.98, 1.02);
			
			
			TLegend *legend = new TLegend(0.7,0.9,0.9,0.8);
			
			legend->AddEntry(hist_divided_exp, "experimental", ",l");
			legend->AddEntry(hist_divided, "HAL", "l");
			
			canvas.Clear();
			cout<<ynum<<" "<<ptnum<<endl;
			hist_divided->GetXaxis()->SetRangeUser(1000,1600);
			hist_divided->SetStats(0);
			hist_divided->Sumw2(false);
			hist_divided->Draw();  
			hist_divided_exp->GetXaxis()->SetRangeUser(1000,1600);
			hist_divided_exp->Sumw2(false);
			hist_divided_exp->Draw("SAME");
			legend->Draw();
			canvas.Print("./outputs/HAL and exp pim.pdf");    
		}
	}
	
	
	
	canvas.Print("./outputs/HAL and exp pim.pdf]");    
	return 0;
}





  