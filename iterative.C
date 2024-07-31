#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"


double my_integral(TH1D* h, int limit)
{
	double integral=0;
	
	for(int i=1; i<=h->GetNbinsX(); i++)
	{
		if(h->GetBinContent(i)>0 && h->GetBinLowEdge(i)<=2000)
			integral+=h->GetBinContent(i);
	}
	
	return integral;
}


void iterative() 
{
	
	
	
	TF1 *function = new TF1("f","[0]+[1]*x",1000,2500);
	
	
		
//	for(int i=1;i<=20;i++)
	//{
	int	i=20;
	char fileName[100];	
	sprintf(fileName, "./outputs/iterative.pdf[");
	TCanvas canvas("canvas");
	canvas.Print(fileName); 
	
	TFile *pimp_file = new TFile("everything_fullstat_pip_c010.root", "read");
	TFile *pimp_file_copy = new TFile("everything_fullstat_pip_c010.root", "read");
	TFile *pimp_file_norp = new TFile("everything_fullstat_pip_norp_c010.root", "read");
	
	TFile *HALflow = new TFile("./files/file_pip_nodelta_flow_010.root", "read");
	TFile *HALnoflow = new TFile("./files/file_pip_nodelta_noflow_010.root", "read");
	TFile *HALdelta = new TFile("./files/onlydelta_file.root", "read");
	
	TFile *MC_file = new TFile("everything_pip_reco_c010.root","read");
	
	TFile *MC_file_norp = new TFile("everything_pip_reco_c010.root","read");
	
	TFile *output = new TFile("./files/delta++_010.root","recreate");
	
	//double weight = (double)i/1000;
	
	for(int ynum=0; ynum<8; ynum++)
	{		
		for(int ptnum=0; ptnum<8; ptnum++)
		{
			
	
	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
			
			
			
			char histName[100];	
			
			sprintf(histName, "same-event %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D *sameevent_exp = (TH1D*)pimp_file->Get(histName);
			cout<<histName<<endl;
			
			sprintf(histName, "bckg %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* mixed_exp = (TH1D*)pimp_file->Get(histName);
			TH1D* mixed_exp_copy = (TH1D*)pimp_file_copy->Get(histName);
			
			sprintf(histName, "bckg %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* mixed_exp_norp = (TH1D*)pimp_file_norp->Get(histName);
			TH1D* copy = (TH1D*)pimp_file_norp->Get(histName);
			
			sprintf(histName, "M_{inv} %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* flow_hist = (TH1D*)HALflow->Get(histName);
			
			sprintf(histName, "M_{inv} %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* noflow_hist = (TH1D*)HALnoflow->Get(histName);
			
			sprintf(histName, "M_{inv} %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* delta_hist = (TH1D*)HALdelta->Get(histName);
			
			sprintf(histName, "same-event %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D *sameevent_MC = (TH1D*)MC_file->Get(histName);
			
			sprintf(histName, "bckg %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* mixed_MC = (TH1D*)MC_file->Get(histName);
			
			
			sameevent_exp->Rebin(10);
			mixed_exp->Rebin(10);
			mixed_exp_copy->Rebin(10);
			mixed_exp_norp->Rebin(10);
			flow_hist->Rebin(10);
			noflow_hist->Rebin(10);
			delta_hist->Rebin(10);
			sameevent_MC->Rebin(10);
			mixed_MC->Rebin(10);
			
			
			sameevent_exp->Sumw2();
			mixed_exp->Sumw2();
			mixed_exp_copy->Sumw2();
			mixed_exp_norp->Sumw2();
			flow_hist->Sumw2();
			noflow_hist->Sumw2();
			delta_hist->Sumw2();
			sameevent_MC->Sumw2();
			mixed_MC->Sumw2();
			

			double weight = 0.04;
			
			mixed_exp->Scale(sameevent_exp->GetEntries()/mixed_exp->GetEntries(), "nosw2");
			
			TH1D* UDetectorFlow = flow_hist;
			UDetectorFlow->Divide(noflow_hist);
			UDetectorFlow->Multiply(mixed_exp);
		//	cout<<0.001*sameevent_exp->Integral()/UDetectorFlow->Integral()<<endl;
		//	UDetectorFlow->Scale(0.001*sameevent_exp->Integral()/UDetectorFlow->Integral());
			UDetectorFlow->Scale(0.001*my_integral(sameevent_exp,2000)/my_integral(UDetectorFlow,2000));
		//	UDetectorFlow->Scale(weight);
			
			
			TH1D* UDelta = delta_hist;
			UDelta->Divide(noflow_hist);
			UDelta->Multiply(mixed_exp);
			
			sameevent_exp->Add(mixed_exp, -1);
			
			sameevent_exp->Add(UDetectorFlow, 1);			
			

				
			UDelta->Scale(weight*my_integral(sameevent_exp,2000)/my_integral(UDelta,2000));
			sameevent_exp->Add(UDelta);	
			
			

			
			
			mixed_exp->Scale(sameevent_exp->GetEntries()/mixed_exp->GetEntries(), "nosw2");
			
			UDetectorFlow = flow_hist;
			UDetectorFlow->Divide(noflow_hist);
			UDetectorFlow->Multiply(mixed_exp);
		//	cout<<0.001*sameevent_exp->Integral()/UDetectorFlow->Integral()<<endl;
		//	UDetectorFlow->Scale(0.001*sameevent_exp->Integral()/UDetectorFlow->Integral());
			UDetectorFlow->Scale(0.001*my_integral(sameevent_exp,2000)/my_integral(UDetectorFlow,2000));
		//	UDetectorFlow->Scale(weight);
			
			sameevent_exp->Add(UDetectorFlow, 1);	
			

			
			
			TH1D *results = sameevent_exp;
			double diff=1;
			double lim=0.001;
			double previous=1;
			double current=2;
			int iterations=0;
			
			
		//	for(int i=0;i<2;i++)
			while(diff>lim)
			{
			mixed_exp->Scale(sameevent_exp->GetEntries()/mixed_exp->GetEntries(), "nosw2");
				UDetectorFlow = flow_hist;
				UDetectorFlow->Divide(noflow_hist);
				UDetectorFlow->Multiply(mixed_exp);
			//	cout<<0.001*sameevent_exp->Integral()/UDetectorFlow->Integral()<<endl;
			//	UDetectorFlow->Scale(0.001*sameevent_exp->Integral()/UDetectorFlow->Integral());
				UDetectorFlow->Scale(0.001*my_integral(sameevent_exp,2000)/my_integral(UDetectorFlow,2000), "nosw2");
				
				sameevent_exp->Add(UDetectorFlow, 1);	
				
				
				UDelta = delta_hist;
				UDelta->Divide(noflow_hist);
				UDelta->Multiply(mixed_exp);
				

				cout<<"UDelta scale 1st "<<my_integral(sameevent_exp,2000)/my_integral(UDelta,2000)<<endl;
				
				UDelta->Scale(weight*my_integral(sameevent_exp,2000)/my_integral(UDelta,2000), "nosw2");
				
				sameevent_exp->Add(UDelta, 1);
			
			
				cout<<"UDelta scale 2nd "<<my_integral(sameevent_exp,2000)/my_integral(UDelta,2000)<<endl;
				
				current=weight*my_integral(sameevent_exp,2000)/my_integral(UDelta,2000);
				UDelta->Scale(weight*my_integral(sameevent_exp,2000)/my_integral(UDelta,2000), "nosw2");
				
				diff = TMath::Abs(current-previous)/previous;
				cout<<"diff "<<diff<<endl;
				previous=current;
				iterations++;
			}
			cout<<endl<<"iterations: "<<iterations<<endl;
			
			UDelta->Scale(my_integral(sameevent_exp,2000)/my_integral(UDelta,2000));
			
			
			
			
			
			sprintf(histName, "delta %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			sameevent_exp->SetNameTitle(histName,histName);
			
			
			
		//	UDelta->Scale(my_integral(sameevent_exp,2000)/my_integral(UDelta,2000));
			
			sameevent_exp->SetStats(0);
			copy->SetStats(0);
			canvas.Clear();
			//canvas.SetLogy();
			UDelta->GetXaxis()->SetRangeUser(1000,2000);
			sameevent_exp->GetXaxis()->SetRangeUser(1000,2000);
		//	sameevent_exp->GetYaxis()->SetRangeUser(1.1*mixed_exp->GetMinimum(),1.1*sameevent_exp->GetMaximum());
			//sameevent_exp->Add(UDelta, 1);
			//sameevent_exp->Draw("hist");
			mixed_exp->SetLineColor(kCyan+3);
			
			
			sameevent_exp->Scale(1/(sameevent_exp->GetBinWidth(1)));
			
		//	cout<<my_integral(sameevent_exp, 2000)/	(UDelta, 2000)<<endl;
			
			UDetectorFlow->SetLineColor(kGreen+2);
			UDelta->SetLineColor(kRed);
			
		//	sameevent_exp->Sumw2(false);
			
			
			TLegend *legend = new TLegend(0.6,0.9,0.9,0.8);
			legend->AddEntry(sameevent_exp, "same-mixed+U(det+flow)+U(delta)", ",l");
			legend->AddEntry(UDetectorFlow, "U(detector,flow)", ",l");
			legend->AddEntry(UDelta, "U(#Delta)", ",l");
		//	legend->AddEntry(copy, "mixed bckg exp", ",l");
		//	legend->AddEntry(mixed_exp, "mixed rp", ",l");
		//	legend->AddEntry(multiplied_weight, "mixed rp*((mixed rp/mixed no rp)-1)", ",l");
		//	result->SetLineColor(kMagenta);
		//	result->SetLineColor(kRed);
		//	sameevent_exp->Draw("hist");
			sameevent_exp->Draw("");
		//	UDetectorFlow->Draw("same");
		//	UDelta->Draw("same");
		//	multiplied_weight->Draw("");
	//		mixed_exp->Draw("histsame");
		//	copy->Draw("same");
	//	mixed_weight->Draw();	
	//	mixed_exp->Draw("hist");	
	//		legend->Draw();
			sprintf(fileName, "./outputs/iterative.pdf");
			canvas.Print(fileName);  
			cout<<endl<<endl<<endl;
			
			
			output->cd();
			sameevent_exp->Write();
			
		}
	}
//	pimp_file->Close();
//	HALflow->Close();
//	HALnoflow->Close();
//	HALdelta->Close();
//	MC_file->Close();

	sprintf(fileName, "./outputs/iterative.pdf]");
	canvas.Print(fileName);  
//	}
	output->Save();
	
	return 0;
}