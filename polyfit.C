#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"


using namespace std;


int polyfit()
{
	
	TFile *file_fits = new TFile("./outputs/poly4_fits_HAL.root","recreate");
	
	
	TFile *nodelta_file = new TFile("./files/file_nodelta_flow_010.root", "read");
	TFile *nodelta0f_file = new TFile("./files/file_nodelta_noflow_010.root", "read");
	
	TF1* fits[8][8];

	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};


	TCanvas canvas("canvas");
	canvas.Print("./outputs/fits HAL.pdf["); 

	for(int ynum=0; ynum<8; ynum++)
	{
		for(int ptnum=0; ptnum<8; ptnum++)
		{
				
				int lower=1000;
				int upper=1800;
				
				
				
				
				char fitname[100];
				sprintf(fitname, "fit HAL %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				//fits[ynum][ptnum] = new TF1(fitname, "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4", lower, upper);
				fits[ynum][ptnum] = new TF1(fitname, "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4", lower, upper);
				
			
				char histName[100];
				sprintf(histName, "M_{inv} %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				cout<<histName<<endl;
				TH1D* hist_nodelta = (TH1D*)nodelta_file->Get(histName);
				TH1D* hist_nodelta0f = (TH1D*)nodelta0f_file->Get(histName);
			
				hist_nodelta->Rebin(10);
				hist_nodelta0f->Rebin(10);
			
				hist_nodelta->Scale(1/hist_nodelta->Integral());
				hist_nodelta0f->Scale(1/hist_nodelta0f->Integral());
				
			
			
				sprintf(histName, "M_{inv} divided %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				
				TH1D* hist_divided = new TH1D(*hist_nodelta);
				hist_divided->Divide(hist_nodelta0f);
				
			//	double parmin, parmax;
				
			//	fits[ynum][ptnum]->GetParLimits(0,parmin,parmax);
				
			//	cout<<endl<<endl<<endl<<parmin<<" "<<parmax<<endl<<endl;
				
				
				
		//		fits[ynum][ptnum]->SetParameters(1,0,0,0,0);
		//		fits[ynum][ptnum]->SetParError(0, 100);
		//		fits[ynum][ptnum]->SetParError(1, 100);
		//		fits[ynum][ptnum]->SetParError(2, 100);
		//		fits[ynum][ptnum]->SetParError(3, 100);
		//		fits[ynum][ptnum]->SetParError(4, 100);
			
			//	fits[ynum][ptnum]->SetParLimits(0,-100,200);
			//	fits[ynum][ptnum]->SetParLimits(1,0-100,200);
			//	fits[ynum][ptnum]->SetParLimits(2,-100,200);
			//	fits[ynum][ptnum]->SetParLimits(3,-100,1);
			//	fits[ynum][ptnum]->SetParLimits(4,-200,200);
			
			
				ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000); 
				hist_divided->Fit(fitname, "M", "", lower, upper);
				file_fits->cd();
				fits[ynum][ptnum]->Write();
				cout<<endl<<endl<<endl;
				cout<<"nPar "<<fits[ynum][ptnum]->GetNpar()<<endl;
				cout<<"NDF "<<fits[ynum][ptnum]->GetNDF()<<endl;
				cout<<"chi^2 "<<fits[ynum][ptnum]->GetChisquare()<<endl;
				cout<<"prob "<<fits[ynum][ptnum]->GetProb()<<endl;
				cout<<endl<<endl<<endl;
				
				canvas.Clear();
				hist_divided->GetYaxis()->SetRangeUser(0, 2);
			//	hist_divided->SetMarkerSize(2);
				hist_divided->Draw();
				fits[ynum][ptnum]->Draw("same");
				
	
				
				canvas.Print("./outputs/fits HAL.pdf");    
				
		}
	}
	
	
	canvas.Print("./outputs/fits HAL.pdf]");    

/*


	TFile *file_fits_exp = new TFile("./outputs/poly4_fits_exp.root","recreate");
	
	
	TFile *sameevent_file = new TFile("sameevent_pim_010_100k.root", "read");
	TFile *background_file = new TFile("background_pim_010_10k.root", "read");
	
	TF1* fits_exp[8][8];

	TCanvas canvas1("canvas");
	canvas.Print("./outputs/fits exp.pdf["); 



	for(int ynum=0; ynum<8; ynum++)
	{
		for(int ptnum=0; ptnum<8; ptnum++)
		{
				char fitname1[100];
				sprintf(fitname1, "fit exp %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));

			
				char histName1[100];
				sprintf(histName1, "same-event %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				cout<<histName1<<endl;
				TH1D* hist_sameevent = (TH1D*)sameevent_file->Get(histName1);
				
				sprintf(histName1, "bckg %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				
				TH1D* hist_bckg = (TH1D*)background_file->Get(histName1);
			
				
				
				hist_sameevent->Sumw2();
				hist_bckg->Sumw2();
				
			
				hist_sameevent->Rebin(10);
				hist_bckg->Rebin(10);
			
				hist_sameevent->Scale(1/hist_sameevent->Integral());
				hist_bckg->Scale(1/hist_bckg->Integral());
				
				
				TH1D* hist_divided = new TH1D(*hist_sameevent);
				hist_divided->Divide(hist_bckg);
				
				
				
				auto f = new TF1("f", "([0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4)*([5]+[6]*x+[7]*x^2)", 1070, 1800);
				//auto f = new TF1("f", "([0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4)", 1070, 1800);
				
				f->SetParameter(0, fits[ynum][ptnum]->GetParameter(0));
				f->SetParameter(1, fits[ynum][ptnum]->GetParameter(1));
				f->SetParameter(2, fits[ynum][ptnum]->GetParameter(2));
				f->SetParameter(3, fits[ynum][ptnum]->GetParameter(3));
				f->SetParameter(4, fits[ynum][ptnum]->GetParameter(4));
				
				
				ROOT::Fit::DataRange range;
				range.AddRange(0, 1000, 1100);
				range.AddRange(0, 1400, 1800);
				auto fitFunc = [&](double *x, double *p) {
					if (!range.IsInside(x[0])) {
						TF1::RejectPoint();
						return 0.0;
					}
					return f->EvalPar(x,p);
				};
				auto f1 = new TF1(fitname1,fitFunc,1070,1800,f->GetNpar());
				
				
			//	f1->FixParameter(0, fits[ynum][ptnum]->GetParameter(0));
			//	f1->FixParameter(1, fits[ynum][ptnum]->GetParameter(1));
			//	f1->FixParameter(2, fits[ynum][ptnum]->GetParameter(2));
			//	f1->FixParameter(3, fits[ynum][ptnum]->GetParameter(3));
			//	f1->FixParameter(4, fits[ynum][ptnum]->GetParameter(4));
			
				f1->FixParameter(0, fits[ynum][ptnum]->GetParameter(0));
				f1->FixParameter(1, fits[ynum][ptnum]->GetParameter(1));
				f1->FixParameter(2, fits[ynum][ptnum]->GetParameter(2));
				f1->FixParameter(3, fits[ynum][ptnum]->GetParameter(3));
				f1->FixParameter(4, fits[ynum][ptnum]->GetParameter(4));
				
				hist_divided->Fit(f1, "N");
				
				
				TF1 *tmp = new TF1("tmp", "([0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4)*([5]+[6]*x+[7]*x^2)", 1070, 2400);
				
				tmp->SetParameter(0, f1->GetParameter(0));
				tmp->SetParameter(1, f1->GetParameter(1));
				tmp->SetParameter(2, f1->GetParameter(2));
				tmp->SetParameter(3, f1->GetParameter(3));
				tmp->SetParameter(4, f1->GetParameter(4));
				tmp->SetParameter(5, f1->GetParameter(5));
				tmp->SetParameter(6, f1->GetParameter(6));
				tmp->SetParameter(7, f1->GetParameter(7));
				
				
				TF1 *multiplied = new TF1("multiplied", "[5]+[6]*x+[7]*x^2", 1070, 1800);
				multiplied->SetParameter(5, f1->GetParameter(5));
				multiplied->SetParameter(6, f1->GetParameter(6));
				multiplied->SetParameter(7, f1->GetParameter(7));
				
	
				file_fits->cd();
			//	fits_exp[ynum][ptnum]->Write();
				file_fits_exp->cd();
				f1->Write();
				canvas1.Clear();
				hist_divided->GetYaxis()->SetRangeUser(0, 2);
				hist_divided->SetMarkerSize(2);
				hist_divided->SetStats(0);
				hist_divided->Draw();
				tmp->SetLineColor(kGreen);
				tmp->Draw("SAME");
				f1->Draw("SAME");
				multiplied->SetLineColor(kMagenta);
				multiplied->Draw("SAME");
				
				TLegend *legend = new TLegend(0.7,0.9,0.9,0.8);
			
				legend->AddEntry(f1, "fit", ",l");
				legend->AddEntry(tmp, "whole function", "l");
				legend->AddEntry(multiplied, "multiplication factor", "l");
				
				legend->Draw();
				
				
				cout<<endl<<endl<<endl;
				cout<<"nPar "<<f1->GetNpar()<<endl;
				cout<<"NDF "<<f1->GetNDF()<<endl;
				cout<<"chi^2 "<<f1->GetChisquare()<<endl;
				cout<<"prob "<<f1->GetProb()<<endl;
				cout<<endl<<endl<<endl;
				
				
				
				
				//fits[ynum][ptnum]->SetLineColor(kGreen);
				//fits[ynum][ptnum]->Draw("same");
				
				canvas1.Print("./outputs/fits exp.pdf");    
				
		}
	}
	
	
	canvas1.Print("./outputs/fits exp.pdf]");    
*/



	
	return 0;
}