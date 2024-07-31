#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

using namespace std;


int flowhist()
{
	//TFile *all_file = new TFile("./files/file_all.root", "read");
	TFile *nodelta_file = new TFile("./files/file_nodelta_flow_1020.root", "read");
	TFile *output = new TFile("./outputs/flow_pim_1020.root", "recreate");
	TFile *nodelta0f_file = new TFile("./files/file_nodelta_noflow_1020.root", "read");
	
	//TFile *MC_bckg = new TFile("./background_pim_sim_010.root", "read");
	//TFile *MC_same = new TFile("./sameevent_pim_sim_010.root", "read");
	
	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	
	TH1D *plots[8][8];
	TCanvas* canvases[8];
	TCanvas canvas("canvas");
	//canvas.Print("nodelta divided by all HAL.pdf["); 
	canvas.Print("./outputs/flow_pim_1020.pdf["); 
	int filenum=0;
	for(int ynum=0; ynum<8; ynum++)
	{
		for(int ptnum=0; ptnum<8; ptnum++)
		{
			char histName[100];
			sprintf(histName, "M_{inv} %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			cout<<histName<<endl;
			//TH1D* hist_all = (TH1D*)all_file->Get(histName);
			TH1D* num = (TH1D*)nodelta_file->Get(histName);
			sprintf(histName, "M_{inv} %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* den = (TH1D*)nodelta0f_file->Get(histName);
			
			
			
			//hist_all->Rebin(10);
			num->Rebin(10);
			den->Rebin(10);
			
			num->Scale(1/num->Integral());
			den->Scale(1/den->Integral());
			
			TH1D* hist_divided = new TH1D(*num);
			//hist_divided->Divide(hist_all);
			hist_divided->Divide(den);
			
			
			sprintf(histName, "divided %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			hist_divided->SetNameTitle(histName,  histName);
			//output->cd();
		//	hist_divided->GetYaxis()->SetRangeUser(0, 2);
			
			//hist_divided->Write();
			
			canvas.Clear();
			cout<<ynum<<" "<<ptnum<<endl;
			hist_divided->Sumw2(false);
			hist_divided->GetYaxis()->SetRangeUser(0.8, 1.2);
			hist_divided->Draw();
			//num->Draw();
			//den->Draw();
			//canvas.Print("nodelta divided by all HAL.pdf");    
			canvas.Print("./outputs/flow_pim_1020.pdf");    
		}
	}
	
	/*
	TFile *file1= new TFile("./output_nodelta/superpack_0/th1_8/histo1d.root", "read");
	TCanvas *c1 = (TCanvas*) file1->Get("canvas");
	TH1D* hist1 = (TH1D*)c1->GetPrimitive("num");
	
	
	TFile *file2= new TFile("./output_all/superpack_0/th1_8/histo1d.root", "read");
	TCanvas *c2 = (TCanvas*) file2->Get("canvas");
	TH1D* hist2 = (TH1D*)c2->GetPrimitive("num");
	
	hist1->Rebin(8);
	hist2->Rebin(8);
	
	hist1->Scale(1/hist1->Integral(), "nosw2");
	hist2->Scale(1/hist2->Integral(), "nosw2");
	
	hist1->Divide(hist2);
	hist1->GetYaxis()->SetRangeUser(0.5, 1.5);
	hist1->SetNameTitle("divided all","divided all");
	output->cd();
	hist1->Write();
	canvas.cd();
	canvas.Clear();
	hist1->Draw("hist");
	canvas.Print("nodelta divided by all HAL HAL.pdf");  
	//canvas.Clear();
	//hist2->Draw("hist");
	//canvas.Print("divided HAL.pdf");  
	
	canvas.Print("nodelta divided by all HAL.pdf]");
	*/
	
	
	
	canvas.Print("./outputs/flow_pim_1020.pdf]");
	
	//output->Save();
	//
	//all_file->Close();
	//nodelta_file->Close();
	//output->Close();
	
	
	
	
	
	
	/*
	
	TLegend *legends[8];
	vector<int> colors = {kBlue+2, kRed+1, kMagenta+3, kOrange+10, kGreen+3, kCyan+3, kYellow+3, kCyan-6};
	
	
	TFile *divided = new TFile("divided_inverted.root","read");
	TFile *constypt = new TFile("constant_ypt.root", "recreate");
	
	TCanvas canvasY("canvasY");
	canvasY.Print("constant y.pdf["); 
	char histName[100];
	
	
	for(int ynum=0; ynum<8; ynum++)
	{
		char canName [100];
		sprintf(canName, "%f<y<%f", ylim.at(ynum), ylim.at(ynum+1));
		canvases[ynum] = new TCanvas(canName, canName,800,800);
		legends[ynum] = new TLegend(0.5,0.95,0.9,0.75);
		for(int ptnum=0; ptnum<8; ptnum++)
		{
			sprintf(histName, "divided %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			TH1D* hist_divided = (TH1D*)divided->Get(histName);
			canvases[ynum]->cd();
			//hist_divided->Draw();
			plots[ynum][ptnum] = new TH1D(*hist_divided);
			canvases[ynum]->cd();
			plots[ynum][ptnum]->SetStats(0);
			//plots[ynum][ptnum]->GetYaxis()->SetRangeUser(0.95, 1.05);
			plots[ynum][ptnum]->GetXaxis()->SetRangeUser(1000, 1600);
			plots[ynum][ptnum]->SetLineColor(colors.at(ptnum));
			plots[ynum][ptnum]->SetLineWidth(2);
			plots[ynum][ptnum]->Draw("same");
			
			legends[ynum]->AddEntry(plots[ynum][ptnum], histName, "l");
		}
		
		sprintf(histName, "%f<y<%f, varying pT", ylim.at(ynum), ylim.at(ynum+1));
		plots[ynum][0]->SetNameTitle(histName,histName);
		legends[ynum]->Draw();
		constypt->cd();
		canvases[ynum]->Write();
		
		
		
			canvasY.Clear();
			canvases[ynum]->DrawClonePad();
			canvasY.Print("constant y.pdf");    
		
	}
	
	canvasY.Print("constant y.pdf]");
	
	
	TH1D *plotspt[8][8];
	TCanvas* canvasespt[8];
	
	
	TCanvas canvasPt("canvasPt");
	canvasPt.Print("constant Pt.pdf["); 
	
	filenum=0;
	
	
	for(int ptnum1=0; ptnum1<8; ptnum1++)
	{
		char canName [100];
		sprintf(canName, "%f<pT<%f", pTlim.at(ptnum1), pTlim.at(ptnum1+1));
		canvasespt[ptnum1] = new TCanvas(canName, canName,800,800);
		legends[ptnum1] = new TLegend(0.5,0.95,0.9,0.75);
		for(int ynum1=0; ynum1<8; ynum1++)
		{
			sprintf(histName, "divided %f<y<%f,  %f<pT<%f", ylim.at(ynum1), ylim.at(ynum1+1), pTlim.at(ptnum1), pTlim.at(ptnum1+1));
			TH1D* hist_divided = (TH1D*)divided->Get(histName);
			canvasespt[ptnum1]->cd();
			//hist_divided->Draw();
			plotspt[ptnum1][ynum1] = new TH1D(*hist_divided);
			canvasespt[ptnum1]->cd();
			plotspt[ptnum1][ynum1]->SetStats(0);
			//plotspt[ptnum1][ynum1]->GetYaxis()->SetRangeUser(0.95, 1.05);
			plotspt[ptnum1][ynum1]->GetXaxis()->SetRangeUser(1000, 1600);
			plotspt[ptnum1][ynum1]->SetLineColor(colors.at(ynum1));
			plotspt[ptnum1][ynum1]->SetLineWidth(2);
			plotspt[ptnum1][ynum1]->Draw("same");
			
			legends[ptnum1]->AddEntry(plotspt[ptnum1][ynum1], histName, "l");
		}
		
		sprintf(histName, "%f<pT<%f, varying y", pTlim.at(ptnum1), pTlim.at(ptnum1+1));
		plotspt[ptnum1][0]->SetNameTitle(histName,histName);
		legends[ptnum1]->Draw();
		constypt->cd();
		canvasespt[ptnum1]->Write();
		
		
			canvasPt.Clear();
			canvasespt[ptnum1]->DrawClonePad();
			canvasPt.Print("constant Pt.pdf");    
	}
	
	
		canvasPt.Print("constant Pt.pdf]");    
	
	
	
	
	constypt->Save();
	
	*/
	
	return 0;
}
