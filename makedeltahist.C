#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"


int makedeltahist()
{
	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	
	
	TFile *file = new TFile("./files/delta++_3040.root","read");
	TFile *file1 = new TFile("./files/fits_d++_3040.root","read");
	
	
	for(int ynum=0; ynum<8; ynum++)
	{
	//	int ynum=0;
	//	int ptnum=0;
	
		char cName[100];
		sprintf(cName, "./images/pip/c3040/delta++_1 %f<y<%f.png", ylim.at(ynum), ylim.at(ynum+1));
		TCanvas *c = new TCanvas("c","c",1600,1600);
		c->Divide(2,2);
		for(int ptnum=0; ptnum<4; ptnum++)
		{
				char histName[100];
				sprintf(histName, "delta %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				TH1D *delta = (TH1D*)file->Get(histName);
				delta->GetXaxis()->SetRangeUser(1000,1800);
				sprintf(histName, "delta %0.1f<y<%0.1f, %0.1f<pT<%0.1f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				delta->SetNameTitle(histName,histName);
				
				c->cd(ptnum+1);
				
				gPad->SetBottomMargin(0.15);
				gPad->SetTopMargin(0.15);
				gPad->SetLeftMargin(0.15);
				gPad->SetRightMargin(0.15);
				
				
				
				delta->SetMarkerStyle(8);
				delta->SetMarkerSize(1);
				delta->SetMarkerColor(delta->GetLineColor());
				
				delta->SetTitleSize(10,"t");
				delta->GetYaxis()->SetTitle("1/N dN/dM_{inv}");
				delta->GetXaxis()->SetTitle("M_{inv} (MeV/c^{2})");
				delta->GetXaxis()->SetTitleOffset(1.0);
				delta->GetYaxis()->SetTitleOffset(1.3);
				delta->GetXaxis()->SetTitleSize(0.06);
				delta->GetYaxis()->SetTitleSize(0.06);
				delta->GetXaxis()->SetLabelSize(0.06);
				delta->GetYaxis()->SetLabelSize(0.06);
						
				for(int i=0; i<delta->GetNbinsX(); i++)
				{
					if(delta->GetBinContent(i+1)<0)
						delta->SetBinContent(i+1,0);
				}
				
				
				char fitName[100];
				sprintf(fitName, "breit_wigner %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				TF1 *fit = (TF1*)file1->Get(fitName);
				
				fit->SetLineWidth(2);
				
				
				
				delta->Draw("");
				fit->Draw("same");
		}
		
		
		

		char cName1[100];
		sprintf(cName1, "./images/pip/c3040/delta++_2 %f<y<%f.png", ylim.at(ynum), ylim.at(ynum+1));
		TCanvas *c1 = new TCanvas("c1","c1",1600,1600);
		c1->Divide(2,2);
		
		for(int ptnum=4; ptnum<8; ptnum++)
		{
				char histName[100];
				sprintf(histName, "delta %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				TH1D *delta = (TH1D*)file->Get(histName);
				delta->GetXaxis()->SetRangeUser(1000,1800);
				sprintf(histName, "#Delta^{++} %0.1f<y<%0.1f, %0.1f<pT<%0.1f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				delta->SetNameTitle(histName,histName);
				
				
				c1->cd(ptnum-3);
				
				gPad->SetBottomMargin(0.15);
				gPad->SetTopMargin(0.15);
				gPad->SetLeftMargin(0.15);
				gPad->SetRightMargin(0.15);
				
				
				
				delta->SetMarkerStyle(8);
				delta->SetMarkerSize(1);
				delta->SetMarkerColor(delta->GetLineColor());
				
				delta->SetTitleSize(10,"t");
				delta->GetYaxis()->SetTitle("1/N dN/dM_{inv}");
				delta->GetXaxis()->SetTitle("M_{inv} (MeV/c^{2})");
				delta->GetXaxis()->SetTitleOffset(1.0);
				delta->GetYaxis()->SetTitleOffset(1.3);
				delta->GetXaxis()->SetTitleSize(0.06);
				delta->GetYaxis()->SetTitleSize(0.06);
				delta->GetXaxis()->SetLabelSize(0.06);
				delta->GetYaxis()->SetLabelSize(0.06);
				
						
				for(int i=0; i<delta->GetNbinsX(); i++)
				{
					if(delta->GetBinContent(i+1)<0)
						delta->SetBinContent(i+1,0);
				}
				
				
				char fitName[100];
				sprintf(fitName, "breit_wigner %f<y<%f, %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
				TF1 *fit = (TF1*)file1->Get(fitName);
				
				fit->SetLineWidth(2);
				
				
				delta->Draw("");
				fit->Draw("same");
		}
		
		c->SaveAs(cName);
		c1->SaveAs(cName1);
	}
	
	
	TH2D* masses = (TH2D*)file1->Get("masses");
	TH2D* widths = (TH2D*)file1->Get("widths");
	TH2D* masses_errors = (TH2D*)file1->Get("masses_errors");
	TH2D* widths_errors = (TH2D*)file1->Get("widths_errors");
	TH2D* chi2ndf = (TH2D*)file1->Get("chi2ndf");
	
	
	vector<TLine*> lines;
	vector<double> center = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	
	for(int i=0; i<center.size(); i++)
	{
		TLine *l = new TLine(center.at(i),0,center.at(i),1.6);
		lines.push_back(l);
	}
	
	for(int i=0; i<center.size(); i++)
	{
		TLine *l = new TLine(0,center.at(i),1.6,center.at(i));
		lines.push_back(l);
	}
	

	
	
	TCanvas *c = new TCanvas("c","c",800,800);
	
	c->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.08);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.2);
	
	masses->GetYaxis()->SetTitle("p_{T} (MeV/c)");
	masses->GetXaxis()->SetTitle("y");
	masses->GetYaxis()->SetTitleSize(0.042);
	masses->GetXaxis()->SetTitleSize(0.042);
	
	masses->GetYaxis()->SetLabelSize(0.04);
	masses->GetXaxis()->SetLabelSize(0.04);
	
	masses->GetZaxis()->SetRangeUser(1000, masses->GetMaximum());
	TPaletteAxis *palette = (TPaletteAxis*)masses->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.85);
    palette->SetX2NDC(0.9);
    palette->SetY1NDC(0.15);
    palette->SetY2NDC(0.9);
    gPad->Modified();
    gPad->Update();

	
	masses->GetXaxis()->SetTitleOffset(1.6);
	masses->GetYaxis()->SetTitleOffset(1.6);
	
	masses->SetNameTitle("#Delta^{++} masses","#Delta^{++} masses");
	
	masses->Draw("colzTEXT");
	
	
	for(int i=0; i<lines.size(); i++)
	{
		lines.at(i)->Draw("same");
	}
	
	c->SaveAs("./images/pip/c3040/delta++_masses_c3040.png");
	
	
	
	
	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	
	c1->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.08);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.2);
	
	widths->GetYaxis()->SetTitle("p_{T} (MeV/c)");
	widths->GetXaxis()->SetTitle("y");
	widths->GetYaxis()->SetTitleSize(0.042);
	widths->GetXaxis()->SetTitleSize(0.042);
	
	widths->GetYaxis()->SetLabelSize(0.04);
	widths->GetXaxis()->SetLabelSize(0.04);
	
	widths->GetZaxis()->SetRangeUser(0, widths->GetMaximum());
	palette = (TPaletteAxis*)widths->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.85);
    palette->SetX2NDC(0.9);
    palette->SetY1NDC(0.15);
    palette->SetY2NDC(0.9);
    gPad->Modified();
    gPad->Update();

	
	widths->GetXaxis()->SetTitleOffset(1.6);
	widths->GetYaxis()->SetTitleOffset(1.6);
	
	widths->SetNameTitle("#Delta^{++} widths","#Delta^{++} widths");
	
	widths->Draw("colzTEXT");
	
	
	for(int i=0; i<lines.size(); i++)
	{
		lines.at(i)->Draw("same");
	}
	
	c1->SaveAs("./images/pip/c3040/delta++_widths_c3040.png");
	
	
	
	
	TCanvas *c2 = new TCanvas("c2","c2",800,800);
	
	c2->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.08);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.2);
	
	masses_errors->GetYaxis()->SetTitle("p_{T} (MeV/c)");
	masses_errors->GetXaxis()->SetTitle("y");
	masses_errors->GetYaxis()->SetTitleSize(0.042);
	masses_errors->GetXaxis()->SetTitleSize(0.042);
	
	masses_errors->GetYaxis()->SetLabelSize(0.04);
	masses_errors->GetXaxis()->SetLabelSize(0.04);
	
	masses_errors->GetZaxis()->SetRangeUser(0, masses_errors->GetMaximum());
	palette = (TPaletteAxis*)masses_errors->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.85);
    palette->SetX2NDC(0.9);
    palette->SetY1NDC(0.15);
    palette->SetY2NDC(0.9);
    gPad->Modified();
    gPad->Update();

	
	masses_errors->GetXaxis()->SetTitleOffset(1.6);
	masses_errors->GetYaxis()->SetTitleOffset(1.6);
	
	masses_errors->SetNameTitle("#Delta^{++} masses estimated errors","#Delta^{++} masses estimated errors");
	
	masses_errors->Draw("colzTEXT");
	
	
	for(int i=0; i<lines.size(); i++)
	{
		lines.at(i)->Draw("same");
	}
	
	c2->SaveAs("./images/pip/c3040/delta++_masses_errors_c3040.png");
	
	
	
	
	TCanvas *c3 = new TCanvas("c3","c3",800,800);
	
	c3->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.08);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.2);
	
	widths_errors->GetYaxis()->SetTitle("p_{T} (MeV/c)");
	widths_errors->GetXaxis()->SetTitle("y");
	widths_errors->GetYaxis()->SetTitleSize(0.042);
	widths_errors->GetXaxis()->SetTitleSize(0.042);
	
	widths_errors->GetYaxis()->SetLabelSize(0.04);
	widths_errors->GetXaxis()->SetLabelSize(0.04);
	
	widths_errors->GetZaxis()->SetRangeUser(0, widths_errors->GetMaximum());
	palette = (TPaletteAxis*)widths_errors->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.85);
    palette->SetX2NDC(0.9);
    palette->SetY1NDC(0.15);
    palette->SetY2NDC(0.9);
    gPad->Modified();
    gPad->Update();

	
	widths_errors->GetXaxis()->SetTitleOffset(1.6);
	widths_errors->GetYaxis()->SetTitleOffset(1.6);
	
	widths_errors->SetNameTitle("#Delta^{++} widths estimated errors","#Delta^{++} widths estimated errors");
	
	widths_errors->Draw("colzTEXT");
	
	
	for(int i=0; i<lines.size(); i++)
	{
		lines.at(i)->Draw("same");
	}
	
	c3->SaveAs("./images/pip/c3040/delta++_widths_errors_c3040.png");
	
	
	
	
	TCanvas *c4 = new TCanvas("c4","c4",800,800);
	
	c4->cd();
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.08);
	gPad->SetLeftMargin(0.18);
	gPad->SetRightMargin(0.2);
	
	chi2ndf->GetYaxis()->SetTitle("p_{T} (MeV/c)");
	chi2ndf->GetXaxis()->SetTitle("y");
	chi2ndf->GetYaxis()->SetTitleSize(0.05);
	chi2ndf->GetXaxis()->SetTitleSize(0.05);
	
	chi2ndf->GetYaxis()->SetLabelSize(0.05);
	chi2ndf->GetXaxis()->SetLabelSize(0.05);
	
	chi2ndf->GetZaxis()->SetRangeUser(0, chi2ndf->GetMaximum());
	palette = (TPaletteAxis*)chi2ndf->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.85);
    palette->SetX2NDC(0.9);
    palette->SetY1NDC(0.15);
    palette->SetY2NDC(0.9);
    gPad->Modified();
    gPad->Update();

	
	chi2ndf->GetXaxis()->SetTitleOffset(1.6);
	chi2ndf->GetYaxis()->SetTitleOffset(1.6);
	
	chi2ndf->SetNameTitle("#Delta^{++} #chi^{2}/ndf of fits","#Delta^{++} #chi^{2}/ndf of fits");
	
	chi2ndf->Draw("colzTEXT");
	
	
	for(int i=0; i<lines.size(); i++)
	{
		lines.at(i)->Draw("same");
	}
	
	c4->SaveAs("./images/pip/c3040/delta++_chi2ndf_c3040.png");
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	return 0;
}