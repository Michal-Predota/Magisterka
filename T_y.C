


int T_y()
{
	
	TFile *file1 = new TFile("./files/temp_010_pip.root","read");
	TFile *file2 = new TFile("./files/temp_1020_pip.root","read");
	TFile *file3 = new TFile("./files/temp_2030_pip.root","read");
	TFile *file4 = new TFile("./files/temp_3040_pip.root","read");
	
	
	TH1D *temp1 = (TH1D*)file1->Get("#Delta^{++}");
	TH1D *temp2 = (TH1D*)file2->Get("#Delta^{++}");
	TH1D *temp3 = (TH1D*)file3->Get("#Delta^{++}");
	TH1D *temp4 = (TH1D*)file4->Get("#Delta^{++}");
	
	
	
	TCanvas *c = new TCanvas("c","c",800,800);
	c->cd();
	gPad->SetBottomMargin(0.12);
	gPad->SetTopMargin(0.08);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.12);
	
	temp1->SetStats(0);
	temp1->SetMarkerStyle(21);
	temp1->SetMarkerSize(2);
	temp1->SetMarkerColor(temp1->GetLineColor());
	
	temp1->GetXaxis()->SetTitle("y");
	temp1->GetYaxis()->SetTitle("T_{eff} (MeV)");
	temp1->SetTitle("#Delta^{++} effective temperature");
	temp1->GetYaxis()->SetTitleSize(0.05);
	temp1->GetXaxis()->SetTitleSize(0.05);
	temp1->GetYaxis()->SetLabelSize(0.04);
	temp1->GetXaxis()->SetLabelSize(0.04);
	
	temp1->GetYaxis()->SetRangeUser(80,200);
	
	
	temp1->Draw();
	
	
	temp2->SetLineColor(kRed+2);
	temp2->SetMarkerStyle(21);
	temp2->SetMarkerSize(2);
	temp2->SetMarkerColor(temp2->GetLineColor());
	temp2->Draw("same");
	
	
	temp3->SetLineColor(kGreen+2);
	temp3->SetMarkerStyle(21);
	temp3->SetMarkerSize(2);
	temp3->SetMarkerColor(temp3->GetLineColor());
	temp3->Draw("same");
	
	
	temp4->SetLineColor(kCyan-3);
	temp4->SetMarkerStyle(21);
	temp4->SetMarkerSize(2);
	temp4->SetMarkerColor(temp4->GetLineColor());
	temp4->Draw("same");
	
	TLegend* legend = new TLegend(0.75, 0.6, 0.98, 0.92);
	
	legend->AddEntry(temp1, "0-10% centrality", "lp");
	legend->AddEntry(temp2, "10-20% centrality", "lp");
	legend->AddEntry(temp3, "20-30% centrality", "lp");
	legend->AddEntry(temp4, "30-40% centrality", "lp");
	
	legend->Draw();
	
	
	c->SaveAs("./images/pip/temperature/temperature.png");
	
	
	return 0;
}