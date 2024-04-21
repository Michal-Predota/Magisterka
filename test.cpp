#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>

#include "/home/ytcvc/HAL/analysis/spectra/Decay.h"



int test()
{
	
	TF1* breit_wigner = new TF1("breit_wigner", "[0]*(x*x*0.117)/(TMath::Power((1.232*1.232-x*x), 2)+(x*x*0.117*0.117))", 0, 2);
	breit_wigner->SetParameter(0, 1);
	double norm = breit_wigner->Integral(0, 2);
	//breit_wigner->SetParameter(0, 1/norm);
	
	
	TH1D* h = new TH1D("h","h", 200, 0, 2000);
	
	
	/*for(int i=0; i<10000; i++)
	{
		h->Fill(breit_wigner->GetRandom()*1000);
	}
	
	//h->Draw();
	


	TH1D* dau = new TH1D("dau","dau", 500, 0, 2);
	
	
	Hal::DecayChannel ch1(2212, 211, 10);
	Hal::Decay decay(2224);
	decay.AddDecayChannel(ch1);
	decay.Init();
	std::vector<Hal::McTrack*> tracks;
	for (int i = 0; i < 3; i++) {
	tracks.push_back(new Hal::McTrack());
	}
	Hal::McTrack Delta;
		
	int n_deltas=0;
		
	while(n_deltas<100000)
	{
		Double_t px = gRandom->Gaus(0, 0.5);
		Double_t py = gRandom->Gaus(0, 0.5);
		Double_t pz = gRandom->Gaus(0, 0.5);
		TLorentzVector sum;
		double m = breit_wigner->GetRandom();
		sum.SetXYZM(px, py, pz, m);
		Delta.SetMomentum(px, py, pz, sum.E());
		decay.DecayParticle(Delta, tracks, false, m);
		
		
		dau->Fill((tracks[0]->GetMomentum()+tracks[1]->GetMomentum()).M());
		
		//std::cout<<"Delta mass: "<<m<<"	mass from tracks: "<<(tracks[0]->GetMomentum()+tracks[1]->GetMomentum()).M()<<std::endl;
		n_deltas++;
	}
	*/
	
	
	/*TCanvas* c = new TCanvas("c","c",800,800);
	c->cd();
	
	
	
	TFile* f1 = new TFile("e5000d2.root");
	TFile* f2 = new TFile("nodeltas_noflow.root");
	
	TH1D* deltas = (TH1D*)f1->Get("den");
	TH1D* nodeltas = (TH1D*)f2->Get("den");
	
	
	norm = deltas->Integral(deltas->FindBin(1000), deltas->FindBin(1150))/nodeltas->Integral(deltas->FindBin(1000), deltas->FindBin(1150));
	std::cout<<deltas->Integral(deltas->FindBin(1100), deltas->FindBin(1150))<<std::endl;
	
	nodeltas->Scale(norm);
	//deltas->Draw();
	//nodeltas->SetLineColor(kRed);
	//nodeltas->Draw("SAME");
	
	
	deltas->Add(nodeltas, -1);
	deltas->Draw();*/
	
	
	
	/*TCanvas* c = new TCanvas("c","c",800,800);
	c->cd();
	
	TFile* f1 = new TFile("results.root");
	TH1D* experimental = (TH1D*)f1->Get("invmass");
	experimental->Sumw2(false);
	//experimental->SetBins(1000, 1000, 2500);
	
	
	//minvhal.root - 10k events
	int nEvents = 10000;
	TFile* f2 = new TFile("bckg_1000bin.root");
	TH1D* hal = (TH1D*)f2->Get("num");
	hal->SetLineColor(kRed);
	
	
//	TFile* f3 = new TFile("analysis_pim_no_cut.root");
//	TH1D* experimental_nocut = (TH1D*)f3->Get("#pi^{-} p");
//	experimental_nocut->SetLineColor(kGreen);
	//experimental_nocut->SetBins(1000, 1000, 2500);
	
	
	
	experimental->Scale(1/(experimental->Integral()));
//	experimental_nocut->Scale(1/(experimental_nocut->Integral()), "nosw2");
	hal->Scale(1/(hal->Integral()));
	
	cout<<hal->GetNbinsX()<<"	"<<experimental->GetNbinsX()<<endl;
	
	cout<<hal->Integral()<<"	"<<experimental->Integral()<<endl;

	
	//experimental->GetYaxis()->SetTitle("1/N dN/dm");
	gPad->SetLeftMargin(0.15);
	experimental->SetStats(000);
	hal->SetStats(000);
	
	TLegend *legend = new TLegend(0.9,0.9,0.7,0.8);
	legend->AddEntry(hal,"num", "l");
	legend->AddEntry(experimental,"pim p", "l");
//	legend->AddEntry(experimental_nocut,"pim p no cuts", "l");
	


	//experimental_nocut->Draw();
	hal->Draw("SAME");
	experimental->Draw("SAME");
	
//	experimental_nocut->Add(hal, -1);
//	experimental_nocut->Sumw2(false);
//	experimental_nocut->SetLineColor(kBlue);
	
//	experimental_nocut->Draw();
	
	legend->Draw();*/
	
	
	
	
	
	TFile* f1 = new TFile("pty_p.root");
	TFile* f2 = new TFile("pty_pi.root");
	
	TFile* f3 = new TFile("analysis_alt.root");
	TFile* f4 = new TFile("spec_pim.root");
	
	TH2D* ptyHAL_p = (TH2D*)f1->Get("pty proton");
	TH2D* ptyHAL_pi = (TH2D*)f2->Get("pty pion");
	
	TH2D* ptyEXP_p = (TH2D*)f3->Get("pty_proton");
	TH2D* ptyEXP_pi = (TH2D*)f3->Get("pty_pim");
	
	ptyHAL_p->Scale(1/(10000*ptyHAL_p->GetXaxis()->GetBinWidth(1)*ptyHAL_p->GetYaxis()->GetBinWidth(1)));
	ptyHAL_pi->Scale(1/(10000*ptyHAL_pi->GetXaxis()->GetBinWidth(1)*ptyHAL_pi->GetYaxis()->GetBinWidth(1)));
	
//	ptyEXP_p->Scale(1/ptyEXP_p->Integral());
//	ptyEXP_pi->Scale(1/ptyEXP_pi->Integral());
	
//	cout<<ptyEXP_pi->GetXaxis()->GetLast()<<endl;
	
	//ptyEXP_pi->Draw("colz");
	
//	ptyEXP_pi->Divide(ptyHAL_pi);
//	ptyEXP_pi->Draw("COLZ");
	
	
	
	
	TCanvas *cp = new TCanvas("cp","cp",800,800);
	cp->Divide(1,2);
	cp->cd(1);
	ptyHAL_p->Draw("colz");
	cp->cd(2);
	ptyEXP_p->DrawCopy("colz");
	
	TCanvas *cpi = new TCanvas("cpi","cpi",800,800);
	cpi->Divide(1,2);
	cpi->cd(1);
	ptyHAL_pi->Draw("colz");
	cpi->cd(2);
	ptyEXP_pi->DrawCopy("colz");
	
	
	TCanvas *c = new TCanvas("c","c",800,800);
	c->Divide(1,2);
	c->cd(1);
	ptyEXP_pi->Divide(ptyHAL_pi);
	ptyEXP_pi->SetName("pi_ratio");
	ptyEXP_pi->SetTitle("exp/HAL pions");
	ptyEXP_pi->GetYaxis()->SetTitle("pT (MeV)");
	ptyEXP_pi->GetXaxis()->SetTitle("y");
	ptyEXP_pi->SetStats(000);
	ptyEXP_pi->Draw("colz");
	
	c->cd(2);
	ptyEXP_p->Divide(ptyHAL_p);
	ptyEXP_p->SetName("p_ratio");
	ptyEXP_p->SetTitle("exp/HAL protons");
	ptyEXP_p->GetYaxis()->SetTitle("pT (MeV)");
	ptyEXP_p->GetXaxis()->SetTitle("y");
	ptyEXP_p->SetStats(000);
	ptyEXP_p->Draw("colz");
	
	
	return 0;
}