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
	
	
	TH1D* h = new TH1D("h","h", 500, 0, 2000);
	
	
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
	
	TCanvas* c = new TCanvas("c","c",800,800);
	c->cd();
	dau->Draw();*/
	//breit_wigner->Draw("SAME");
	
	TFile* f1 = new TFile("deltas_noflow.root");
	TFile* f2 = new TFile("nodeltas_noflow.root");
	
	TH1D* deltas = (TH1D*)f1->Get("den");
	TH1D* nodeltas = (TH1D*)f2->Get("den");
	
	deltas->Add(nodeltas, -1);

	deltas->Draw();
	
	
	return 0;
}