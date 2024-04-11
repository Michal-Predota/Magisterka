/*
 * starflow.C
 *
 *  Created on: 22 paź 2023
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

/*
 * starflow.C
 *
 *  Created on: 1 sie 2023
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek@gmail.com
 *      Warsaw University of Technology, Faculty of Physics
 */


#include <TFile.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#ifndef __CLING__
#include "Header.h"
#include "Decay.h"

#endif

class FlowGenerator : public TObject {
  Double_t fMag;

public:
  FlowGenerator() : fMag(0) {};
  void SetFlow(Double_t mag) { fMag = mag; }
  virtual ~FlowGenerator() {}

  virtual Double_t GenerateWeight(Hal::TwoTrack* pair) {
    auto rphi = [](Double_t dphi) {
      if (dphi < -0.5 * TMath::Pi()) dphi += TMath::TwoPi();
      if (dphi > 1.5 * TMath::Pi()) dphi -= TMath::TwoPi();
      return dphi;
    };
    Double_t phi1   = pair->GetTrack1()->GetMomentum().Phi();
    Double_t phi2   = pair->GetTrack2()->GetMomentum().Phi();
    Double_t eta1   = pair->GetTrack1()->GetMomentum().Eta();
    Double_t eta2   = pair->GetTrack2()->GetMomentum().Eta();
    Double_t dphi   = phi1 - phi2;
    Double_t weight = 1.0 + fMag * TMath::Cos(dphi);
    Double_t deta   = eta1 - eta2;
    return weight;
  }
};

class MinvAna : public Hal::TwoTrackAna {
  TH1D* fMinvNum;
  TH1D* fMinvDen;
  TH1D* phi_pion;
  TH1D* phi_proton;
  TH2D* pty_proton;
  TH2D* pty_pion;
  
  FlowGenerator* fFlow;

protected:
  virtual Hal::Package* Report() const {  // generuje zapis do pliku
    auto report = Hal::TwoTrackAna::Report();
    report->AddObject(fMinvNum);
    report->AddObject(fMinvDen);
    report->AddObject(phi_pion);
    report->AddObject(phi_proton);
    report->AddObject(pty_proton);
    report->AddObject(pty_pion);
    return report;
  }
 
 
 /*
TFile* f1 = new TFile("protons_ratio.root");	
TFile* f2 = new TFile("pions_ratio.root");
	
TH2D* ratio_p = (TH2D*)f1->Get("p_ratio");
TH2D* ratio_pi = (TH2D*)f2->Get("pi_ratio"); 
  */
 virtual void ProcessPair() {  // wywoływane dla każdej pary zwykłej
 
 //TCanvas *mydummycanvas=new TCanvas();



	//cout<<ratio_pi->GetNbinsX()<<endl;
 
	Double_t minv;
	Double_t weight;
    TLorentzVector mom1 = fCurrentSignalPair->GetTrack1()->GetMomentum();
    TLorentzVector mom2 = fCurrentSignalPair->GetTrack2()->GetMomentum();
	//cout<<mom1.M()<<endl;
	//cout<<mom2.M()<<endl;
	
	
	//if(TMath::Abs(mom1.M()-mom2.M())>1)
	//{
		/*
		if(mom1.M()<500)
			mom1*(ratio_pi->GetBinContent(ratio_pi->FindBin(mom1.Rapidity(), mom1.Pt())));
		if(mom1.M()>500)
			mom1*(ratio_p->GetBinContent(ratio_p->FindBin(mom1.Rapidity(), mom1.Pt())));
		
		if(mom2.M()<500)
			mom2*(ratio_pi->GetBinContent(ratio_pi->FindBin(mom2.Rapidity(), mom2.Pt())));
		if(mom2.M()>500)
			mom2*(ratio_p->GetBinContent(ratio_p->FindBin(mom2.Rapidity(), mom2.Pt())));
		*/
		minv       = (mom1 + mom2).M();
		weight     = fFlow->GenerateWeight(fCurrentSignalPair);
	
		fMinvNum->Fill(minv, weight);
		//cout<<minv<<endl;
		
	//	if(fCurrentSignalPair->GetTrack1()->GetMomentum().M()<200)
		//{
			phi_pion->Fill(fCurrentSignalPair->GetTrack2()->GetMomentum().Phi());
			pty_pion->Fill(fCurrentSignalPair->GetTrack2()->GetMomentum().Rapidity(), fCurrentSignalPair->GetTrack1()->GetMomentum().Pt());
		//}
		
		//if(fCurrentSignalPair->GetTrack2()->GetMomentum().M()<200)
	//	{
			phi_pion->Fill(fCurrentSignalPair->GetTrack1()->GetMomentum().Phi());
			pty_pion->Fill(fCurrentSignalPair->GetTrack1()->GetMomentum().Rapidity(), fCurrentSignalPair->GetTrack1()->GetMomentum().Pt());
	//	}
		
		/*if(fCurrentSignalPair->GetTrack2()->GetMomentum().M()>200)
		{
			phi_proton->Fill(fCurrentSignalPair->GetTrack2()->GetMomentum().Phi());
			pty_proton->Fill(fCurrentSignalPair->GetTrack2()->GetMomentum().Rapidity(), fCurrentSignalPair->GetTrack2()->GetMomentum().Pt());
		}
		
		if(fCurrentSignalPair->GetTrack1()->GetMomentum().M()>200)
		{
			phi_proton->Fill(fCurrentSignalPair->GetTrack1()->GetMomentum().Phi());
			pty_proton->Fill(fCurrentSignalPair->GetTrack1()->GetMomentum().Rapidity(), fCurrentSignalPair->GetTrack1()->GetMomentum().Pt());
		}*/
	//}
	
//	f1->Close();
	//f2->Close();
	
  }
  
  
  virtual void ProcessPair_Perfect() {  // wywoływane dla każdej zwykłej (ale bez wag od flow)
	Double_t minv;
	Double_t weight;
    TLorentzVector mom1 = fCurrentBackgroundPair->GetTrack1()->GetMomentum();
    TLorentzVector mom2 = fCurrentBackgroundPair->GetTrack2()->GetMomentum();
	
	if(TMath::Abs(mom1.M()-mom2.M())>1)
	{
		minv       = (mom1 + mom2).M();
		
		//std::cout<<minv<<std::endl;
		weight     = fFlow->GenerateWeight(fCurrentSignalPair);
	
		fMinvDen->Fill(minv);
	}
	
  }

public:
  MinvAna(Int_t bins, Double_t min, Double_t max) : fFlow(nullptr) {
    fMinvNum = new TH1D("num", "num:m_{Inv} [MeV/c];m", bins, min, max);
    fMinvDen = new TH1D("den", "den;m_{Inv} [MeV/c];m", bins, min, max);
	phi_pion = new TH1D("phi pion", "phi, pion", 100, -2*TMath::Pi(), 2*TMath::Pi());
	phi_proton = new TH1D("phi proton", "phi, proton", 100, -2*TMath::Pi(), 2*TMath::Pi());
	pty_pion = new TH2D("pty pion", "pty pion", 100, 0, 2, 100, 0, 1500);
	pty_proton = new TH2D("pty proton", "pty proton", 100, 0, 2, 100, 0, 1500);
	
  };
  void SetWeight(FlowGenerator* flow) { fFlow = flow; }
  virtual ~MinvAna() {};
};


void minflow() {
  gStyle->SetPalette(kRainBow);
  Hal::AnalysisManager* run = new Hal::AnalysisManager();
  HalOTF::Source* source    = new HalOTF::Source(10000);
  /**wczytywanie spektr i ich ustawianie**/
  TString path = "spec_pim.root";
  TFile* fx    = new TFile(path);
  TFile *fp = new TFile("spec_p.root");
  
  auto w       = new FlowGenerator();
  w->SetFlow(0.4);
  TH2D* h                = (TH2D*) fx->Get("HALpim");
  TH2D *hp = (TH2D*) fp->Get("HALp");
  
  

  
  
  HalOTF::Reader* reader = new HalOTF::Reader();
  reader->SetSpiecies(*hp, 2212, 100);
  run->SetSource(source);
  run->AddTask(reader);
  
  reader = new HalOTF::Reader();
  reader->SetSpiecies(*h, -211, 100);
  run->AddTask(reader);
  // analiza
  MinvAna* ana = new MinvAna(1000, 1000, 2500);
  ana->SetFormat(new HalOTF::ComplexEvent());


  Hal::TrackPdgCut pid1, pid2;
  pid1.SetMinAndMax(2212);
  pid2.SetMinAndMax(-211);
  ana->AddCut(pid1, "{0}+im");
  ana->AddCut(pid2, "{1}+im");
  
  ana->EnableNonIdentical();
  
  

  ana->SetWeight(w);
  ana->SetOption(Hal::TwoTrackAna::BackgroundOptionPerfect());
  run->AddTask(ana);
  run->SetOutput("flow_minv.root");
  run->Init();

  run->Run();
}
