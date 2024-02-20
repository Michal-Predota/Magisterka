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
  
  FlowGenerator* fFlow;

protected:
  virtual Hal::Package* Report() const {  // generuje zapis do pliku
    auto report = Hal::TwoTrackAna::Report();
    report->AddObject(fMinvNum);
    report->AddObject(fMinvDen);
    report->AddObject(phi_pion);
    report->AddObject(phi_proton);
    return report;
  }
  
  virtual void ProcessPair() {  // wywoływane dla każdej pary zwykłej
	Double_t minv;
	Double_t weight;
    TLorentzVector mom1 = fCurrentSignalPair->GetTrack1()->GetMomentum();
    TLorentzVector mom2 = fCurrentSignalPair->GetTrack2()->GetMomentum();
	
	if(mom1.M()>100&&mom1.M()<200&&mom2.M()>900)
	{
		minv       = (mom1 + mom2).M();
		weight     = fFlow->GenerateWeight(fCurrentSignalPair);
	
		fMinvNum->Fill(minv, weight);
		phi_pion->Fill(fCurrentSignalPair->GetTrack1()->GetMomentum().Phi());
		phi_proton->Fill(fCurrentSignalPair->GetTrack2()->GetMomentum().Phi());
	}
	
	
//	std::cout<<fCurrentSignalPair->GetTrack1()->GetMomentum().M()<<"	"<<fCurrentSignalPair->GetTrack2()->GetMomentum().M()<<std::endl;
	
	
  }
  
  
  virtual void ProcessPair_Perfect() {  // wywoływane dla każdej zwykłej (ale bez wag od flow)
	Double_t minv;
	Double_t weight;
    TLorentzVector mom1 = fCurrentBackgroundPair->GetTrack1()->GetMomentum();
    TLorentzVector mom2 = fCurrentBackgroundPair->GetTrack2()->GetMomentum();
	
	if(mom1.M()>100&&mom1.M()<200&&mom2.M()>900)
	{
		minv       = (mom1 + mom2).M();
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
  };
  void SetWeight(FlowGenerator* flow) { fFlow = flow; }
  virtual ~MinvAna() {};
};


void minflow() {
  gStyle->SetPalette(kRainBow);
  Hal::AnalysisManager* run = new Hal::AnalysisManager();
  HalOTF::Source* source    = new HalOTF::Source(10000);
  /**wczytywanie spektr i ich ustawianie**/
  TString path = "spec_pip.root";
  TFile* fx    = new TFile(path);
  TFile *fp = new TFile("spec_p.root");
  
  auto w       = new FlowGenerator();
  w->SetFlow(0.4);
  TH2D* h                = (TH2D*) fx->Get("HALpip");
  TH2D *hp = (TH2D*) fp->Get("HALp");
  
  
 
  
  
  
  HalOTF::Reader* reader = new HalOTF::Reader();
  reader->SetSpiecies(*hp, 2212, 100);
  run->SetSource(source);
  run->AddTask(reader);
  
  reader = new HalOTF::Reader();
  reader->SetSpiecies(*h, 211, 100);
  run->AddTask(reader);
  // analiza
  MinvAna* ana = new MinvAna(500, 1000, 1500);
  ana->SetFormat(new HalOTF::ComplexEvent());


  Hal::TrackPdgCut pid1, pid2;
  pid1.SetValue(2212);
  pid2.SetValue(211);
  ana->AddCut(pid1, "{0}");
  ana->AddCut(pid2, "{1}");
  
  ana->EnableNonIdentical();
  
  ana->SetOption("id");
  

  ana->SetWeight(w);
  ana->SetOption(Hal::TwoTrackAna::BackgroundOptionPerfect());
  run->AddTask(ana);
  run->SetOutput("flow_minv.root");
  run->Init();

  run->Run();
}
