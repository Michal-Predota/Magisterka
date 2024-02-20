/*
 * OTFReader.cxx
 *
 *  Created on: 28 maj 2022
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "OTFReader.h"
#include "OTFMcEvent.h"
#include "OTFRecoEvent.h"

#include "Cout.h"
#include "DataManager.h"
#include "Event.h"
#include "McTrack.h"
#include "OTFData.h"
#include "Std.h"
#include "/home/ytcvc/HAL/analysis/spectra/Decay.h"

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TString.h>
#include <TF1.h>
#include <TVector.h>

#include <TRandom.h>

namespace HalOTF {
  Reader::Reader() :
    fSpectras(nullptr),
    fOwner(kFALSE),
    fRegister(kFALSE),
    fMultiplicity(1),
    fPids(211),
    fCharge(1),
    fMass(0),
    fSmear(0.01),
    fMcEvent(nullptr),
    fRecoEvent(nullptr) {}

  Hal::Task::EInitFlag Reader::Init() {
    fMcEvent              = new OTF::McEvent();
    fRecoEvent            = new OTF::RecoEvent();
    Hal::DataManager* mng = Hal::DataManager::Instance();
    if (mng->GetObject("OTF::McEvent.")) {  // branch exists
      fMcEvent   = (OTF::McEvent*) mng->GetObject("OTF::McEvent.");
      fRecoEvent = (OTF::RecoEvent*) mng->GetObject("OTF::RecoEvent.");
    } else {
      fOwner = kTRUE;
      mng->Register("OTF::McEvent.", "HalEvents", fMcEvent, fRegister);
      mng->Register("OTF::RecoEvent.", "HalEvents", fRecoEvent, fRegister);
    }
    return Hal::Task::EInitFlag::kSUCCESS;
  }

  void Reader::SetSpiecies(const TH2D& h, Int_t pid, Double_t w) {
    TH2D* copy         = (TH2D*) h.Clone();
    TDatabasePDG* pdg  = TDatabasePDG::Instance();
    TParticlePDG* part = pdg->GetParticle(pid);
    if (part == nullptr) {
      Hal::Cout::PrintInfo(Form("Cannot add particle with PID = %i", pid), Hal::EInfo::kWarning);
      return;
    }
    copy->SetDirectory(nullptr);
    fSpectras     = copy;
    fPids         = pid;
    fMass         = part->Mass();
    fMultiplicity = w;
    fCharge       = part->Charge() * 3;
  }
  
  
  
	int n_deltas=0;
  void Reader::Exec(Option_t* /*opt*/) {
	  
	  
    PrepareTables();
    Int_t shift = fMcEvent->GetNTracks();
	
	
	int tmp=0;
	int i=0;
	
	
	TF1* curve = new TF1("curve", "1+2*[0]*TMath::Cos(x)", -TMath::Pi(), TMath::Pi());
	
	TF1* breit_wigner = new TF1("breit_wigner", "[0]*(x*x*0.117)/(TMath::Power((1.232*1.232-x*x), 2)+(x*x*0.117*0.117))", 0, 2000);
	breit_wigner->SetParameter(0, 1);
	double norm = breit_wigner->Integral(0, 2000);
	breit_wigner->SetParameter(0, 1/norm);
	
	
	while(i<fMultiplicity)
	{
		
		Double_t Psi = gRandom->Uniform(0, 2*TMath::Pi());
		
		
		Double_t pt, y;
		fSpectras->GetRandom2(y, pt);
	
		//std::cout<<pt/1000<<std::endl;
		
		if(pt/1000<=0.2)
			curve->SetParameter(0, -0.05);
		else if(0.2<pt/1000<=0.25)
			curve->SetParameter(0, -0.075);
		else if(0.25<pt/1000<=0.3)
			curve->SetParameter(0, -0.09);
		else if(0.3<pt/1000<=0.35)
			curve->SetParameter(0, -0.1);
		else if(0.35<pt/1000<=0.4)
			curve->SetParameter(0, -0.12);
		else if(0.4<pt/1000<=0.45)
			curve->SetParameter(0, -0.14);
		else if(0.45<pt/1000<=0.8)
			curve->SetParameter(0, -0.15);
		else
			curve->SetParameter(0, -0.2);
		
		
		
		
		
		Double_t phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
		Double_t check = gRandom->Uniform(0, 1);
		
		//std::cout<<check<<"	"<<curve->Eval(phi)<<"	"<<curve->GetParameter(0)<<std::endl;

	
		if(check<curve->Eval(phi))
		{
			
			Double_t mt  = TMath::Sqrt(pt * pt + fMass * fMass);
			Double_t px  = pt * TMath::Cos(phi);
			Double_t py  = pt * TMath::Sin(phi);
			Double_t pz  = mt * TMath::SinH(y);
	
			 OTF::McTrack tr;
			TLorentzVector p;
			p.SetXYZM(px, py, pz, fMass);
			p.Rotate(Psi, TVector3(1,1,1));
			tr.SetMomentum(p);
			tr.SetPdgCode(fPids);
			TLorentzVector xr(gRandom->Gaus(0, 1), gRandom->Gaus(0, 1), gRandom->Gaus(0), 0);
			tr.SetFreezout(xr);
			fMcEvent->AddTrack(tr);
	
			OTF::RecoTrack rtr;
			px         = px + gRandom->Gaus(0, fSmear) * px;
			py         = py + gRandom->Gaus(0, fSmear) * py;
			pz         = pz + gRandom->Gaus(0, fSmear) * pz;
			Double_t e = TMath::Sqrt(px * px + py * py + pz * pz + fMass * fMass);
			rtr.SetMom(px*1000, py*1000, pz*1000, e*1000);
			
			//std::cout<<"c"<<TLorentzVector(px, py, pz, e).M()<<std::endl;
			
			rtr.SetNHits(5);
			rtr.SetCharge(fCharge);
			rtr.SetMcIndex(i + shift);
			fRecoEvent->AddTrack(rtr);
			tmp=i + shift;
			i++;
		}
	}
	
	 //decay
		Hal::DecayChannel ch1(2212, 211, 10);
		Hal::Decay decay(2224);
		decay.AddDecayChannel(ch1);
		decay.Init();
		std::vector<Hal::McTrack*> tracks;
		for (int i = 0; i < 3; i++) {
		tracks.push_back(new Hal::McTrack());
		}
		Hal::McTrack Delta;
		
		
	while(n_deltas<1000)
	{
		Double_t px = gRandom->Gaus(0, 0.5);
		Double_t py = gRandom->Gaus(0, 0.5);
		Double_t pz = gRandom->Gaus(0, 0.5);
		TLorentzVector sum;
		sum.SetXYZM(px, py, pz, breit_wigner->GetRandom());
		Delta.SetMomentum(px, py, pz, sum.E());
		decay.DecayParticle(Delta, tracks);
		
		
	
		OTF::McTrack tr1, tr2;
		
		TLorentzVector p1, p2;
		p1.SetXYZM(tracks[0]->GetMomentum().Px()*1000, tracks[0]->GetMomentum().Py()*1000, tracks[0]->GetMomentum().Pz()*1000, tracks[0]->GetMomentum().M()*1000);
		tr1.SetMomentum(p1);
		tr1.SetPdgCode(2212);
		TLorentzVector xr1(gRandom->Gaus(0, 1), gRandom->Gaus(0, 1), gRandom->Gaus(0), 0);
		tr1.SetFreezout(xr1);
		fMcEvent->AddTrack(tr1);
		
		
		
		p2.SetXYZM(tracks[1]->GetMomentum().Px()*1000, tracks[1]->GetMomentum().Py()*1000, tracks[1]->GetMomentum().Pz()*1000, tracks[1]->GetMomentum().M()*1000);
		tr2.SetMomentum(p1);
		tr2.SetPdgCode(211);
		TLorentzVector xr2(gRandom->Gaus(0, 1), gRandom->Gaus(0, 1), gRandom->Gaus(0), 0);
		tr2.SetFreezout(xr2);
		fMcEvent->AddTrack(tr2);
		
		
		
		OTF::RecoTrack rtr1, rtr2;
		rtr1.SetMom(tracks[0]->GetMomentum()*1000);
		//std::cout<<"a"<<(tracks[0]->GetMomentum()*1000).M()<<std::endl;
		rtr1.SetNHits(5);
		rtr1.SetCharge(fCharge);
		rtr1.SetMcIndex(n_deltas + tmp);
		n_deltas++;
		fRecoEvent->AddTrack(rtr1);
		
		rtr2.SetMom(tracks[1]->GetMomentum()*1000);
		//std::cout<<"b"<<(tracks[1]->GetMomentum()*1000).M()<<std::endl;
		rtr2.SetNHits(5);
		rtr2.SetCharge(fCharge);
		rtr2.SetMcIndex(n_deltas + tmp);
		fRecoEvent->AddTrack(rtr2);
		//std::cout<<n_deltas<<std::endl;
	}
  }

  Reader::~Reader() {
    if (fSpectras) delete fSpectras;
    if (fOwner && !fRegister) {
      if (fMcEvent) delete fMcEvent;
      if (fRecoEvent) delete fRecoEvent;
    }
  }

  void Reader::PrepareTables() {
    if (fOwner) {
      fRecoEvent->Clear();
      fMcEvent->Clear();
    }
  }

}  // namespace HalOTF