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
#include "TFile.h"
#include "TH2D.h"

#include <TRandom.h>

#include <vector>

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
  
  
  
  
  
void set_curve_params(double pt, int i, TF1* fun, std::vector<double> lim, std::vector<double> val)
{
	//std::cout<<lim.size()<<"	"<<val.size()<<std::endl;
	for(int j=0; j<lim.size(); j++)
	{
		if(pt<lim.at(j))
		{
			fun->SetParameter(i, val.at(j));
			//std::cout<<i<<"	"<<j<<"	"<<val.at(j)<<std::endl;
			//std::cout<<pt<<"	"<<lim.at(j)<<std::endl;
			break;
		}
	}
}	

//TFile* f1 = new TFile("protons_ratio.root");	
//TFile* f2 = new TFile("pions_ratio.root");
	
//TH2D* ratio_p = (TH2D*)f1->Get("p_ratio");
//TH2D* ratio_pi = (TH2D*)f2->Get("pi_ratio"); 
  
	
TH1D* dau = new TH1D("dau","dau",500,1000,1800);

void Reader::Exec(Option_t* /*opt*/) {
	  
TF1* curve = new TF1("curve", "1+2*[0]*TMath::Cos(x)+2*[1]*TMath::Cos(2*x)+2*[2]*TMath::Cos(3*x)+2*[3]*TMath::Cos(4*x)", -TMath::Pi(), TMath::Pi());
TF1* breit_wigner = new TF1("breit_wigner", "[0]*(x*x*0.117)/(TMath::Power((1.232*1.232-x*x), 2)+(x*x*0.117*0.117))", 0, 2);
breit_wigner->SetParameter(0, 1);
double norm = breit_wigner->Integral(0, 2);
breit_wigner->SetParameter(0, 1/norm);

std::vector<double> v1_lim = {0.22083, 0.272916, 0.32083, 0.36875, 0.41875,0.46875,0.51875,0.56,0.614583,0.664583,0.714583,0.7625,0.810416,0.860416,0.9083,0.960417,1.00625,1.0541, 1.1041, 1.15625,1.2083, 1.25416,1.3,
1.352083,1.402083,1.447916,1.497916,1.54583,1.597916,1.64375,1.69583,1.7416,1.797916,1.839583,1.8916,1.939583};

std::vector<double> v1_val = {-0.0574,-0.0744,-0.08510,-0.0978,-0.1085,-0.1191,-0.1191,-0.1361,-0.1361,-0.1425,-0.1446,-0.1510,-0.1553,-0.1617,-0.1638,-0.1659,-0.1680,-0.1702,-0.1744,-0.1787,-0.1787,-0.1808,-0.1851,
-0.1851,-0.1914,-0.1957,-0.1957,-0.2,-0.2,-0.2085,-0.2127,-0.2191,-0.2127,-0.2255,-0.2340};


std::vector<double> v2_lim = {0.2752,0.3268,0.3763,0.4279,0.4752,0.5247,0.5806,0.6279,0.6752,0.7290,0.7741,0.8258,0.8752,0.9290,0.9806,1.0279,1.0752,1.1268,1.1763,1.2279,1.2752,1.3311,1.3806,1.4258,1.4774,1.5225,1.5741,1.6279,
1.6752,1.7225,1.7806,1.8279,1.8752,1.9247,1.9741};

std::vector<double> v2_val = {-0.0256,-0.0330,-0.0422,-0.0504,-0.0568,-0.0660,-0.0779,-0.0899,-0.0981,-0.1128,-0.1256,-0.1366,-0.1486,-0.1623,-0.1733,-0.1834,-0.1954,-0.2027,-0.2100,-0.2165,-0.2220,-0.2275,-0.2339,-0.2394,
-0.2403,-0.2467,-0.2440,-0.2504,-0.2541,-0.2550,-0.2587,-0.2623,-0.2587,-0.2596,-0.2688};


std::vector<double> v3_lim = {0.2244,0.2734,0.3183,0.3693,0.4183,0.4653,0.5142,0.5591,0.6163,0.6591,0.7081,0.7612,0.8081,0.8571,0.9061,0.9571,1.0061,1.0551,1.0979,1.1489,1.2,1.2489,1.2979,1.3448,1.3979,1.4387,1.4959,1.5408,
1.5897,1.6408,1.6857,1.7387,1.7836,1.8326,1.8836,1.9367};

std::vector<double> v3_val = {-0.0006,0.0019,0.0035,0.0041,0.0064,0.0083,0.0103,0.0122,0.0145,0.0177,0.0190,0.0225,0.0254,0.0280,0.0322,0.0345,0.0377,0.0419,0.0441,0.0483,0.0490,0.0512,0.0577,0.0583,0.0593,0.0593,0.0654,0.0619,
0.0709,0.0683,0.0777,0.0835,0.0683,0.0858,0.0796,0.0735};


std::vector<double> v4_lim = {0.2458,0.34583,0.45,0.5458,0.6479,0.7479,0.8458,0.9437,1.0437,1.1458,1.2395,1.3437,1.4395,1.5395,1.6416,1.7375,1.8395,1.9375};

std::vector<double> v4_val = {0,0.0017,0.0028,0.0033,0.0058,0.0073,0.0089,0.0133,0.015,0.0185,0.0219,0.0251,0.0228,0.0216,0.0292,0.0342,0.0344, 0.0337};


std::vector<double> v5_lim = {0.2474,0.3443,0.4391,0.536,0.6288,0.7319,0.8329,0.9278,1.0247,1.1237,1.2185,1.3237,1.4185,1.5154,1.6123,1.7072,1.8061};

std::vector<double> v5_val = {0.0001,0.0002,-0.0007,-0.0001,-0.0014,-0.0004,-0.0018,-0.0022,-0.0034,-0.0043,-0.0030,-0.0128,-0.0043,-0.0162,-0.0079,-0.0409,-0.0144};


std::vector<double> v6_lim = {0.2494,0.3483,0.4451,0.5354,0.6365,0.7419,0.8387,0.9397,1.0344,1.1333,1.2301,1.3333,1.4344,1.5311,1.6215,1.7247,1.8258};

std::vector<double> v6_val = {0,0.0003,0.0001,-0.0009,-0.0016,-0.0010,-0.0014,-0.0051,0.0018,0.00447,-0.00447,-0.002,-0.0016,0.0197,0.0171,-0.040, 0.0440};



	int tmp=0;
	int i=0;
	int which=1;

	  
    PrepareTables();
    Int_t shift = fMcEvent->GetNTracks();
	
	
	
	
	
	
	
	while(i<fMultiplicity)
	{
		
		Double_t Psi = gRandom->Uniform(0, 2*TMath::Pi());
		
		
		Double_t pt, y;
		fSpectras->GetRandom2(y, pt);
	
		

		
		set_curve_params(pt/1000, 0, curve, v1_lim, v1_val);
		set_curve_params(pt/1000, 1, curve, v2_lim, v2_val);
		set_curve_params(pt/1000, 2, curve, v3_lim, v3_val);
		set_curve_params(pt/1000, 3, curve, v4_lim, v4_val);
		
		
		
		
		Double_t phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
		Double_t check = gRandom->Uniform(0, curve->GetMaximum());
		
	
		if(check<curve->Eval(phi))
		{
			
			Double_t mt  = TMath::Sqrt(pt * pt + fMass * fMass);
			Double_t px  = pt * TMath::Cos(phi);
			Double_t py  = pt * TMath::Sin(phi);
			Double_t pz  = mt * TMath::SinH(y);
	
			OTF::McTrack tr;
			TLorentzVector p;
			p.SetXYZM(px, py, pz, fMass);
			p.Rotate(Psi, TVector3(1, 1, 1));
			
			//std::cout<<ratio_p->GetNbinsX()<<std::endl;
			
			
			
			tr.SetMomentum(p);
			tr.SetPdgCode(fPids);
			TLorentzVector xr(gRandom->Gaus(0, 1), gRandom->Gaus(0, 1), gRandom->Gaus(0), 0);
			tr.SetFreezout(xr);
			fMcEvent->AddTrack(tr);
	
			OTF::RecoTrack rtr;
			px         = px + gRandom->Gaus(0, fSmear) * px;
			py         = py + gRandom->Gaus(0, fSmear) * py;
			pz         = pz + gRandom->Gaus(0, fSmear) * pz;
			Double_t e = TMath::Sqrt(px * px + py * py + pz * pz + fMass * fMass*1000*1000);
			rtr.SetMom(px, py, pz, e);
			
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
		Hal::DecayChannel ch1(2212, -211, 10);
		Hal::Decay decay(2224);
		decay.AddDecayChannel(ch1);
		decay.Init();
		std::vector<Hal::McTrack*> tracks;
		for (int i = 0; i < 3; i++) {
		tracks.push_back(new Hal::McTrack());
		}
		Hal::McTrack Delta;
		
	int n_deltas=0;
		
	while(n_deltas<5)
	{
		Double_t px = gRandom->Gaus(0, 0.5);
		Double_t py = gRandom->Gaus(0, 0.5);
		Double_t pz = gRandom->Gaus(0, 0.5);
		TLorentzVector sum;
		double m = breit_wigner->GetRandom();
		sum.SetXYZM(px, py, pz, m);
		Delta.SetMomentum(px, py, pz, sum.E());
		decay.DecayParticle(Delta, tracks, false, m);
		
		if(TMath::IsNaN(tracks[0]->GetMomentum().M())||TMath::IsNaN(tracks[1]->GetMomentum().M()))
			continue;
	
		else{
			OTF::McTrack tr1, tr2;
			OTF::RecoTrack rtr1, rtr2;
			
			double e1 = TMath::Sqrt(tracks[0]->GetMomentum().Px()*1000 * tracks[0]->GetMomentum().Px()*1000 + tracks[0]->GetMomentum().Py()*1000 * tracks[0]->GetMomentum().Py()*1000 + tracks[0]->GetMomentum().Pz()*1000 * tracks[0]->GetMomentum().Pz()*1000 + tracks[0]->GetMomentum().M()*1000 * tracks[0]->GetMomentum().M()*1000);
			double e2 = TMath::Sqrt(tracks[1]->GetMomentum().Px()*1000 * tracks[1]->GetMomentum().Px()*1000 + tracks[1]->GetMomentum().Py()*1000 * tracks[1]->GetMomentum().Py()*1000 + tracks[1]->GetMomentum().Pz()*1000 * tracks[1]->GetMomentum().Pz()*1000 + tracks[1]->GetMomentum().M()*1000 * tracks[1]->GetMomentum().M()*1000);
		
		
			TLorentzVector p1, p2;
			p1.SetXYZM(tracks[0]->GetMomentum().Px()*1000, tracks[0]->GetMomentum().Py()*1000, tracks[0]->GetMomentum().Pz()*1000, e1);
		
			p2.SetXYZM(tracks[1]->GetMomentum().Px()*1000, tracks[1]->GetMomentum().Py()*1000, tracks[1]->GetMomentum().Pz()*1000, e2);
			
			
			tr1.SetMomentum(p2);
			tr1.SetPdgCode(-211);
			TLorentzVector xr1(gRandom->Gaus(0, 1), gRandom->Gaus(0, 1), gRandom->Gaus(0), 0);
			tr1.SetFreezout(xr1);
			fMcEvent->AddTrack(tr1);
			
			rtr2.SetMom(tracks[1]->GetMomentum()*1000);
			//std::cout<<"b"<<(tracks[1]->GetMomentum()*1000).M()<<std::endl;
			rtr2.SetNHits(5);
			rtr2.SetCharge(fCharge);
			rtr2.SetMcIndex(which + tmp);
			fRecoEvent->AddTrack(rtr2);
			//std::cout<<"tmp	"<<tmp<<std::endl;
			//std::cout<<"which	"<<which<<std::endl;
			//std::cout<<which + tmp<<std::endl;
			
			which++;
			//std::cout<<tracks[0]->GetMomentum().M()<<std::endl;
			
			//tracks[0]=proton, tracks[1]=pion
			
			tr2.SetMomentum(p1);
			tr2.SetPdgCode(2212);
			TLorentzVector xr2(gRandom->Gaus(0, 1), gRandom->Gaus(0, 1), gRandom->Gaus(0), 0);
			tr2.SetFreezout(xr2);
			fMcEvent->AddTrack(tr2);
			
			rtr1.SetMom(tracks[0]->GetMomentum()*1000);
			//std::cout<<"a"<<(tracks[0]->GetMomentum()*1000).M()<<std::endl;
			rtr1.SetNHits(5);
			rtr1.SetCharge(fCharge);
			rtr1.SetMcIndex(which + tmp);
			fRecoEvent->AddTrack(rtr1);
			//std::cout<<"tmp	"<<tmp<<std::endl;
			//std::cout<<"which	"<<which<<std::endl;
			//std::cout<<which + tmp<<std::endl;
			which++;
			
			n_deltas++;
		//	tmp=tmp+which;
			
			//std::cout<<((tracks[0]->GetMomentum()+tracks[1]->GetMomentum()).M()*1000)<<std::endl;	
		}
	}
	
	//dau->SaveAs("/mnt/c/Users/bumcy/Desktop/dau.root");
	
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