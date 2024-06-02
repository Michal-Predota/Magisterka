
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
  
  
  
  
 
	
Int_t set_curve_params(double pt, int i, TF1* fun, std::vector<double> lim, std::vector<double> val)
{
	//std::cout<<lim.size()<<"	"<<val.size()<<std::endl;
	for(int j=0; j<lim.size(); j++)
	{
		if(pt<lim.at(j))
		{
			fun->SetParameter(i, val.at(j));
			//std::cout<<i<<"	"<<j<<"	"<<val.at(j)<<std::endl;
			//std::cout<<pt<<"	"<<lim.at(j)<<std::endl;
			return val.at(j);
		}
	}
}	


Int_t set_flow(double y, double pt, int v_n, TF1* fun)//ustawianie flow zgodnie z parametrami od Behruza
{
	 //v1

std::vector<double> y_neg005_005_v1_pt ={0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_neg005_005_v1 ={0.00555, -0.00572, -0.0123, -0.016, -0.0162, -0.0152, -0.0136, -0.0124, -0.011, -0.00975, -0.00831, -0.00671, -0.00564, -0.00541, -0.00448, -0.00411, -0.00467, -0.00485, -0.00335, -0.00358, -0.00233, -0.00124, -0.000238, -0.00163, -0.00197, -0.0022, -0.00229, -0.00385, -0.00486, -0.00726, -0.00531, -0.0113, -0.0107, -0.0102, -0.0101};

std::vector<double> y_015_025_v1_pt ={0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_015_025_v1 ={-0.0536, -0.0658, -0.0661, -0.0682, -0.0719, -0.0769, -0.0821, -0.0879, -0.0937, -0.0992, -0.104, -0.109, -0.112, -0.116, -0.119, -0.122, -0.124, -0.125, -0.127, -0.128, -0.128, -0.128, -0.128, -0.129, -0.129, -0.128, -0.127, -0.123, -0.124, -0.129, -0.124, -0.125, -0.126, -0.129};

std::vector<double> y_neg025_015_v1_pt ={0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_neg025_neg015_v1 ={-0.0427, -0.0585, -0.0678, -0.0744, -0.0803, -0.0846, -0.0884, -0.0933, -0.0961, -0.097, -0.097, -0.0982, -0.0997, -0.102, -0.105, -0.108, -0.11, -0.113, -0.115, -0.117, -0.119, -0.122, -0.124, -0.125, -0.128, -0.128, -0.132, -0.132, -0.134, -0.139, -0.141, -0.144, -0.147, -0.148, -0.155, -0.148};

std::vector<double> y_035_045_v1_pt ={0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_035_045_v1 ={-0.139, -0.169, -0.181, -0.191, -0.2, -0.208, -0.215, -0.223, -0.23, -0.236, -0.242, -0.245, -0.248, -0.25, -0.252, -0.254, -0.256, -0.257, -0.258, -0.258, -0.259, -0.258, -0.262, -0.257, -0.258, -0.257, -0.258, -0.253, -0.261, -0.244, -0.253, -0.258};

std::vector<double> y_neg045_neg035_v1_pt ={0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_neg045_neg035_v1 ={-0.0841, -0.0974, -0.109, -0.119, -0.127, -0.133, -0.142, -0.154, -0.166, -0.178, -0.188, -0.197, -0.205, -0.212, -0.218, -0.224, -0.228, -0.233, -0.237, -0.24, -0.244, -0.247, -0.25, -0.254, -0.256, -0.259, -0.262, -0.26, -0.267, -0.27, -0.27, -0.276, -0.281, -0.277, -0.286, -0.288, -0.294};

std::vector<double> y_055_065_v1_pt ={0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775};
std::vector<double> y_055_065_v1 ={-0.281, -0.314, -0.333, -0.343, -0.352, -0.358, -0.362, -0.367, -0.371, -0.374, -0.376, -0.378, -0.38, -0.382, -0.382, -0.38, -0.38, -0.382, -0.375, -0.378, -0.378, -0.377, -0.375, -0.374, -0.378};

std::vector<double> y_neg065_neg055_v1_pt ={0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775};
std::vector<double> y_neg065_neg055_v1 ={-0.152, -0.168, -0.186, -0.203, -0.222, -0.24, -0.258, -0.274, -0.288, -0.301, -0.312, -0.322, -0.331, -0.338, -0.344, -0.349, -0.354, -0.357, -0.361, -0.364, -0.366, -0.369, -0.372, -0.375, -0.379, -0.379, -0.38, -0.382, -0.383, -0.383, -0.383, -0.386};

//v2

std::vector<double> y_neg005_005v2_pt ={0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_neg005_005v2 ={-0.0167, -0.0234, -0.0262, -0.0272, -0.0291, -0.032, -0.0347, -0.0376, -0.04, -0.0431, -0.0464, -0.0494, -0.0535, -0.0573, -0.0599, -0.0637, -0.0659, -0.0685, -0.0712, -0.0741, -0.076, -0.081, -0.0811, -0.0836, -0.081, -0.0844, -0.0866, -0.0844, -0.0904, -0.0906, -0.0911, -0.0872, -0.0966, -0.0886, -0.102};

std::vector<double> y_015_025v2_pt ={0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_015_025v2 ={-0.0159, -0.0191, -0.0245, -0.03, -0.0338, -0.036, -0.0385, -0.0405, -0.0434, -0.046, -0.0496, -0.0518, -0.0543, -0.0565, -0.0591, -0.0616, -0.0636, -0.065, -0.0669, -0.0683, -0.0684, -0.0727, -0.073, -0.075, -0.0773, -0.0821, -0.0814, -0.0838, -0.0734, -0.0694, -0.0844, -0.105, -0.0709, -0.0855};

std::vector<double> y_neg025_neg015v2_pt ={0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_neg025_neg015v2 ={-0.0182, -0.0188, -0.0175, -0.0177, -0.0181, -0.019, -0.0203, -0.023, -0.0265, -0.0292, -0.0334, -0.0367, -0.0399, -0.0429, -0.0463, -0.0495, -0.052, -0.0545, -0.0584, -0.06, -0.0614, -0.0641, -0.0655, -0.0678, -0.0684, -0.0724, -0.0721, -0.0741, -0.0756, -0.0716, -0.076, -0.0675, -0.0711, -0.0543, -0.076, -0.0801};

std::vector<double> y_035_045v2_pt ={0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_035_045v2 ={-0.013, -0.011, -0.0134, -0.0155, -0.0183, -0.0205, -0.0218, -0.0233, -0.0237, -0.0234, -0.0253, -0.0256, -0.0258, -0.0291, -0.0311, -0.0308, -0.0352, -0.034, -0.0363, -0.0391, -0.0403, -0.0414, -0.0413, -0.0536, -0.0403, -0.0483, -0.0426, -0.0517, -0.0473, -0.0425, -0.0556, -0.0285};

std::vector<double> y_neg045_neg035v2_pt ={0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_neg045_neg035v2 ={-0.0128, -0.00907, -0.00525, -0.00504, -0.00627, -0.00861, -0.00999, -0.012, -0.0124, -0.0135, -0.015, -0.0159, -0.0175, -0.0184, -0.0201, -0.0223, -0.0234, -0.0249, -0.0254, -0.0286, -0.0283, -0.0298, -0.0316, -0.0322, -0.0293, -0.0285, -0.03, -0.034, -0.0298, -0.0322, -0.037, -0.031, -0.0242, -0.0343, -0.0259, -0.0302, -0.024};

std::vector<double> y_055_065v2_pt ={0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775};
std::vector<double> y_055_065v2 ={0.0128, 0.0244, 0.0289, 0.0318, 0.0329, 0.0331, 0.032, 0.0311, 0.0309, 0.0289, 0.0294, 0.0261, 0.0268, 0.0312, 0.0306, 0.0226, 0.0285, 0.0222, 0.0243, 0.022, 0.0236, 0.0227, 0.00524, 0.0449, 0.012};

std::vector<double> y_neg065_neg055v2_pt ={0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775};
std::vector<double> y_neg065_neg055v2 ={0.00364, 0.00339, 0.00212, 0.00205, 0.00238, 0.00404, 0.0068, 0.00962, 0.013, 0.0157, 0.0184, 0.0202, 0.0216, 0.0238, 0.0246, 0.0244, 0.0239, 0.0266, 0.0259, 0.0257, 0.0248, 0.0242, 0.0248, 0.0362, 0.0247, 0.017, 0.0225, 0.0334, 0.0329, 0.0301, 0.00759, 0.0247};

//v3

std::vector<double> y_neg005_005v3_pt ={0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_neg005_005v3 ={-0.0029, -0.00517, -0.00266, -0.00249, 0.00024, 0.00164, 0.00041, 0.000784, 0.000637, 0.00022, 0.000448, 0.00164, 0.000167, 0.00157, -0.000136, 0.00317, -0.000127, 0.00206, -0.00037, 0.00277, 0.00451, -0.000922, -0.00511, 0.0015, 0.001, -0.0017, 0.00331, 0.00959, -0.0043, -0.00208, 0.00201, 0.0201, -0.0152, 0.0219, -0.013};

std::vector<double> y_015_025v3_pt ={0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_015_025v3 ={0.00592, 0.00873, 0.00982, 0.0109, 0.0109, 0.00895, 0.0085, 0.00761, 0.00984, 0.0123, 0.0121, 0.013, 0.015, 0.0175, 0.0183, 0.0197, 0.0195, 0.019, 0.0226, 0.0197, 0.0227, 0.0216, 0.0226, 0.0316, 0.016, 0.00987, 0.014, 0.0229, 0.011, 0.0239, 0.0373, 0.0111, 0.0472, 0.0101};

std::vector<double> y_neg025_neg015v3_pt ={0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_neg025_neg015v3 ={-0.0021, -0.000588, 0.00184, 0.00364, 0.0044, 0.00486, 0.00561, 0.00663, 0.00596, 0.00683, 0.00614, 0.00802, 0.00901, 0.00869, 0.00948, 0.0108, 0.00993, 0.0131, 0.0157, 0.018, 0.0181, 0.0219, 0.0125, 0.0206, 0.026, 0.0223, 0.0229, 0.0224, 0.0134, 0.0183, 0.0224, 0.0304, 0.0296, 0.0305, 0.0306, 0.0569};

std::vector<double> y_035_045v3_pt ={0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_035_045v3 ={0.00907, 0.0114, 0.013, 0.0157, 0.0209, 0.0217, 0.0225, 0.0259, 0.0277, 0.028, 0.0295, 0.0259, 0.0277, 0.0331, 0.0347, 0.0343, 0.0508, 0.035, 0.0353, 0.0408, 0.0424, 0.034, 0.0408, 0.0469, 0.0422, 0.039, 0.0436, 0.0332, 0.0752, 0.0199, 0.0745, 0.0741};

std::vector<double> y_neg045_neg035v3_pt ={0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775, 1.825, 1.875, 1.925, 1.975};
std::vector<double> y_neg045_neg035v3 ={0.000701, 0.0056, 0.00256, 0.00474, 0.00433, 0.00364, 0.00557, 0.00737, 0.00663, 0.00772, 0.00936, 0.0108, 0.0129, 0.0118, 0.0154, 0.0197, 0.0207, 0.021, 0.0245, 0.0244, 0.0255, 0.0278, 0.0252, 0.0223, 0.0356, 0.0318, 0.0312, 0.0345, 0.037, 0.04, 0.0308, 0.0406, 0.0169, 0.0467, 0.051, 0.0863, -0.00637};

std::vector<double> y_055_065v3_pt ={0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775};
std::vector<double> y_055_065v3 ={0.0128, 0.016, 0.0183, 0.0204, 0.0228, 0.0255, 0.0283, 0.0346, 0.0402, 0.0361, 0.0353, 0.0501, 0.038, 0.0265, 0.0493, 0.0514, 0.0492, 0.0736, 0.0225, 0.0383, 0.044, 0.0697, 0.0361, 0.0242, 0.0444};

std::vector<double> y_neg065_neg055v3_pt ={0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.075, 1.125, 1.175, 1.225, 1.275, 1.325, 1.375, 1.425, 1.475, 1.525, 1.575, 1.625, 1.675, 1.725, 1.775};
std::vector<double> y_neg065_neg055v3 ={0.000451, 0.00314, 0.0059, 0.00509, 0.00456, 0.00623, 0.00598, 0.00799, 0.00939, 0.0101, 0.0117, 0.0151, 0.0142, 0.0169, 0.0209, 0.0197, 0.0225, 0.0238, 0.023, 0.0251, 0.034, 0.0314, 0.0383, 0.0237, 0.0359, 0.0572, 0.05, 0.0238, 0.03, 0.0503, 0.00463, 0.0479};

//v4

std::vector<double> y_neg005_005v4_pt ={0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95};
std::vector<double> y_neg005_005v4 ={-0.00326, 0.00241, 0.00181, 0.00196, 0.004, 0.00136, 0.00118, 0.00104, 0.000172, -0.00013, 0.00202, 0.0119, 0.0184, 0.0096, 0.0365, 0.0113, 0.0451, -0.0252};

std::vector<double> y_015_025v4_pt ={0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95};
std::vector<double> y_015_025v4 ={-0.00227, 0.000673, 0.0024, 0.00487, 0.00148, -0.00102, 0.00167, 0.00078, 0.00613, -0.00885, 0.0107, -0.00334, 0.0192, 0.0257, -0.0199, 0.00122, -0.041};

std::vector<double> y_neg025_neg015v4_pt ={0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95};
std::vector<double> y_neg025_neg015v4 ={-0.00085, 0.000197, 0.000188, -0.000722, 0.000948, 0.000792, 0.000395, 0.002, -0.000969, -0.000684, 0.0013, 0.00629, 0.00211, 0.0182, -0.0271, 0.00902, -0.029, 0.00506};

std::vector<double> y_035_045v4_pt ={0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95};
std::vector<double> y_035_045v4 ={-0.000368, -0.0074, -0.00254, 0.00158, -0.00371, -0.0118, -7.94e-05, 0.000621, -0.0154, -0.00687, 0.0111, -0.0158, -0.0468, -0.0864, 0.0212, -0.0523};

std::vector<double> y_neg045_neg035v4_pt ={0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95};
std::vector<double> y_neg045_neg035v4 ={-0.000794, -0.00204, -0.00107, 0.000155, -0.00445, -0.00257, -0.00354, -0.0048, -0.00592, -0.00392, -0.00354, -0.00902, 0.00355, -0.00323, 0.0214, -0.0211, -0.041, -0.0576, -0.0439};

std::vector<double> y_055_065v4_pt ={0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75};
std::vector<double> y_055_065v4 ={-0.00486, -0.00289, -0.0231, -0.0168, -0.0157, -0.00463, -0.0374, 0.0143, -0.0724, -0.0139, -0.000916, 0.000751, -0.0705};

std::vector<double> y_neg065_neg055v4_pt ={0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75};
std::vector<double> y_neg065_neg055v4 ={-0.00474, -0.00489, -0.00148, 0.000285, -0.000415, -0.00416, -0.0033, -0.00407, -0.00896, -0.0114, -0.009, 0.00282, 0.0347, -0.0659, -0.0369, -0.0304};

	
	
	
	
	//v1
	if(v_n==1)
	{
		
		if(-0.05<y-1 && y-1<0.05)
		{
			for(int i=0; i<y_neg005_005_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_neg005_005_v1_pt.at(i)-pt)<TMath::Abs(y_neg005_005_v1_pt.at(1)-y_neg005_005_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, y_neg005_005_v1.at(i));
					//std::cout<<y_neg005_005_v1.at(i)<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						
					}
				}
			}
		}
		
		if(0.05<y-1 && y-1<0.15)
		{
			for(int i=0; i<y_neg005_005_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_neg005_005_v1_pt.at(i)-pt)<TMath::Abs(y_neg005_005_v1_pt.at(1)-y_neg005_005_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, (y_neg005_005_v1.at(i)+y_015_025_v1.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.15<y-1 && y-1<0.25)
		{
			for(int i=0; i<y_015_025_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_015_025_v1_pt.at(i)-pt)<TMath::Abs(y_015_025_v1_pt.at(1)-y_015_025_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, y_015_025_v1.at(i));
					//std::cout<<y_015_025_v1.at(i)<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.25<y-1 && y-1<0.35)
		{
			for(int i=0; i<y_035_045_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_035_045_v1_pt.at(i)-pt)<TMath::Abs(y_035_045_v1_pt.at(1)-y_035_045_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, (y_015_025_v1.at(i)+y_035_045_v1.at(i))/2);
					//std::cout<<(y_015_025_v1.at(i)+y_035_045_v1.at(i))/2<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.35<y-1 && y-1<0.45)
		{
			for(int i=0; i<y_035_045_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_035_045_v1_pt.at(i)-pt)<TMath::Abs(y_035_045_v1_pt.at(1)-y_035_045_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, y_035_045_v1.at(i));
					//std::cout<<y_035_045_v1.at(i)<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.45<y-1 && y-1<0.55)
		{
			for(int i=0; i<y_055_065_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_055_065_v1_pt.at(i)-pt)<TMath::Abs(y_055_065_v1_pt.at(1)-y_055_065_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, (y_035_045_v1.at(i)+y_055_065_v1.at(i))/2);
					//std::cout<<(y_035_045_v1.at(i)+y_055_065_v1.at(i))/2<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.55<y-1 && y-1<0.65)
		{
			for(int i=0; i<y_055_065_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_055_065_v1_pt.at(i)-pt)<TMath::Abs(y_055_065_v1_pt.at(1)-y_055_065_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, y_055_065_v1.at(i));
					//std::cout<<y_055_065_v1.at(i)<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		
		
		
		
		if(-0.15<y-1 && y-1<-0.05)
		{
			for(int i=0; i<y_neg025_015_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_015_v1_pt.at(i)-pt)<TMath::Abs(y_neg025_015_v1_pt.at(1)-y_neg025_015_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, (y_neg025_neg015_v1.at(i)+y_neg005_005_v1.at(i))/2);
					//std::cout<<(y_neg025_neg015_v1.at(i)+y_neg005_005_v1.at(i))/2<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		
		if(-0.25<y-1 && y-1<-0.15)
		{
			for(int i=0; i<y_neg025_015_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_015_v1_pt.at(i)-pt)<TMath::Abs(y_neg025_015_v1_pt.at(1)-y_neg025_015_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, y_neg025_neg015_v1.at(i));
					//std::cout<<y_neg025_neg015_v1.at(i)<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.35<y-1 && y-1<-0.25)
		{
			for(int i=0; i<y_neg025_015_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_015_v1_pt.at(i)-pt)<TMath::Abs(y_neg025_015_v1_pt.at(1)-y_neg025_015_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, (y_neg025_neg015_v1.at(i)+y_neg045_neg035_v1.at(i))/2);
					//std::cout<<(y_neg025_neg015_v1.at(i)+y_neg045_neg035_v1.at(i))/2<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.45<y-1 && y-1<-0.35)
		{
			for(int i=0; i<y_neg045_neg035_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_neg045_neg035_v1_pt.at(i)-pt)<TMath::Abs(y_neg045_neg035_v1_pt.at(1)-y_neg045_neg035_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, y_neg045_neg035_v1.at(i));
					//std::cout<<y_neg045_neg035_v1.at(i)<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.55<y-1 && y-1<-0.45)
		{
			for(int i=0; i<y_neg045_neg035_v1_pt.size(); i++)
			{
				//std::cout<<y_neg045_neg035_v1.size()<<" "<<y_neg065_neg055_v1.size()<<std::endl;
				if(TMath::Abs(y_neg045_neg035_v1_pt.at(i)-pt)<TMath::Abs(y_neg045_neg035_v1_pt.at(1)-y_neg045_neg035_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, (y_neg045_neg035_v1.at(i)+y_neg065_neg055_v1.at(i))/2);
					//std::cout<<(y_neg045_neg035_v1.at(i)+y_neg065_neg055_v1.at(i))/2<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		
		if(-0.65<y-1 && y-1<-0.55)
		{
			for(int i=0; i<y_neg065_neg055_v1_pt.size(); i++)
			{
				if(TMath::Abs(y_neg065_neg055_v1_pt.at(i)-pt)<TMath::Abs(y_neg065_neg055_v1_pt.at(1)-y_neg065_neg055_v1_pt.at(2))/2)
				{
					try{
					fun->SetParameter(0, y_neg065_neg055_v1.at(i));
					//std::cout<<y_neg065_neg055_v1.at(i)<<std::endl;
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
	}
	//v2
	if(v_n==2)
	{
		if(-0.05<y-1 && y-1<0.05)
		{
			for(int i=0; i<y_neg005_005v2_pt.size(); i++)
			{
				if(TMath::Abs(y_neg005_005v2_pt.at(i)-pt)<TMath::Abs(y_neg005_005v2_pt.at(1)-y_neg005_005v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, y_neg005_005v2.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.05<y-1 && y-1<0.15)
		{
			for(int i=0; i<y_neg005_005v2_pt.size(); i++)
			{
				if(TMath::Abs(y_neg005_005v2_pt.at(i)-pt)<TMath::Abs(y_neg005_005v2_pt.at(1)-y_neg005_005v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, (y_neg005_005v2.at(i)+y_015_025v2.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.15<y-1 && y-1<0.25)
		{
			for(int i=0; i<y_015_025v2_pt.size(); i++)
			{
				if(TMath::Abs(y_015_025v2_pt.at(i)-pt)<TMath::Abs(y_015_025v2_pt.at(1)-y_015_025v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, y_015_025v2.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.25<y-1 && y-1<0.35)
		{
			for(int i=0; i<y_035_045v2_pt.size(); i++)
			{
				if(TMath::Abs(y_035_045v2_pt.at(i)-pt)<TMath::Abs(y_035_045v2_pt.at(1)-y_035_045v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, (y_015_025v2.at(i)+y_035_045v2.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.35<y-1 && y-1<0.45)
		{
			for(int i=0; i<y_035_045v2_pt.size(); i++)
			{
				if(TMath::Abs(y_035_045v2_pt.at(i)-pt)<TMath::Abs(y_035_045v2_pt.at(1)-y_035_045v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, y_035_045v2.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.45<y-1 && y-1<0.55)
		{
			for(int i=0; i<y_055_065v2_pt.size(); i++)
			{
				if(TMath::Abs(y_055_065v2_pt.at(i)-pt)<TMath::Abs(y_055_065v2_pt.at(1)-y_055_065v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, (y_035_045v2.at(i)+y_055_065v2.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.55<y-1 && y-1<0.65)
		{
			for(int i=0; i<y_055_065v2_pt.size(); i++)
			{
				if(TMath::Abs(y_055_065v2_pt.at(i)-pt)<TMath::Abs(y_055_065v2_pt.at(1)-y_055_065v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, y_055_065v2.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
					
			}
		}
		
		
		
		
		
		if(-0.15<y-1 && y-1<-0.05)
		{
			for(int i=0; i<y_neg025_neg015v2_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_neg015v2_pt.at(i)-pt)<TMath::Abs(y_neg025_neg015v2_pt.at(1)-y_neg025_neg015v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, (y_neg025_neg015v2.at(i)+y_neg005_005v2.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.25<y-1 && y-1<-0.15)
		{
			for(int i=0; i<y_neg025_neg015v2_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_neg015v2_pt.at(i)-pt)<TMath::Abs(y_neg025_neg015v2_pt.at(1)-y_neg025_neg015v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, y_neg025_neg015v2.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.35<y-1 && y-1<-0.25)
		{
			for(int i=0; i<y_neg025_neg015v2_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_neg015v2_pt.at(i)-pt)<TMath::Abs(y_neg025_neg015v2_pt.at(1)-y_neg025_neg015v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, (y_neg025_neg015v2.at(i)+y_neg045_neg035v2.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.45<y-1 && y-1<-0.35)
		{
			for(int i=0; i<y_neg045_neg035v2_pt.size(); i++)
			{
				if(TMath::Abs(y_neg045_neg035v2_pt.at(i)-pt)<TMath::Abs(y_neg045_neg035v2_pt.at(1)-y_neg045_neg035v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, y_neg045_neg035v2.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.55<y-1 && y-1<-0.45)
		{
			for(int i=0; i<y_neg045_neg035v2_pt.size(); i++)
			{
				if(TMath::Abs(y_neg045_neg035v2_pt.at(i)-pt)<TMath::Abs(y_neg045_neg035v2_pt.at(1)-y_neg045_neg035v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, (y_neg045_neg035v2.at(i)+y_neg065_neg055v2.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		
		if(-0.65<y-1 && y-1<-0.55)
		{
			for(int i=0; i<y_neg065_neg055v2_pt.size(); i++)
			{
				if(TMath::Abs(y_neg065_neg055v2_pt.at(i)-pt)<TMath::Abs(y_neg065_neg055v2_pt.at(1)-y_neg065_neg055v2_pt.at(2))/2)
				{
					try{
					fun->SetParameter(1, y_neg065_neg055v2.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
	}
	
	//v3
	if(v_n==3)
	{
		if(-0.05<y-1 && y-1<0.05)
		{
			for(int i=0; i<y_neg005_005v3_pt.size(); i++)
			{
				if(TMath::Abs(y_neg005_005v3_pt.at(i)-pt)<TMath::Abs(y_neg005_005v3_pt.at(1)-y_neg005_005v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, y_neg005_005v3.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.05<y-1 && y-1<0.15)
		{
			for(int i=0; i<y_neg005_005v3_pt.size(); i++)
			{
				if(TMath::Abs(y_neg005_005v3_pt.at(i)-pt)<TMath::Abs(y_neg005_005v3_pt.at(1)-y_neg005_005v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, (y_neg005_005v3.at(i)+y_015_025v3.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.15<y-1 && y-1<0.25)
		{
			for(int i=0; i<y_015_025v3_pt.size(); i++)
			{
				if(TMath::Abs(y_015_025v3_pt.at(i)-pt)<TMath::Abs(y_015_025v3_pt.at(1)-y_015_025v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, y_015_025v3.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.25<y-1 && y-1<0.35)
		{
			for(int i=0; i<y_035_045v3_pt.size(); i++)
			{
				if(TMath::Abs(y_035_045v3_pt.at(i)-pt)<TMath::Abs(y_035_045v3_pt.at(1)-y_035_045v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, (y_015_025v3.at(i)+y_035_045v3.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.35<y-1 && y-1<0.45)
		{
			for(int i=0; i<y_035_045v3_pt.size(); i++)
			{
				if(TMath::Abs(y_035_045v3_pt.at(i)-pt)<TMath::Abs(y_035_045v3_pt.at(1)-y_035_045v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, y_035_045v3.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.45<y-1 && y-1<0.55)
		{
			for(int i=0; i<y_055_065v3_pt.size(); i++)
			{
				if(TMath::Abs(y_055_065v3_pt.at(i)-pt)<TMath::Abs(y_055_065v3_pt.at(1)-y_055_065v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, (y_035_045v3.at(i)+y_055_065v3.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.55<y-1 && y-1<0.65)
		{
			for(int i=0; i<y_055_065v3_pt.size(); i++)
			{
				if(TMath::Abs(y_055_065v3_pt.at(i)-pt)<TMath::Abs(y_055_065v3_pt.at(1)-y_055_065v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, y_055_065v3.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		
		
		
		if(-0.15<y-1 && y-1<-0.05)
		{
			for(int i=0; i<y_neg025_neg015v3_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_neg015v3_pt.at(i)-pt)<TMath::Abs(y_neg025_neg015v3_pt.at(1)-y_neg025_neg015v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, (y_neg025_neg015v3.at(i)+y_neg005_005v3.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.25<y-1 && y-1<-0.15)
		{
			for(int i=0; i<y_neg025_neg015v3_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_neg015v3_pt.at(i)-pt)<TMath::Abs(y_neg025_neg015v3_pt.at(1)-y_neg025_neg015v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, y_neg025_neg015v3.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.35<y-1 && y-1<-0.25)
		{
			for(int i=0; i<y_neg025_neg015v3_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_neg015v3_pt.at(i)-pt)<TMath::Abs(y_neg025_neg015v3_pt.at(1)-y_neg025_neg015v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, (y_neg025_neg015v3.at(i)+y_neg045_neg035v3.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.45<y-1 && y-1<-0.35)
		{
			for(int i=0; i<y_neg045_neg035v3_pt.size(); i++)
			{
				if(TMath::Abs(y_neg045_neg035v3_pt.at(i)-pt)<TMath::Abs(y_neg045_neg035v3_pt.at(1)-y_neg045_neg035v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, y_neg045_neg035v3.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.55<y-1 && y-1<-0.45)
		{
			for(int i=0; i<y_neg045_neg035v3_pt.size(); i++)
			{
				if(TMath::Abs(y_neg045_neg035v3_pt.at(i)-pt)<TMath::Abs(y_neg045_neg035v3_pt.at(1)-y_neg045_neg035v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, (y_neg045_neg035v3.at(i)+y_neg065_neg055v3.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		
		if(-0.65<y-1 && y-1<-0.55)
		{
			for(int i=0; i<y_neg065_neg055v3_pt.size(); i++)
			{
				if(TMath::Abs(y_neg065_neg055v3_pt.at(i)-pt)<TMath::Abs(y_neg065_neg055v3_pt.at(1)-y_neg065_neg055v3_pt.at(2))/2)
				{
					try{
					fun->SetParameter(2, y_neg065_neg055v3.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
	}
	
	//v4
	if(v_n==4)
	{
		if(-0.05<y-1 && y-1<0.05)
		{
			for(int i=0; i<y_neg005_005v4_pt.size(); i++)
			{
				if(TMath::Abs(y_neg005_005v4_pt.at(i)-pt)<TMath::Abs(y_neg005_005v4_pt.at(1)-y_neg005_005v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, y_neg005_005v4.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.05<y-1 && y-1<0.15)
		{
			for(int i=0; i<y_neg005_005v4_pt.size(); i++)
			{
				if(TMath::Abs(y_neg005_005v4_pt.at(i)-pt)<TMath::Abs(y_neg005_005v4_pt.at(1)-y_neg005_005v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, (y_neg005_005v4.at(i)+y_015_025v4.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.15<y-1 && y-1<0.25)
		{
			for(int i=0; i<y_015_025v4_pt.size(); i++)
			{
				if(TMath::Abs(y_015_025v4_pt.at(i)-pt)<TMath::Abs(y_015_025v4_pt.at(1)-y_015_025v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, y_015_025v4.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.25<y-1 && y-1<0.35)
		{
			for(int i=0; i<y_035_045v4_pt.size(); i++)
			{
				if(TMath::Abs(y_035_045v4_pt.at(i)-pt)<TMath::Abs(y_035_045v4_pt.at(1)-y_035_045v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, (y_015_025v4.at(i)+y_035_045v4.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.35<y-1 && y-1<0.45)
		{
			for(int i=0; i<y_035_045v4_pt.size(); i++)
			{
				if(TMath::Abs(y_035_045v4_pt.at(i)-pt)<TMath::Abs(y_035_045v4_pt.at(1)-y_035_045v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, y_035_045v4.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.45<y-1 && y-1<0.55)
		{
			for(int i=0; i<y_055_065v4_pt.size(); i++)
			{
				if(TMath::Abs(y_055_065v4_pt.at(i)-pt)<TMath::Abs(y_055_065v4_pt.at(1)-y_055_065v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, (y_035_045v4.at(i)+y_055_065v4.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(0.55<y-1 && y-1<0.65)
		{
			for(int i=0; i<y_055_065v4_pt.size(); i++)
			{
				if(TMath::Abs(y_055_065v4_pt.at(i)-pt)<TMath::Abs(y_055_065v4_pt.at(1)-y_055_065v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, y_055_065v4.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		
		
		
		
		if(-0.15<y-1 && y-1<-0.05)
		{
			for(int i=0; i<y_neg025_neg015v4_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_neg015v4_pt.at(i)-pt)<TMath::Abs(y_neg025_neg015v4_pt.at(1)-y_neg025_neg015v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, (y_neg025_neg015v4.at(i)+y_neg005_005v4.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.25<y-1 && y-1<-0.15)
		{
			for(int i=0; i<y_neg025_neg015v4_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_neg015v4_pt.at(i)-pt)<TMath::Abs(y_neg025_neg015v4_pt.at(1)-y_neg025_neg015v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, y_neg025_neg015v4.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.35<y-1 && y-1<-0.25)
		{
			for(int i=0; i<y_neg025_neg015v4_pt.size(); i++)
			{
				if(TMath::Abs(y_neg025_neg015v4_pt.at(i)-pt)<TMath::Abs(y_neg025_neg015v4_pt.at(1)-y_neg025_neg015v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, (y_neg025_neg015v4.at(i)+y_neg045_neg035v4.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.45<y-1 && y-1<-0.35)
		{
			for(int i=0; i<y_neg045_neg035v4_pt.size(); i++)
			{
				if(TMath::Abs(y_neg045_neg035v4_pt.at(i)-pt)<TMath::Abs(y_neg045_neg035v4_pt.at(1)-y_neg045_neg035v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, y_neg045_neg035v4.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		if(-0.55<y-1 && y-1<-0.45)
		{
			for(int i=0; i<y_neg045_neg035v4_pt.size(); i++)
			{
				if(TMath::Abs(y_neg045_neg035v4_pt.at(i)-pt)<TMath::Abs(y_neg045_neg035v4_pt.at(1)-y_neg045_neg035v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, (y_neg045_neg035v4.at(i)+y_neg065_neg055v4.at(i))/2);
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
		
		
		if(-0.65<y-1 && y-1<-0.55)
		{
			for(int i=0; i<y_neg065_neg055v4_pt.size(); i++)
			{
				if(TMath::Abs(y_neg065_neg055v4_pt.at(i)-pt)<TMath::Abs(y_neg065_neg055v4_pt.at(1)-y_neg065_neg055v4_pt.at(2))/2)
				{
					try{
					fun->SetParameter(3, y_neg065_neg055v4.at(i));
					return 1;
					}
					catch(const std::out_of_range& e){
						return 0;
					}
				}
			}
		}
	}
	return 0;
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
/*
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

*/

	int tmp=0;
	int i=0;
	int which=1;

	  
    PrepareTables();
    Int_t shift = fMcEvent->GetNTracks();
	


	
		
	Double_t Psi = gRandom->Uniform(0, 2*TMath::Pi());
		
	while(i<fMultiplicity)
	{
		
		Double_t pt, y;
		fSpectras->GetRandom2(y, pt);
	
		//TFile *vfile = new TFile("/mnt/c/Users/bumcy/Desktop/Magisterka/HAL/decayer/vtest.root", "UPDATE");

//TH1D* v1hist = (TH1D*)vfile->Get("v1");
//TH1D* v2hist = (TH1D*)vfile->Get("v2");
//TH1D* v3hist = (TH1D*)vfile->Get("v3");
//TH1D* v4hist = (TH1D*)vfile->Get("v4");

		
		//set_curve_params(pt, 0, curve, v1_lim, v1_val);
		//v1hist->Fill(set_curve_params(pt, 0, curve, v1_lim, v1_val));
		//set_curve_params(pt, 1, curve, v2_lim, v2_val);
		///v2hist->Fill(set_curve_params(pt, 0, curve, v2_lim, v2_val));
		//set_curve_params(pt, 2, curve, v3_lim, v3_val);
		//v3hist->Fill(set_curve_params(pt, 0, curve, v3_lim, v3_val));
		//set_curve_params(pt, 3, curve, v4_lim, v4_val);
		//v4hist->Fill(set_curve_params(pt, 0, curve, v4_lim, v4_val));
		
		int bool1;
		int bool2;
		int bool3;
		int bool4;
			
		if(fMass==0.938272)
		{
			bool1 = set_flow(y, pt, 1, curve);
			bool2 = set_flow(y, pt, 2, curve);
			bool3 = set_flow(y, pt, 3, curve);
			bool4 = set_flow(y, pt, 4, curve);
		}
		if(fMass==0.13957)
		{
			bool1 = 1;
			bool2 = 1;
			bool3 = 1;
			bool4 = 1;
		}
		
	//	vfile->Close();
		
		Double_t phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
		Double_t check = gRandom->Uniform(0, curve->GetMaximum());
		
		//if(fMass==0.938272)
		//	std::cout<<curve->GetParameter(0)<<" "<<curve->GetParameter(1)<<" "<<curve->GetParameter(2)<<" "<<curve->GetParameter(3)<<std::endl;
		
	//	std::cout<<curve->GetMinimum()<<std::endl;
	
		
		//std::cout<<bool1<<" "<<bool2<<" "<<bool3<<" "<<bool4<<std::endl;
		
		//std::cout<<(check<curve->Eval(phi))<<std::endl;
		
		//if(!bool1 || !bool2 || !bool3 || !bool4)
		//	std::cout<<pt<<" "<<y-1<<" "<<bool1<<" "<<bool2<<" "<<bool3<<" "<<bool4<<std::endl;
	
		if(check<curve->Eval(phi) && bool1 && bool2 && bool3 && bool4)
		{
			
			if(fMass!=0.13957 && fMass!=0.938272)
				std::cout<<fMass<<std::endl;
			
			
			Double_t mt  = TMath::Sqrt(pt * pt + fMass * fMass);
			Double_t px  = pt * TMath::Cos(phi);
			Double_t py  = pt * TMath::Sin(phi);
			Double_t pz  = mt * TMath::SinH(y);
			
			
	
			OTF::McTrack tr;
			TLorentzVector p;
			
			p.SetXYZM(px, py, pz, fMass);
			
	//		std::cout<<"2 "<<p.Phi()<<std::endl;
			
			p.Rotate(Psi, TVector3(1,1,1));
			
	//		std::cout<<"3 "<<p.Phi()<<std::endl;
		
			//std::cout<<(p*1000).M()<<std::endl;
			
			//std::cout<<px<<std::endl;
			
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
			
			TLorentzVector vec = TLorentzVector(px,py,pz,e);
			
			rtr.SetMom(vec);
		//	std::cout<<px<<std::endl;
		//	std::cout<<(p*1000).M()<<std::endl;
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
		Hal::DecayChannel ch1(2212, -211, 1);
		Hal::Decay decay(2114);
		decay.AddDecayChannel(ch1);
		decay.Init();
		std::vector<Hal::McTrack*> tracks;
		for (int j = 0; j < 3; j++) {
		tracks.push_back(new Hal::McTrack());
		}
		Hal::McTrack Delta;
		
	int n_deltas=0;
		
	while(n_deltas<25)
	{
		//std::cout<<n_deltas<<std::endl;
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
			
			
			double e1 = TMath::Sqrt(tracks[0]->GetMomentum().Px() * tracks[0]->GetMomentum().Px() + tracks[0]->GetMomentum().Py() * tracks[0]->GetMomentum().Py() + tracks[0]->GetMomentum().Pz() * tracks[0]->GetMomentum().Pz() + tracks[0]->GetMomentum().M() * tracks[0]->GetMomentum().M());
			double e2 = TMath::Sqrt(tracks[1]->GetMomentum().Px() * tracks[1]->GetMomentum().Px()* + tracks[1]->GetMomentum().Py() * tracks[1]->GetMomentum().Py() + tracks[1]->GetMomentum().Pz() * tracks[1]->GetMomentum().Pz() + tracks[1]->GetMomentum().M() * tracks[1]->GetMomentum().M());
		
		
			TLorentzVector p1, p2;
			p1.SetXYZM(tracks[0]->GetMomentum().Px(), tracks[0]->GetMomentum().Py(), tracks[0]->GetMomentum().Pz(), e1);
		
			p2.SetXYZM(tracks[1]->GetMomentum().Px(), tracks[1]->GetMomentum().Py(), tracks[1]->GetMomentum().Pz(), e2);
			
			int bool1;
			int bool2;
			int bool3;
			int bool4;

			bool1 = set_flow(p1.Rapidity(), p1.Pt(), 1, curve);
			bool2 = set_flow(p1.Rapidity(), p1.Pt(), 2, curve);
			bool3 = set_flow(p1.Rapidity(), p1.Pt(), 3, curve);
			bool4 = set_flow(p1.Rapidity(), p1.Pt(), 4, curve);
			
			//std::cout<<p1.Rapidity()-1<<" "<<p1.Pt()<<" "<<bool1<<" "<<bool2<<" "<<bool3<<" "<<bool4<<std::endl;
			
			Double_t check = gRandom->Uniform(0, curve->GetMaximum());
			if(check<curve->Eval(p1.Phi()) && bool1 && bool2 && bool3 && bool4)
			{
				tr1.SetMomentum(p2);
				tr1.SetPdgCode(-211);
				TLorentzVector xr1(gRandom->Gaus(0, 1), gRandom->Gaus(0, 1), gRandom->Gaus(0), 0);
				tr1.SetFreezout(xr1);
				fMcEvent->AddTrack(tr1);
			
				//std::cout<<p1.M()<<" " <<p2.M()<<std::endl;
				
				rtr2.SetMom(tracks[1]->GetMomentum());
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
				
				rtr1.SetMom(tracks[0]->GetMomentum());
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
			}
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
