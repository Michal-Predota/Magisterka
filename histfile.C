#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"


using namespace std;

int histfile()
{
	TFile* output;
	TFile *input;
	
	
	output = new TFile("./files/file_pip_nodelta_noflow_3040.root", "RECREATE");
	
	
	
	vector<double> pTlim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	vector<double> ylim = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
	
	int file_counter=0;
	for(int ynum=0; ynum<8; ynum++)
	{
		for(int ptnum=0; ptnum<8; ptnum++)
		{
			char fileName [100];
			sprintf(fileName, "./output_pip_noflow_3040/superpack_0/th1_%d/histo1d.root", 17+file_counter);
			//cout<<fileName<<endl;
			input= new TFile(fileName, "read");
			TCanvas *c = (TCanvas*) input->Get("canvas");
			
			char histName[100];
			sprintf(histName, "M_{inv} %f<y<%f,  %f<pT<%f", ylim.at(ynum), ylim.at(ynum+1), pTlim.at(ptnum), pTlim.at(ptnum+1));
			cout<<histName<<endl;
			TH1D* hist = (TH1D*)c->GetPrimitive(histName);
			
			
			//cout<<hist->GetEntries()<<" "<<hist->GetBinContent(100000)<<endl;
			//hist->Rebin(10);
			//hist->Sumw2(false);
			
			//hist->Scale(1/hist->Integral(), "nosw2");
			//cout<<hist->GetEntries()<<endl;
			output->cd();
			hist->Write();
			
			file_counter++;
		}
	}
	
	output->Save();
	output->Close();
	
	
	
	
	
	
	return 0;
}