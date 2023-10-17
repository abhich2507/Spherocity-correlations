#include <iostream>
//#include <vector>
#include "Pythia8/Pythia.h"
//#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
using namespace Pythia8;

int main(){
	

	
	TFile *output=new TFile("tut4.root", "recreate");
	TTree *tree = new TTree("tr","tr");
	int ntrackmax=10000;
	float phi[ntrackmax]={0};
	float eta[ntrackmax]={0};
	float pt[ntrackmax]={0};
	float sp;
	int mult=0;

	int nTrack;
	//branches
	
	tree->Branch("nTrack",&nTrack,"nTrack/I");
	tree->Branch("phi",phi, "phi[nTrack]/F");
	tree->Branch("eta",eta, "eta[nTrack]/F");
	tree->Branch("pt",pt, "pt[nTrack]/F");
	tree->Branch("mult",&mult, "mult/I");
	tree->Branch("sp",&sp, "sp/F");
	

	int nevents =10000;
	 
	Sphericity sph;
	Pythia8::Pythia pythia;
	pythia.readString("Beams:idA= 2212");
	pythia.readString("Beams:idB= 2212");
	pythia.readString("Beams:eCM= 13000");
	pythia.readString("SoftQCD:all=on");//for min-bias
	pythia.readString("Tune:pp=14");
	pythia.readString("PartonLevel:MPI=on");// for MPI
	//pythia.readString("PhaseSpace:pTHatMin =5");
	//pythia.readString("SoftQCD:nonDiffractive = on"); 
	


	 
	pythia.init();
	for (int i=0; i<nevents; i++)
	
	
	{	
	//nTrack[i]=(pythia.event.size());
	 
	if(!pythia.next()) continue;
		
		 if (sph.analyze( pythia.event )) {
      			if (nevents < 3) sph.list();
    			  sp=sph.sphericity();  		
		}
		 nTrack =pythia.event.size();
		std::cout<< "Event:" << i << std::endl;
		std::cout<< "Event size: " <<  nTrack << std::endl;
    		for (int k = 0; k< pythia.event.size(); ++k)
    		
    		{
    		phi[k]= pythia.event[k].phi();
    		eta[k]= pythia.event[k].eta();
    		pt[k]= pythia.event[k].pT();
      		if ( pythia.event[k].isFinal() && pythia.event[k].isCharged() )
      		mult++;  
      		
      		
		
		
		 

		
		}
		cout<<mult<<endl;
		tree->Fill();
	    }
	    
	pythia.stat();
	output->Write();
	
	output->Close();
	

	return 0;
}
	
