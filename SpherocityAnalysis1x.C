#include <iostream>
#include <fstream>
#include "TFile.h"
#include "THnSparse.h"
#include <TNtuple.h>
#include <TTree.h>
#include "TSystem.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TAxis.h"
#include <vector>
#include "TApplication.h"
#include <unistd.h>
using namespace std;

void SpherocityAnalysis1a(){
  
  double gg= 0.25;//eta width
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813; 

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");//Reading Output from Pythia.
  TFile *output=new TFile("analysis_sf.root", "recreate"); 
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 1500000, 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta<0 ) continue;
     // cout<<Eta<<endl;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
     
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
       
       
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
      
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
 
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    ////
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1f.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1b(){
  double gg=0.25;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5) continue;
      if( Eta<-gg || Eta>0) continue;
	
      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}

///////////////////////////

//////////

////////////

/*void addingtree()
{
TChain chain("tr1");
chain.Add("analysis_s01.root");
chain.Add("analysis_s0.root"); 
chain.Scan( "s01:s0");





}
void hist(){
auto c = new TCanvas();
TFile *f1 = TFile::Open("AnalysisOutput1.root");
hSpherocity->Draw();
TFile *f2 = TFile::Open("AnalysisOutput1b.root");
hSpherocity->SetLineColorAlpha(kRed, 0.35);
hSpherocity->Draw("SAME");



y); //read only this branch
   if (fChain == 0) return;
   double sumf=0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      sumf+=S_b;
      // if (Cut(ientry) < 0) continue;
   }
   cout<<sumf<<endl;
}
*/


void addingtree1()
{

TFile *f=new TFile("analysis_sb.root");
TTree *tr1fb=(TTree*)f->Get("tr1b");

tr1fb->AddFriend("tr1","analysis_sf.root");
//declaration of variables
double sumfb=0;
double sumf=0;
double sumb=0;
double sumf2=0 ;
 TLeaf *leaf_f = tr1fb->GetLeaf("s0"); 
 TLeaf *leaf_b = tr1fb->GetLeaf("s01"); 
 TBranch *brf=leaf_f->GetBranch();
 TBranch *brb=leaf_b->GetBranch();
 
int nentries= tr1fb->GetEntries();
for (int i=0;i<nentries; i++) {
tr1fb->GetEntry(i);
brf->GetEntry(i) ;
brb->GetEntry(i) ;

 double valuef = leaf_f->GetValue();
 double valueb = leaf_b->GetValue();
 sumfb+=valuef*valueb;
 sumf+=valuef;
 sumb+=valueb;
 sumf2+=pow(valuef,2);
 
 }
 double sumfb_avg= sumfb/nentries;
 double sumf_avg=sumf/nentries;
 double sumf_avg_2=pow(sumf_avg,2);
 double sumb_avg=sumb/nentries;
 double sumf_b_avg=sumf_avg*sumb_avg;
 double sumf2_avg=sumf2/nentries;
 double b= (sumfb_avg - sumf_b_avg)/(sumf2_avg - sumf_avg_2);

cout<<"Pearson correlation coefficient= "<<b<<endl;


 //sumf+=tr1b->s0;
//cout<<sumf<<endl;}
//tree->Print();
//tr1fb->Scan("s0:s01");
//tr1fb->Draw("s0*s01");
//tr1b->Draw("s0"); 
//tr1b->Draw("pow(s0,2)");

//tr1b->MakeClass("analysis_s0x");
/*int N=tr1b->GetEntries();
float sum=0;
for(int i=0; i<N; i++ )
{ tr1b->GetEntry(i);
 float Pt

 //cout<< "S02 = " <<pow(s0,2)<< endl;*/

f->Write();



}
//////////////
////

////


//2nd loop



void SpherocityAnalysis1c(){
  double gg= 0.5;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sf.root", "recreate");
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  //tree->Branch("Spherocity2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 10, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta <0) continue;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1d(){
  double gg=0.5;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5 ) continue;
      if( Eta <-gg || Eta >0) continue;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}

////
//loop 3rd
void SpherocityAnalysis1e(){
  double gg= 0.75;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sf.root", "recreate");
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  //tree->Branch("Spherocity2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 10, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta <0) continue;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1f(){
  double gg=0.75;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5 ) continue;
      if( Eta <-gg || Eta >0) continue;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}

///
//loop4th

void SpherocityAnalysis1g(){
  double gg= 1;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sf.root", "recreate");
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  //tree->Branch("Spherocity2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 10, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta <0) continue;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1h(){
  double gg=1;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5 ) continue;
      if( Eta <-gg || Eta >0) continue;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////
//

//loop5th
void SpherocityAnalysis1i(){
  double gg= 1.25;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sf.root", "recreate");
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  //tree->Branch("Spherocity2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 10, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta <0) continue;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1j(){
  double gg=1.25;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5 ) continue;
      if( Eta <-gg || Eta >0) continue;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
//////

///



////

void SpherocityAnalysis1k(){
ROOT::EnableImplicitMT();
  double gg= 1.5;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sf.root", "recreate");
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  //tree->Branch("Spherocity2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 1500000, 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta<0 ) continue;
     // cout<<Eta<<endl;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
     
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
       
       
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
      
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
 
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    ////
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1f.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1l(){
  double gg=1.5;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5) continue;
      if( Eta<-gg || Eta>0) continue;
	
      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}

///////////////////////////

//////////

////////////

/*void addingtree()
{
TChain chain("tr1");
chain.Add("analysis_s01.root");
chain.Add("analysis_s0.root"); 
chain.Scan( "s01:s0");





}
void hist(){
auto c = new TCanvas();
TFile *f1 = TFile::Open("AnalysisOutput1.root");
hSpherocity->Draw();
TFile *f2 = TFile::Open("AnalysisOutput1b.root");
hSpherocity->SetLineColorAlpha(kRed, 0.35);
hSpherocity->Draw("SAME");



y); //read only this branch
   if (fChain == 0) return;
   double sumf=0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      sumf+=S_b;
      // if (Cut(ientry) < 0) continue;
   }
   cout<<sumf<<endl;
}
*/



///
///






void SpherocityAnalysis1m(){

  double gg= 1.75;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sf.root", "recreate");
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  //tree->Branch("Spherocity2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 1500000, 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta<0 ) continue;
     // cout<<Eta<<endl;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
     
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
       
       
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
      
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
 
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    ////
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1f.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1n(){
  double gg=1.75;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5) continue;
      if( Eta<-gg || Eta>0) continue;
	
      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}

///////////////////////////

//////////

////////////

/*void addingtree()
{
TChain chain("tr1");
chain.Add("analysis_s01.root");
chain.Add("analysis_s0.root"); 
chain.Scan( "s01:s0");





}
void hist(){
auto c = new TCanvas();
TFile *f1 = TFile::Open("AnalysisOutput1.root");
hSpherocity->Draw();
TFile *f2 = TFile::Open("AnalysisOutput1b.root");
hSpherocity->SetLineColorAlpha(kRed, 0.35);
hSpherocity->Draw("SAME");



y); //read only this branch
   if (fChain == 0) return;
   double sumf=0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      sumf+=S_b;
      // if (Cut(ientry) < 0) continue;
   }
   cout<<sumf<<endl;
}
*/





//

void SpherocityAnalysis1o(){

  double gg= 2;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sf.root", "recreate");
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  //tree->Branch("Spherocity2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 1500000, 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta<0 ) continue;
     // cout<<Eta<<endl;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
     
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
       
       
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
      
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
 
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    ////
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1f.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1p(){
  double gg=2;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5) continue;
      if( Eta<-gg || Eta>0) continue;
	
      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}

///////////////////////////

//////////

////////////

/*void addingtree()
{
TChain chain("tr1");
chain.Add("analysis_s01.root");
chain.Add("analysis_s0.root"); 
chain.Scan( "s01:s0");





}
void hist(){
auto c = new TCanvas();
TFile *f1 = TFile::Open("AnalysisOutput1.root");
hSpherocity->Draw();
TFile *f2 = TFile::Open("AnalysisOutput1b.root");
hSpherocity->SetLineColorAlpha(kRed, 0.35);
hSpherocity->Draw("SAME");



y); //read only this branch
   if (fChain == 0) return;
   double sumf=0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      sumf+=S_b;
      // if (Cut(ientry) < 0) continue;
   }
   cout<<sumf<<endl;
}
*/




void SpherocityAnalysis1q(){

  double gg= 2.25;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sf.root", "recreate");
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  //tree->Branch("Spherocity2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 1500000, 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta<0 ) continue;
     // cout<<Eta<<endl;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
     
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
       
       
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
      
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
 
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    ////
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1f.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1r(){
  double gg=2.25;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5) continue;
      if( Eta<-gg || Eta>0) continue;
	
      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}

///////////////////////////

//////////

////////////

/*void addingtree()
{
TChain chain("tr1");
chain.Add("analysis_s01.root");
chain.Add("analysis_s0.root"); 
chain.Scan( "s01:s0");





}
void hist(){
auto c = new TCanvas();
TFile *f1 = TFile::Open("AnalysisOutput1.root");
hSpherocity->Draw();
TFile *f2 = TFile::Open("AnalysisOutput1b.root");
hSpherocity->SetLineColorAlpha(kRed, 0.35);
hSpherocity->Draw("SAME");



y); //read only this branch
   if (fChain == 0) return;
   double sumf=0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      sumf+=S_b;
      // if (Cut(ientry) < 0) continue;
   }
   cout<<sumf<<endl;
}
*/




void SpherocityAnalysis1s(){

  double gg= 2.5;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s0;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sf.root", "recreate");
  TTree *tree = new TTree("tr1","Spherocity_tree");
  tree->Branch("S_f",&s0,"s0/D");
  //tree->Branch("Spherocity2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 1500000., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 1500000, 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if( Pt<0.5 ) continue;
      if( Eta >gg || Eta<0 ) continue;
     // cout<<Eta<<endl;

      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
     
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	

	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
       
       
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
      
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s0=Spherocity;
 
    s02=pow(Spherocity,2);
    tree->Fill();
   
  //  cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    ////
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1f.root", "recreate");
 // hSphericity->Write();
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}
////////////////////////////////////
////////////////////////
//////////////
















void SpherocityAnalysis1t(){
  double gg=2.5;
  const Int_t kmaxtrack = 8000;
  const Int_t kCbin = 90;
  const Int_t kptbin = 12;
  const Int_t ptdim = 24; //12*2
  const Int_t kEta = 15;
  const Int_t kMinMult = 5; //minimum multiplicty ---
  Float_t pmass = 0.9382720813;

  Int_t total_events = 0.;
  double s01;
  double s02;
  TFile *file =  TFile::Open("tut4.root");
  TFile *output=new TFile("analysis_sb.root", "recreate");
  TTree *tree1 = new TTree("tr1b","Spherocity_tree");
  tree1->Branch("S_b",&s01,"s01/D");
  //tree->Branch("S2",&s02,"s02/D");
  

  Int_t mult;
  Int_t nTrack;
  //Float_t sp;
  //Int_t charge[kmaxtrack];
  
  Float_t pt[kmaxtrack];
  Float_t phi[kmaxtrack];
  Float_t eta[kmaxtrack];
  //Int_t motherID[kmaxtrack];
  
  TTree *t = (TTree*)file->Get("tr");
  
  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("nTrack", &nTrack);
  //t->SetBranchAddress("charge", charge);
  t->SetBranchAddress("pt", pt);
  t->SetBranchAddress("eta", eta);
  t->SetBranchAddress("phi", phi);
 // t->SetBranchAddress("sp",&sp);
  //->Scan();
  //t->SetBranchAddress("motherID", motherID);

  TH1D *hSpherocity = new TH1D("hSpherocity", "Spherocity", 100, 0., 1.);
//  hSpherocity->GetXaxis()->SetTitle("S_{0}");
   TH1D *hSphericity = new TH1D("hSphericity", "Sphericity", 100, 1e9, 2e9);
  hSphericity->GetXaxis()->SetTitle("S_{0}");

  TProfile *hMultSAll   = new TProfile("hMultSAll","MinBias",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSJetty = new TProfile("hMultSJetty","Jetty",   150, 0., 100., 0.0, 2.5);
  TProfile *hMultSIso   = new TProfile("hMultSIso","Isotropic", 150, 0., 100., 0.0, 2.5);
  TH2D *h2NchPt = new TH2D("h2NchPt","h2NchPt", 150, 0., 150., 30, 0.0, 1.5);
  TH1D *histoPt = new TH1D("hPt","Pt", 50, 0., 10.);
  
  Int_t nevents = (Int_t) t->GetEntries();
  
  cout << "Events " << nevents << endl;
  
  Int_t sumEve = 0;
  
  for( Int_t iev = 0; iev < nevents; iev++){
    ///cout<< "iev\n" << iev;
    t->GetEntry(iev);
     //Float_t sphi=sp;
    Int_t ChMult = mult;
   // cout<< "mult\n" << mult;
    Int_t Ntracks = nTrack;
    
    vector <Double_t> vecPx;
    vector <Double_t> vecPy;
    vector <Double_t> SphCrossProd;
    
    Double_t SumTrack = 0., SumPt = 0.,  AvPt = 0.;
    
    for(Int_t itrk = 0; itrk < Ntracks; itrk++){
      
      Float_t Pt = pt[itrk];
      Float_t Eta = eta[itrk];
      Float_t Phi = phi[itrk];
      //Int_t Charge = charge[itrk];
      
      Double_t Px = Pt*TMath::Cos(Phi);
      Double_t Py = Pt*TMath::Sin(Phi);
      
      //Track cuts-----
      if(Pt<0.5) continue;
      if( Eta<-gg || Eta>0) continue;
	
      vecPx.push_back(Px);
      vecPy.push_back(Py);
      
      SumPt += Pt;
      SumTrack += 1.;
      
      histoPt->Fill( Pt );
      
    }//Track loop--itrack
    
    if( SumTrack < kMinMult ) continue;
    
    AvPt = SumPt/SumTrack;
    h2NchPt->Fill( SumTrack, AvPt );
    
    for(Int_t itrk = 0; itrk < SumTrack; itrk++){
      
      TVector3 vPTi;
      vPTi.SetXYZ( vecPx[itrk], vecPy[itrk], 0 );
      
      Double_t SumCrosProd = 0.;
      for(Int_t jtrk = 0; jtrk < SumTrack; jtrk++){
	
	TVector3 vPTj;
	vPTj.SetXYZ( vecPx[jtrk], vecPy[jtrk], 0. );
	TVector3 vecCross = vPTj.Cross( vPTi.Unit() );
	SumCrosProd += vecCross.Mag(); //pt(j)Xnhat(i)

      }//jtrk---

      Double_t RatioSquared = TMath::Power((SumCrosProd/SumPt), 2);
      
      SphCrossProd.push_back( RatioSquared );

    }//itrk------
    
    Double_t *SpheroArray;
    Int_t track_size = SphCrossProd.size();
    if( SumTrack != track_size ) cout <<"Something is wrong here " << endl;
    
    SpheroArray = new Double_t[track_size];
    
    for(Int_t ii = 0; ii < track_size; ii++) SpheroArray[ii] = SphCrossProd[ii];
    
    Double_t minSphero = TMath::MinElement(track_size, SpheroArray); 
    
    Double_t Spherocity = (TMath::Pi()*TMath::Pi()/4.)*minSphero;
    s01=Spherocity;
    s02=pow(Spherocity,2);
    tree1->Fill();
   
    //cout<<  "Tracks= "<<Ntracks<<" Spherocity ="<<Spherocity<<endl;
    

    //cout << " Spherocity " << Spherocity << " dNchdEta " << ChMult << " AvPt =" << AvPt << endl;

    //Clear the vectors, array here----
    vecPx.clear();
    vecPy.clear();
    SphCrossProd.clear();
    delete [] SpheroArray;
    
    hSpherocity->Fill( Spherocity );
    //hSphericity->Fill(sphi);
    hMultSAll->Fill( ChMult, AvPt, 1 );
    if( Spherocity < 0.1 ) hMultSJetty->Fill( ChMult, AvPt, 1 );
    if( Spherocity > 0.9 ) hMultSIso->Fill( ChMult, AvPt, 1 );
    
    
    
    sumEve += 1;
    
  }//event loop--
  
  cout << "Total event number " << sumEve << endl;
  histoPt->Scale(1./sumEve);
  //tree->Scan();
  TFile *fout = new TFile("AnalysisOutput1b.root", "recreate");
 // hSphericity->Write();
 
  hSpherocity->Write();
  hMultSAll->Write();
  hMultSJetty->Write();
  hMultSIso->Write();
  h2NchPt->Write();
  histoPt->Write();
  output->Write();
  output->Close();					   
  fout->Write();
  fout->Close();
  file->Close();
}

/////






void SpherocityAnalysis1x()
{//loop 1
SpherocityAnalysis1a();// forward 
SpherocityAnalysis1b();//backward 
addingtree1();	//calculating Pearson correlation coefficient for giving eta width
//loop 2
SpherocityAnalysis1c();// forward 
SpherocityAnalysis1d();//backward 
addingtree1();	//calculating Pearson correlation coefficient for giving eta width
//loop 3
SpherocityAnalysis1e();// forward 
SpherocityAnalysis1f();//backward 
addingtree1();//calculating Pearson correlation coefficient for giving eta width
//loop 4
SpherocityAnalysis1g();// forward 
SpherocityAnalysis1h();//backward 
addingtree1();//calculating Pearson correlation coefficient for giving eta width
//loop 5
SpherocityAnalysis1i();// forward 
SpherocityAnalysis1j();//backward 
addingtree1();//calculating Pearson correlation coefficient for giving eta width
//loop 6
SpherocityAnalysis1k();// forward 
SpherocityAnalysis1l();//backward 
addingtree1();//calculating Pearson correlation coefficient for giving eta width
//loop 7
SpherocityAnalysis1m();// forward 
SpherocityAnalysis1n();//backward 
addingtree1();//calculating Pearson correlation coefficient for giving eta width
//loop 8
SpherocityAnalysis1o();// forward 
SpherocityAnalysis1p();//backward 
addingtree1();//calculating Pearson correlation coefficient for giving eta width
//loop 9
SpherocityAnalysis1q();// forward 
SpherocityAnalysis1r();//backward 
addingtree1();//calculating Pearson correlation coefficient for giving eta width
//loop 10
SpherocityAnalysis1s();// forward 
SpherocityAnalysis1t();//backward 
addingtree1();//calculating Pearson correlation coefficient for giving eta width
}








//---------------functions-----------------------------------------------
