#include "gsort.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <TVector3.h>
#include <Math/Interpolator.h>

typedef ROOT::Math::Interpolator interp;

using namespace std;

interp *zner;
interp *znre;
interp *caer;
interp *care;

#include "masslcut.cpp"
#include "massrcut.cpp"

#define P0 3886.37 // updated 5/10/2018
//#define P0 3406.37 // updated 5/10/2018
//#define P0 3841.2
#define Mtot 118.0
#define beamE 168.8 // updated 5/10/2018 beamE=175,target thickness 1mg/cm2
//#define beamE 165.0
#define beamA 48.0
#define targetA 70.0
#define QNum 7
#define GTPosOffset 19.4
#define atomic_mass 931.4940954 // unit MeV/c2
#define v_c 299792458000 // mm/s

bool BuildTable(const char*filename,interp* &er,interp* &re){
  ifstream in(filename);
  if(!in.is_open()) {
    cout<<"Cannot open file: "<<filename<<endl;
    return false;
  }
  Double_t energy[200],range[200];
  int cnt = 0;
  Double_t x,y,z;
  while(1){
    in>>energy[cnt]>>range[cnt];
    if(!in.good()){
      energy[cnt]=0;
      range[cnt]=0;
      break;
    }
    cnt++;
    if(cnt>=200) break;
  }
  
  er = new interp(cnt,ROOT::Math::Interpolation::kCSPLINE);
  re = new interp(cnt,ROOT::Math::Interpolation::kCSPLINE);
  er->SetData(cnt,energy,range);
  re->SetData(cnt,range,energy);
  return true;
}

Double_t updatep(Double_t pin,Double_t thick, int index=0){
  if(pin<=0) return 0;
  interp* er,*re;
  Double_t mass;
  if(index==0) {
    er = zner;
    re = znre;
    mass = 70;
  }else{
    er = caer;
    re = care;
    mass = 48;
  }
  Double_t ein = pin*pin*0.5/(mass*atomic_mass);
  Double_t range = er->Eval(ein)-thick;
  ein = re->Eval(range);
  return sqrt(2*mass*atomic_mass*ein);
}

Float_t offGT[QNum*4],gainGT[QNum*4];

template <class T,int m>
void readp(const char* fn, T (& off)[m], T (& gain)[m]){
  ifstream in(fn);
  if(!in.is_open()){
    cout<<"cannot open file: "<<fn<<endl;
    return;
  }
  int cnt=0;
  float x,y;
  while(1){
    in>>x>>y;
    if(!in.good()) break;
    off[cnt]=x;gain[cnt]=y;
    cnt++;
    if(cnt>=m) break;
  }
}

void gsort::Loop(TTree *opt)
{
//   In a ROOT session, you can do:
//      Root > .L gsort.C
//      Root > gsort t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   if (opt==0) return;
   BranchOpt(opt);

   Long64_t nentries = fChain->GetEntriesFast();


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(!GetMass()) continue;
      gammadc();
      tkgammadc();
      opt->Fill();
      // if (Cut(ientry) < 0) continue;
   }
}

gsort::gsort(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rootfile/run0009.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("rootfile/run0009.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
   masslcut=load_masslcut();
   massrcut=load_massrcut();
   readp("GT_112014_2.cal",offGT,gainGT);
   BuildTable("calcium_range.txt",caer,care);
   BuildTable("zinc_range.txt",zner,znre);
}

gsort::~gsort()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gsort::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gsort::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gsort::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("cts", cts, &b_cts);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("dT", dT, &b_dT);
   fChain->SetBranchAddress("dL", dL, &b_dL);
   fChain->SetBranchAddress("dR", dR, &b_dR);
   fChain->SetBranchAddress("phiL", phiL, &b_phiL);
   fChain->SetBranchAddress("phiR", phiR, &b_phiR);
   fChain->SetBranchAddress("fphiL", fphiL, &b_fphiL);
   fChain->SetBranchAddress("fphiR", fphiR, &b_fphiR);
   fChain->SetBranchAddress("thetaL", thetaL, &b_thetaL);
   fChain->SetBranchAddress("thetaR", thetaR, &b_thetaR);
   fChain->SetBranchAddress("fthetaL", fthetaL, &b_fthetaL);
   fChain->SetBranchAddress("fthetaR", fthetaR, &b_fthetaR);
   fChain->SetBranchAddress("ng", &ng, &b_ng);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("e", e, &b_e);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("t", t, &b_t);
   fChain->SetBranchAddress("offt", offt, &b_offt);
   fChain->SetBranchAddress("t0", t0, &b_t0);
   fChain->SetBranchAddress("dtpg", dtpg, &b_dtpg);
   fChain->SetBranchAddress("cc_id", cc_id, &b_cc_id);
   fChain->SetBranchAddress("cry_id", cry_id, &b_cry_id);
   fChain->SetBranchAddress("ntg", &ntg, &b_ntg);
   fChain->SetBranchAddress("tsTK", &tsTK, &b_tsTK);
   fChain->SetBranchAddress("pad", pad, &b_pad);
   fChain->SetBranchAddress("tracked", tracked, &b_tracked);
   fChain->SetBranchAddress("esum", esum, &b_esum);
   fChain->SetBranchAddress("ndet", ndet, &b_ndet);
   fChain->SetBranchAddress("gtkts", gtkts, &b_gtkts);
   fChain->SetBranchAddress("fom", fom, &b_fom);
   fChain->SetBranchAddress("x0", x0, &b_x0);
   fChain->SetBranchAddress("y0", y0, &b_y0);
   fChain->SetBranchAddress("z0", z0, &b_z0);
   fChain->SetBranchAddress("e0", e0, &b_e0);
   fChain->SetBranchAddress("x1", x1, &b_x1);
   fChain->SetBranchAddress("y1", y1, &b_y1);
   fChain->SetBranchAddress("z1", z1, &b_z1);
   fChain->SetBranchAddress("e1", e1, &b_e1);
   fChain->SetBranchAddress("fhcrID", fhcrID, &b_fhcrID);
   fChain->SetBranchAddress("gid", gid, &b_gid);
   fChain->SetBranchAddress("gofft", gofft, &b_gofft);
   fChain->SetBranchAddress("dtptg", dtptg, &b_dtptg);
   Notify();
}

Bool_t gsort::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gsort::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gsort::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void gsort::BranchOpt(TTree *opt){
  assert(opt!=0);

  opt->Branch("MassT",&MassT,"MassT/D");
  opt->Branch("MassP",&MassP,"MassP/D");
  opt->Branch("bt",&bt,"bt/D");
  opt->Branch("bp",&bp,"bp/D");
  opt->Branch("theT",&theT,"theT/D");
  opt->Branch("phiT",&phiT,"phiT/D");
  opt->Branch("theP",&theP,"theP/D");
  opt->Branch("phiP",&phiP,"phiP/D");
  opt->Branch("qval",&qval,"qval/D");
  opt->Branch("qval2",&qval2,"qval2/D");
  opt->Branch("np",&np,"np/I");
  opt->Branch("cts",&cts,"cts[np]/L");

  opt->Branch("ng",&ng,"ng/I");
  opt->Branch("t",&t,"t[ng]/D");
  opt->Branch("t0",&t0,"t0[ng]/D");
  opt->Branch("x",&x,"x[ng]/F");
  opt->Branch("y",&y,"y[ng]/F");
  opt->Branch("z",&z,"z[ng]/F");
  opt->Branch("ge",&ge,"ge[ng]/D");
  opt->Branch("get",&get,"get[ng]/D");
  opt->Branch("gep",&gep,"gep[ng]/D");
  opt->Branch("gthe",&gthe,"gthe[ng]/D");
  opt->Branch("gphi",&gphi,"gphi[ng]/D");
  opt->Branch("gthet",&gthet,"gthet[ng]/D");
  opt->Branch("gthep",&gthep,"gthep[ng]/D");
  opt->Branch("dtpg",&dtpg,"dtpg[ng]/D");

  opt->Branch("ntg",&ntg,"ntg/I");
  opt->Branch("x0",&x0,"x0[ntg]/F");
  opt->Branch("y0",&y0,"y0[ntg]/F");
  opt->Branch("z0",&z0,"z0[ntg]/F");
  opt->Branch("e0",&e0,"e0[ntg]/F");
  opt->Branch("x1",&x1,"x1[ntg]/F");
  opt->Branch("y1",&y1,"y1[ntg]/F");
  opt->Branch("z1",&z1,"z1[ntg]/F");
  opt->Branch("e1",&e1,"e1[ntg]/F");
  opt->Branch("ndet",&ndet,"ndet[ntg]/I");
  opt->Branch("fom",&fom,"fom[ntg]/F");
  opt->Branch("gtkts",&gtkts,"gtkts[ntg]/L");
  opt->Branch("tge",&tge,"tge[ntg]/D");
  opt->Branch("tget",&tget,"tget[ntg]/D");
  opt->Branch("tgep",&tgep,"tgep[ntg]/D");
  opt->Branch("tgthe",&tgthe,"tgthe[ntg]/D");
  opt->Branch("tgphi",&tgphi,"tgphi[ntg]/D");
  opt->Branch("tgthet",&tgthet,"tgthet[ntg]/D");
  opt->Branch("tgthep",&tgthep,"tgthep[ntg]/D");
  opt->Branch("dtptg",&dtptg,"dtptg[ntg]/D");

}

bool gsort::GetMass(){
  float theta1=0.0,theta2=0.0;
  float phi1=0,phi2=0;
  float L1=0.0,L2=0.0;
  float dt=0;
  float f1,f2,p1,p2;
  if(masslcut->IsInside(fthetaL[0],dT[0])){
    theta1 = fthetaL[0];
    theta2 = fthetaR[0];
    phi1=fphiL[0];
    phi2=fphiR[0];
    L1 = dL[0];
    L2 = dR[0];
//    dt = (dT[0]-2.5)*0.1; // Shaofei's code
    dt = (dT[0]-0.0)*0.1;
  }else if(massrcut->IsInside(fthetaR[0],dT[0])){
    theta1 = fthetaR[0];
    theta2 = fthetaL[0];
    phi1 = fphiR[0];
    phi2 = fphiL[0];
    L1 = dR[0];
    L2 = dL[0];
//    dt = -1.*(dT[0]-2.5)*0.1;
    dt = -1.*(dT[0]-0.0)*0.1;
  }else return false;
//  f1 = 1.1094 - 0.12008/TMath::Cos(theta1) - 0.0096272/TMath::Power(TMath::Cos(theta1),2);
//  f1*=1.06;
//  f2 = 1.96556 - 4.3451*TMath::Cos(theta2) + 4.2348*TMath::Power(TMath::Cos(theta2),2);
  p1 = P0 * TMath::Sin(theta2)/TMath::Sin(theta1+theta2);//*f1;
  p2 = P0 * TMath::Sin(theta1)/TMath::Sin(theta1+theta2);//*f2;
  
//  MassT = (-0.032206*dt + Mtot*L1/p1)/(L1/p1+L2/p2);
//  MassP = (0.032206*dt + Mtot*L2/p2)/(L1/p1+L2/p2);

  L1 = L1 * 10;
  L2 = L2 * 10;
//  p1 = updatep(p1,0.70028011/cos(theta1),1);
//  p2 = updatep(p2,0.70028011/cos(theta2),0);
  p1 = updatep(p1,0.45028011/cos(theta1),1);
  p2 = updatep(p2,0.45028011/cos(theta2),0);
  double temp = (L1/p2/v_c+L2/p1/v_c)*atomic_mass;
  MassP = (dt*1e-9+L2/p2*Mtot*atomic_mass/v_c)/temp;
  MassT = (-dt*1e-9+L1/p1*Mtot*atomic_mass/v_c)/temp;

  bt = p2/MassT/atomic_mass;
  bp = p1/MassP/atomic_mass;
//  bt = p2/targetA/atomic_mass;
//  bp = p1/beamA/atomic_mass;

  Double_t et,et2,ep;
  ep = TMath::Power(TMath::Sin(theta2)/TMath::Sin(theta1+theta2),2)*beamE;
//  bp = TMath::Sqrt(ep/beamA)*0.046;
//  bp *= f1;
  et = TMath::Power(TMath::Sin(theta1)/TMath::Sin(theta1+theta2),2)*beamE*beamA/targetA;
  et2 = TMath::Power(TMath::Sin(theta1)/TMath::Sin(theta1+theta2),2)*beamE*48./70.;
//  bt = sqrt(et/targetA)*0.046;
//  bt *=f2;

  theP = theta1;
  phiP = phi1;
  theT = theta2;
  phiT = phi2;

  qval = ep+et-beamE;
  qval2 = ep+et2-beamE;
  return true;
}

void gsort::Clear(){
  MassT=0;
  MassP=0;
  bt=0;
  bp=0;
  qval=0;
  memset(tget,0,sizeof(tget));
}

void gsort::gammadc(){
  if(ng<=0) return;
  for(int i=0;i<ng;i++){
    int gid = (cc_id[i]-1)*4+cry_id[i];
    assert(gid<28);
    ge[i]=e[i]*gainGT[gid]+offGT[gid];
    TVector3 gp(x[i],y[i],z[i]);
    TVector3 pp,pt;
    pt.SetMagThetaPhi(1,theT,phiT);
    pp.SetMagThetaPhi(1,theP,phiP);
    gthe[i]=gp.Angle(TVector3(0,0,1));
    gphi[i]=gp.Phi();
    gthet[i]=gp.Angle(pt);
    gthep[i]=gp.Angle(pp);
    get[i]=ge[i]*(1-bt*TMath::Cos(gthet[i]))/TMath::Sqrt(1-bt*bt);
    gep[i]=ge[i]*(1-bp*TMath::Cos(gthep[i]))/TMath::Sqrt(1-bp*bp);
  }
}

void gsort::tkgammadc(){
  if(ntg<=0) return;
  TVector3 pp,pt;
  pt.SetMagThetaPhi(1,theT,phiT);
  pp.SetMagThetaPhi(1,theP,phiP);
  for(int i=0;i<ntg;i++){
    tge[i]=esum[i];
    TVector3 gp1(x0[i],y0[i]-GTPosOffset,z0[i]);
    TVector3 gp2(x1[i]-x0[i],y1[i]-y0[i],z1[i]-z0[i]);
    tgthe[i]=gp1.Angle(TVector3(0,0,1));
    tgphi[i]=gp1.Phi();
    tgthet[i]=gp1.Angle(pt);
    tgthep[i]=gp1.Angle(pp);
    tget[i]=tge[i]*(1-bt*TMath::Cos(tgthet[i]))/TMath::Sqrt(1-bt*bt);
    tgep[i]=tge[i]*(1-bp*TMath::Cos(tgthep[i]))/TMath::Sqrt(1-bp*bp);
  }
}
