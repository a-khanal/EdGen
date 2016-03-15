#define analysis_7_cxx
// The class definition in analysis_7.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("analysis_7.C")
// Root > T->Process("analysis_7.C","some options")
// Root > T->Process("analysis_7.C+")
//

#include "analysis_7.h"
#include <TH2.h>
#include <TStyle.h>


void analysis_7::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void analysis_7::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   h1_phi = new TH1F("h1_phi","#phi #pi^{0} distribution",150,-180,360);
   h1_costheta = new TH1F("h1_costheta","cos(#theta) #pi^{0} distribution (W rest frame)",100,-1.,1.);
   h1_phi_cross = new TH1F("h1_phi_cross","#phi #pi^{0} distribution (CROSS SECTION)",150,-180,360);
   h1_costheta_cross = new TH1F("h1_costheta_cross","cos(#theta) #pi^{0} distribution (W rest frame) (CROSS SECTION) /100 ",100,-1.,1.);
   fOutput->Add(h1_phi);
   fOutput->Add(h1_costheta);
   fOutput->Add(h1_phi_cross);
   fOutput->Add(h1_costheta_cross);

}

Bool_t analysis_7::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either analysis_7::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  b_Ef->GetEntry(entry);
  b_px->GetEntry(entry);
  b_py->GetEntry(entry);
  b_pz->GetEntry(entry);
  b_weight->GetEntry(entry);
  b_Ein_beam->GetEntry(entry);


  TLorentzVector p_pr(px[0],py[0],pz[0],Ef[0]);
  TLorentzVector p_pi0(px[1],py[1],pz[1],Ef[1]);

  TLorentzVector p_W = p_pr+p_pi0;
  TVector3 b_3 =  p_W.BoostVector();
  b_3 = -b_3;
  p_pi0.Boost(b_3);
  h1_costheta->Fill(p_pi0.CosTheta());
  h1_phi->Fill(p_pi0.Phi()/TMath::Pi()*180.);
  h1_costheta_cross->Fill(p_pi0.CosTheta(),weight[1]/100);
  h1_phi_cross->Fill(p_pi0.Phi()/TMath::Pi()*180.,weight[1]);

   return kTRUE;
}

void analysis_7::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void analysis_7::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  TFile file_out("analysis_7_output.root","recreate");
  TList *outlist = GetOutputList();
  
  outlist->Write();
  file_out.Write();
  file_out.Close();


}
