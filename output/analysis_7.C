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

   h1_mass_omega  = new TH1F("h1_mass_omega","Mass #omega; GeV",200,0.,2.0);
   h1_mass_omega_pi_pi  = new TH1F("h1_mass_omega_pi_pi","Invariant Mass (#omega + #pi^{+} + #pi^{-}); GeV",200,0.,2.0);
   h2_pimomega_pipomega = new TH2F("h2_pimomega_pipomega","Dalitz ;M^{2}(#pi^{-}#omega);M^{2}(#pi^{+}#omega)",200,0.,2.0,200,0.,2.0);
   h1_omega1  = new TH1F("h1_omega1","Mass #omega (fast #pi^{+} fast #pi^{-}); GeV",200,0.,2.0);
   h1_omega2  = new TH1F("h1_omega2","Mass #omega (fast #pi^{+} slow #pi^{-}); GeV",200,0.,2.0);
   h1_omega3  = new TH1F("h1_omega3","Mass #omega (slow #pi^{+} fast #pi^{-}); GeV",200,0.,2.0);
   h1_omega4  = new TH1F("h1_omega4","Mass #omega (slow #pi^{+} slow #pi^{-}); GeV",200,0.,2.0);

   fOutput->Add(h1_mass_omega);
   fOutput->Add(h1_mass_omega_pi_pi);
   fOutput->Add(h2_pimomega_pipomega);
   fOutput->Add(h1_omega1);
   fOutput->Add(h1_omega2);
   fOutput->Add(h1_omega3);
   fOutput->Add(h1_omega4);
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

  TLorentzVector p_pip(px[2],py[2],pz[2],Ef[2]);
  TLorentzVector p_pim(px[3],py[3],pz[3],Ef[3]);
  TLorentzVector p_pip2(px[5],py[5],pz[5],Ef[5]);
  TLorentzVector p_pim2(px[6],py[6],pz[6],Ef[6]);
  TLorentzVector p_pi0(px[4],py[4],pz[4],Ef[4]);
  TLorentzVector p_omega(px[1],py[1],pz[1],Ef[1]);

  TLorentzVector p_d1 = p_pim+p_omega;
  TLorentzVector p_d2 = p_pip+p_omega;
  TLorentzVector p_tot = p_pip+p_pim+p_omega;

  double beta_pip1 = p_pip.Beta();
  double beta_pim1 = p_pim.Beta();
  double beta_pip2 = p_pip2.Beta();
  double beta_pim2 = p_pim2.Beta();

  double weight_v = weight[0] * weight[4];
  

  h1_mass_omega->Fill(p_omega.M(),weight_v);
  h1_mass_omega_pi_pi->Fill(p_tot.M(),weight_v);
  h2_pimomega_pipomega->Fill(p_d1.M2(),p_d2.M2(),weight_v);

  TLorentzVector p_pip_fast,p_pip_slow,p_pim_fast,p_pim_slow;
  if (beta_pip1 >  beta_pip2) {
    p_pip_fast = p_pip;
    p_pip_slow = p_pip2;
  }
  else {
    p_pip_fast = p_pip2;
    p_pip_slow = p_pip;
  }

  if (beta_pim1 >  beta_pim2) {
    p_pim_fast = p_pim;
    p_pim_slow = p_pim2;
  }
  else {
    p_pim_fast = p_pim2;
    p_pim_slow = p_pim;
  }

  TLorentzVector p_omega1 = p_pip_fast+p_pim_fast+p_pi0;
  TLorentzVector p_omega2 = p_pip_fast+p_pim_slow+p_pi0;
  TLorentzVector p_omega3 = p_pip_slow+p_pim_fast+p_pi0;
  TLorentzVector p_omega4 = p_pip_slow+p_pim_slow+p_pi0;

  h1_omega1->Fill(p_omega1.M(),weight_v);
  h1_omega2->Fill(p_omega2.M(),weight_v);
  h1_omega3->Fill(p_omega3.M(),weight_v);
  h1_omega4->Fill(p_omega4.M(),weight_v);

  
  

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
