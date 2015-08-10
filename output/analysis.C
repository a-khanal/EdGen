#define analysis_cxx
// The class definition in analysis.h has been generated automatically
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
// Root > T->Process("analysis.C")
// Root > T->Process("analysis.C","some options")
// Root > T->Process("analysis.C+")
//

#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>


void analysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void analysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   
   h1_phi = new TH1F("h1_phi","#phi #pi^{+} distribution",150,-180,360);
   h1_costheta = new TH1F("h1_costheta","cos(#theta) #pi^{+} distribution (f^{1} rest frame)",100,-1.,1.);
   h1_mass = new TH1F("h1_mass f1","Mass; GeV",100,0.,1.5);
   h1_costheta2 = new TH1F("h1_costheta2","cos(#theta) f^{1} distribution",100,-1.,1.);
   h1_theta_pim = new TH1F("h1_theta_pim","#theta #pi^{-} distribution",100,0.0,TMath::Pi());
   h1_mass2 = new TH1F("h1_mass f1 as mass of #pi^{+}+#pi^{-}+#eta","Mass; GeV",100,0.,1.5);
   h2_pimeta_pipeta = new TH2F("h2_pimeta_pipeta","Dalitz ;#pi^{-}#eta;#pi^{+}#eta",100,0.,1.5,100,0.,1.5);

   fOutput->Add(h1_phi);
   fOutput->Add(h1_costheta);
   fOutput->Add(h1_mass);
   fOutput->Add(h1_mass2);
   fOutput->Add(h1_costheta2);
   fOutput->Add(h1_theta_pim);
   fOutput->Add(h2_pimeta_pipeta);
}

Bool_t analysis::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either analysis::GetEntry() or TBranch::GetEntry()
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


  TLorentzVector p_f1(px[1],py[1],pz[1],Ef[1]);
  TLorentzVector p_pip(px[2],py[2],pz[2],Ef[2]);
  TLorentzVector p_pim(px[3],py[3],pz[3],Ef[3]);
  TLorentzVector p_eta(px[4],py[4],pz[4],Ef[4]);

  TLorentzVector p_d1 = p_pim+p_eta;
  TLorentzVector p_d2 = p_pip+p_eta;
  
  TVector3 b_3 ;
  b_3 =  p_f1.BoostVector();
  b_3 = -b_3;
  TLorentzVector p_f1_2 = p_pip+p_pim+p_eta;
  p_pip.Boost(b_3);
  
  h1_phi->Fill(p_pip.Phi()/TMath::Pi()*180.);
  h1_costheta->Fill(p_pip.CosTheta());
  h1_mass->Fill(p_f1.M());
  h1_mass2->Fill(p_f1_2.M());
  h1_costheta2->Fill(p_f1.CosTheta());
  h1_theta_pim->Fill(p_pim.Theta());
  h2_pimeta_pipeta->Fill(p_d1.M(),p_d2.M(),weight[4]);
   return kTRUE;
}

void analysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void analysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  TFile file_out("analysis_output.root","recreate");
  TList *outlist = GetOutputList();
  
  outlist->Write();
  file_out.Write();
  file_out.Close();
}
