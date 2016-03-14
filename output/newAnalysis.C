#define newAnalysis_cxx
// The class definition in newAnalysis.h has been generated automatically
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
// Root > T->Process("newAnalysis.C")
// Root > T->Process("newAnalysis.C","some options")
// Root > T->Process("newAnalysis.C+")
//

#include "newAnalysis.h"
#include <TH2.h>
#include <TStyle.h>


void newAnalysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void newAnalysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   h1_phi = new TH1F("h1_phi","#phi #pi^{+} distribution (K^{0}_{S}(1) rest frame)",150,-180,360);
   h1_costheta = new TH1F("h1_costheta","cos(#theta) #pi^{+} distribution (K^{0}_{S}(1) rest frame)",100,-1.,1.);
   h1_mass = new TH1F("h1_mass_f","Mass f'_{2}; GeV",100,0.,2.0);
   h1_mass_K0_1 = new TH1F("h1_mass_K0_1","Mass K^{0}_{S}(1); GeV",100,0.,1.5);
   h1_mass_K0_2 = new TH1F("h1_mass_K0_2","Mass K^{0}_{S}(2); GeV",100,0.,1.5);
   h1_vz_K0_1 = new TH1F("h1_vz_K0_1","Vertex K^{0}_{S}(1); m",530,-0.3,5.0);
   h1_vz_K0_2 = new TH1F("h1_vz_K0_2","Vertex K^{0}_{S}(2); m",530,-0.3,5.0);
   h1_vz_piK0_1 = new TH1F("h1_vz_piK0_1","Vertex #pi^{+} from K^{0}_{S}(1); m",530,-0.3,5.0);
   h1_vz_piK0_2 = new TH1F("h1_vz_piK0_2","Vertex #pi^{+} from K^{0}_{S}(2); m",530,-0.3,5.0);

   h1_mass2_K0_1 = new TH1F("h1_mass2_K0_1","Mass K^{0}_{S}(1) from #pi^{+} + #pi^{-}; GeV",100,0.,1.5);
   h1_mass2_K0_2 = new TH1F("h1_mass2_K0_2","Mass K^{0}_{S}(2) from #pi^{+} + #pi^{-}; GeV",100,0.,1.5);
   h1_costheta2 = new TH1F("h1_costheta2","cos(#theta) f'_{2} distribution",100,-1.,1.);
   h1_theta_pim = new TH1F("h1_theta_pim","#theta #pi^{-} distribution",100,0.0,TMath::Pi());
   h1_mass2 = new TH1F("h1_mass_f"," f'_{2} as mass of K^{0}_{S}(1)+K^{0}_{S}(2); GeV",100,0.,2.0);
   h2_mass_ebeam = new TH2F("h2_mass_ebeam","mass f'_{2} vs Energy beam; Energy beam (GeV); mass f'_{2} (GeV)",100,0.0,4.0,100,0.0,2.0);

   fOutput->Add(h1_phi);
   fOutput->Add(h1_costheta);
   fOutput->Add(h1_mass);
   fOutput->Add(h1_mass_K0_1);
   fOutput->Add(h1_mass_K0_2);
   fOutput->Add(h1_mass2_K0_1);
   fOutput->Add(h1_mass2_K0_2);
   fOutput->Add(h1_vz_K0_1);
   fOutput->Add(h1_vz_K0_2);
   fOutput->Add(h1_vz_piK0_1);
   fOutput->Add(h1_vz_piK0_2);
   fOutput->Add(h1_mass2);
   fOutput->Add(h1_costheta2);
   fOutput->Add(h1_theta_pim);
   fOutput->Add(h2_mass_ebeam);


}

Bool_t newAnalysis::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either newAnalysis::GetEntry() or TBranch::GetEntry()
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

  b_px->GetEntry(entry);
  b_py->GetEntry(entry);
  b_pz->GetEntry(entry);
  b_vz->GetEntry(entry);
  b_Ef->GetEntry(entry);
  b_Ein_beam->GetEntry(entry);

  TLorentzVector p_fp2(px[1],py[1],pz[1],Ef[1]);
  TLorentzVector p_K0_1(px[2],py[2],pz[2],Ef[2]);
  TLorentzVector p_K0_2(px[3],py[3],pz[3],Ef[3]);
  TLorentzVector p_pip_1(px[4],py[4],pz[4],Ef[4]);
  TLorentzVector p_pim_1(px[5],py[5],pz[5],Ef[5]);
  TLorentzVector p_pip_2(px[6],py[6],pz[6],Ef[6]);
  TLorentzVector p_pim_2(px[7],py[7],pz[7],Ef[7]);


  TLorentzVector p_fp2_2 = p_K0_1+p_K0_2;
  TLorentzVector p_K0_1_2 = p_pip_1+p_pim_1;
  TLorentzVector p_K0_2_2 = p_pip_2+p_pim_2;

  TVector3 b_3 ;
  b_3 =  p_K0_1.BoostVector();
  b_3 = -b_3;
  p_pip_1.Boost(b_3);

  h1_phi->Fill(p_pip_1.Phi()/TMath::Pi()*180.);
  h1_costheta->Fill(p_pip_1.CosTheta());
  h1_mass->Fill(p_fp2.M());
  h1_mass2->Fill(p_fp2_2.M());
  h1_mass_K0_1->Fill(p_K0_1.M());
  h1_mass_K0_2->Fill(p_K0_2.M());
  h1_mass2_K0_1->Fill(p_K0_1_2.M());
  h1_mass2_K0_2->Fill(p_K0_2_2.M());
  h1_costheta2->Fill(p_fp2.CosTheta());
  h1_theta_pim->Fill(p_pim_1.Theta());
  h2_mass_ebeam->Fill(Ein_beam,p_fp2.M());
  h1_vz_K0_1->Fill(vz[2]);
  h1_vz_K0_2->Fill(vz[3]);
  h1_vz_piK0_1->Fill(vz[4]);
  h1_vz_piK0_2->Fill(vz[6]);

   return kTRUE;
}

void newAnalysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void newAnalysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.


  TFile file_out("analysis_new_output.root","recreate");
  TList *outlist = GetOutputList();
  
  outlist->Write();
  file_out.Write();
  file_out.Close();

}
