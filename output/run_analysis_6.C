{
  TProof::Open("workers=2");
  TChain *mc_edgen = new TChain("T");
  mc_edgen->Add("output.root");
  mc_edgen->SetProof();
  mc_edgen->Process("analysis_6.C++");


}
