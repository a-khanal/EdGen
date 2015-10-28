{
  TProof::Open("workers=2");
  TChain *mc_edgen = new TChain("T");
  mc_edgen->Add("test_K0L.root");
  mc_edgen->SetProof();
  mc_edgen->Process("analysis.C++");


}
