#include "EdModel.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TMath.h"
#include <iostream>

EdModel::EdModel(EdInput *inp){
  fInp=inp;
  fFermiMomentum=0;
  fIsQF=kFALSE;
  H1_spec=0;
  e_out_min = 0.;
  e_out_max = 0.;
  fRandom = 0;
  histo_Random_set = 0;
  histo_tRandom_set = 0;
  histo_qRandom_set = 0;
  histo_eRandom_set = 0;

    if( inp ){
      int tot_part = 100;
      ifile = inp->GetIfile();
      tfile = inp->GetTfile();
      qfile = inp->GetQfile();
      efile = inp->GetEfile();
      length = inp->Get_length();
      len_x = inp->Get_lenx();
      len_y = inp->Get_leny();
      ph_model = inp->GetModel();
      m_model = inp->GetMassModel();
      beam_pid = inp->GetBeamPID();
	if (ph_model == 2) {
	  Float_t Energy_1, Energy_2, E_counts;
	  TTree *Input_spectrum = new TTree("Hin", "HG Monte Carlo input");
	  Input_spectrum->Branch("Energy_1",&Energy_1,"Energy_1/F");
	  Input_spectrum->Branch("Energy_2",&Energy_2,"Energy_2/F");
	  Input_spectrum->Branch("E_counts",&E_counts,"E_counts/F");
	  printf("Reading input file %s\n",ifile.Data());
	  Input_spectrum->ReadFile(ifile.Data(), "Energy_1:Energy_2:E_counts");
	  H1_spec = new EdHisto("H1_spec","H1_spec",Input_spectrum->GetEntries(),Input_spectrum->GetMinimum("Energy_1"),Input_spectrum->GetMaximum("Energy_2"));
	  Axis_t *new_bins = new Axis_t[Input_spectrum->GetEntries() + 1];	    
	  TAxis *axis = H1_spec->GetXaxis(); 
	  for (int i=0; i< Input_spectrum->GetEntries(); i++) {
	    Input_spectrum->GetEntry(i);
	    new_bins[i] = Energy_1;
	    H1_spec->SetBinContent(i+1,E_counts);
	    if (i+1 == Input_spectrum->GetEntries()) new_bins[i+1] = Energy_2; 
	  }
	  axis->Set(Input_spectrum->GetEntries(), new_bins); 
	  delete new_bins; 
	  delete Input_spectrum;
	}
	if (ph_model == 3) {
	  e_out_min = inp->GetEnergy_min();
	  e_out_max = inp->GetEnergy_max();
	  printf("Set Energy range for beam from %.6f GeV to %.6f GeV \n",e_out_min,e_out_max); 
	}
	if (ph_model == 5) {
	  Float_t Energy_1, Energy_2, E_counts;
	  TTree *Input_spectrum2 = new TTree("Hin", "HG Monte Carlo input");
	  Input_spectrum2->Branch("Energy_1",&Energy_1,"Energy_1/F");
	  Input_spectrum2->Branch("Energy_2",&Energy_2,"Energy_2/F");
	  Input_spectrum2->Branch("E_counts",&E_counts,"E_counts/F");
	  printf("Reading input file for t distribution %s\n",tfile.Data());
	  Input_spectrum2->ReadFile(tfile.Data(), "Energy_1:Energy_2:E_counts");
	  H1_tspec = new EdHisto("H1_tspec","H1_tspec",Input_spectrum2->GetEntries(),Input_spectrum2->GetMinimum("Energy_1"),Input_spectrum2->GetMaximum("Energy_2"));
	  Axis_t *new_bins2 = new Axis_t[Input_spectrum2->GetEntries() + 1];	    
	  TAxis *axis = H1_tspec->GetXaxis(); 
	  for (int i=0; i< Input_spectrum2->GetEntries(); i++) {
	    Input_spectrum2->GetEntry(i);
	    new_bins2[i] = Energy_1;
	    H1_tspec->SetBinContent(i+1,E_counts);
	    if (i+1 == Input_spectrum2->GetEntries()) new_bins2[i+1] = Energy_2; 
	  }
	  axis->Set(Input_spectrum2->GetEntries(), new_bins2); 
	  delete new_bins2; 
	  delete Input_spectrum2;
	  TTree *Input_spectrum3 = new TTree("Hin", "HG Monte Carlo input");
	  Input_spectrum3->Branch("Energy_1",&Energy_1,"Energy_1/F");
	  Input_spectrum3->Branch("Energy_2",&Energy_2,"Energy_2/F");
	  Input_spectrum3->Branch("E_counts",&E_counts,"E_counts/F");
	  printf("Reading input file for q2 distribution %s\n",qfile.Data());
	  Input_spectrum3->ReadFile(qfile.Data(), "Energy_1:Energy_2:E_counts");
	  H1_qspec = new EdHisto("H1_qspec","H1_qspec",Input_spectrum3->GetEntries(),Input_spectrum3->GetMinimum("Energy_1"),Input_spectrum3->GetMaximum("Energy_2"));
	  Axis_t *new_bins3 = new Axis_t[Input_spectrum3->GetEntries() + 1];	    
	  TAxis *axisq = H1_qspec->GetXaxis(); 
	  for (int i=0; i< Input_spectrum3->GetEntries(); i++) {
	    Input_spectrum3->GetEntry(i);
	    new_bins3[i] = Energy_1;
	    H1_qspec->SetBinContent(i+1,E_counts);
	    if (i+1 == Input_spectrum3->GetEntries()) new_bins3[i+1] = Energy_2; 
	  }
	  axisq->Set(Input_spectrum3->GetEntries(), new_bins3); 
	  delete new_bins3; 
	  delete Input_spectrum3;
	  TTree *Input_spectrum4 = new TTree("Hin", "HG Monte Carlo input");
	  Input_spectrum4->Branch("Energy_1",&Energy_1,"Energy_1/F");
	  Input_spectrum4->Branch("Energy_2",&Energy_2,"Energy_2/F");
	  Input_spectrum4->Branch("E_counts",&E_counts,"E_counts/F");
	  printf("Reading input file for Energy e' distribution %s\n",efile.Data());
	  Input_spectrum4->ReadFile(efile.Data(), "Energy_1:Energy_2:E_counts");
	  H1_espec = new EdHisto("H1_espec","H1_espec",Input_spectrum4->GetEntries(),Input_spectrum4->GetMinimum("Energy_1"),Input_spectrum4->GetMaximum("Energy_2"));
	  Axis_t *new_bins4 = new Axis_t[Input_spectrum4->GetEntries() + 1];	    
	  TAxis *axise = H1_espec->GetXaxis(); 
	  for (int i=0; i< Input_spectrum4->GetEntries(); i++) {
	    Input_spectrum4->GetEntry(i);
	    new_bins4[i] = Energy_1;
	    H1_espec->SetBinContent(i+1,E_counts);
	    if (i+1 == Input_spectrum4->GetEntries()) new_bins4[i+1] = Energy_2; 
	  }
	  axise->Set(Input_spectrum4->GetEntries(), new_bins4); 
	  delete new_bins4; 
	  delete Input_spectrum4;
	}
	tg_Z = inp->Get_tg_Z();
	tg_N = inp->Get_tg_N();
	tg_mass = inp->Get_tg_mass();
	energy = inp->Get_eEnergy();
	npart = inp->GetNpart();
	tot_part = npart;
	for (int i=0; i<tot_part; i++) {
	  pid[i] = inp->GetPid(i);
	  theta_min[i] = inp->Get_thetaMin(i);
	  theta_max[i] = inp->Get_thetaMax(i);
	  energy_min[i] = inp->Get_energyMin(i);
	  energy_max[i] = inp->Get_energyMax(i);
	}
	nvertex = inp->GetNvertex();
	tot_part = nvertex;
	for (int i=0; i<tot_part; i++) {
	  npvert[i] = inp->GetNpvert(i);
	  overt[i] = inp->GetOvert(i);
	  v_ratio[i] = inp->GetV_ratio(i);
	  v_type[i] = inp->GetV_type(i);
	}
	offset.SetXYZ( inp->GetTgtXoff(),
		       inp->GetTgtYoff(),
		       inp->GetTgtZoff() );

	if(inp->IsQF()){
	  // Get Fermi Momentum Distribution
	  TDirectory* savedir=gDirectory;//standard ROOT hack for not losing newly created objects when close file
	  TFile* fermidat = new TFile(inp->GetQFFile());
	  savedir->cd();
	  if(!fermidat->IsOpen()) std::cerr<<"EdModel::EdModel(EdInput *inp) : No Quasi Free file found = "<<inp->GetQFFile()<<std::endl;
	  else {
	    std::cout << "QF distribution : " << inp->GetQFFermi()<<std::endl;
	    fFermiMomentum  = (TH1F*)fermidat->Get(inp->GetQFFermi())->Clone("hFermi");
	    if(!fFermiMomentum){ std::cout<<"Did not find "<<inp->GetQFFermi()<<" in "<<inp->GetQFFile()<<std::endl; exit(-1);}
	  }
	  if(fFermiMomentum)   fIsQF=kTRUE;
	  fermidat->Close();
	  delete fermidat;
	  
	}
    }
    
    length = length / 100. ; // conversion distances in m
    len_x = len_x / 100. ; 
    len_y = len_y / 100. ; 
    offset *= 1e-2; // conversion distances in cm to m

    return;
}

EdModel::~EdModel(){
  if(H1_spec) delete H1_spec;
  if(H1_tspec) delete H1_tspec;
  if(fFermiMomentum) delete fFermiMomentum;
  return;
}


double EdModel::GetEnergy(){
  double e_out = 0.;
  if (ph_model == 1 || ph_model==5) { // PhaseSpace Single Energy
    e_out = energy;

  }
  else if (ph_model == 2) { // PhaseSpace Multiple Energy
    //    printf("here 1 %d\n",histo_Random_set);
    if (histo_Random_set == 0 || H1_spec->GetRandom2() == 0) {
      histo_Random_set = 1;
      // printf("here 2 \n");
      H1_spec->SetRandom(fRandom);
      // printf("here 3 \n");
    
    }
    while (TMath::IsNaN(e_out) || e_out ==0) e_out = H1_spec->GetEdRandom();
  }
  else if (ph_model == 3) { // PhaseSpace Flat multiple Energy
    e_out = fRandom->Uniform(e_out_min,e_out_max);
  }
  return e_out;
}

double EdModel::Get_tvalue(){
  double e_out = 0.;
  if (ph_model == 5) { // PhaseSpace Multiple Energy
    //    printf("here 1 %d\n",histo_Random_set);
    if (histo_tRandom_set == 0 || H1_tspec->GetRandom2() == 0) {
      histo_tRandom_set = 1;
      // printf("here 2 \n");
      H1_tspec->SetRandom(fRandom);
      // printf("here 3 \n");
    
    }
    while (TMath::IsNaN(e_out) || e_out ==0) e_out = H1_tspec->GetEdRandom();
  }
  return e_out;
}

double EdModel::Get_qvalue(){
  double e_out = 0.;
  if (ph_model == 5) { // PhaseSpace Multiple Energy
    //    printf("here 1 %d\n",histo_Random_set);
    if (histo_qRandom_set == 0 || H1_qspec->GetRandom2() == 0) {
      histo_qRandom_set = 1;
      // printf("here 2 \n");
      H1_qspec->SetRandom(fRandom);
      // printf("here 3 \n");
    
    }
    while (TMath::IsNaN(e_out) || e_out ==0) e_out = H1_qspec->GetEdRandom();
  }
  return e_out;
}

double EdModel::Get_evalue(){
  double e_out = 0.;
  if (ph_model == 5) { // PhaseSpace Multiple Energy
    //    printf("here 1 %d\n",histo_Random_set);
    if (histo_eRandom_set == 0 || H1_espec->GetRandom2() == 0) {
      histo_eRandom_set = 1;
      // printf("here 2 \n");
      H1_espec->SetRandom(fRandom);
      // printf("here 3 \n");
    
    }
    while (TMath::IsNaN(e_out) || e_out ==0) e_out = H1_espec->GetEdRandom();
  }
  return e_out;
}

const char * EdModel::GetMassModelString(){
  if (m_model == 1) return "Breit-Wigner";
  else if (m_model == 2) return "Flat";
  else if (m_model == 3) return "m=mass";
  else return "Sorry: No mass model supported. Check your input file";

}

