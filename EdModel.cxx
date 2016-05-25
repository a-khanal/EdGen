#include "EdModel.h"
#include "TFile.h"
#include "TDirectory.h"
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

    if( inp ){
      int tot_part = 100;
      ifile = inp->GetIfile();
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
  if(fFermiMomentum) delete fFermiMomentum;
  return;
}


double EdModel::GetEnergy(){
  double e_out = 0.;
  if (ph_model == 1) { // PhaseSpace Single Energy
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
    while (isnan(e_out) || e_out ==0) e_out = H1_spec->GetEdRandom();
  }
  else if (ph_model == 3) { // PhaseSpace Flat multiple Energy
    e_out = fRandom->Uniform(e_out_min,e_out_max);
  }
  return e_out;
}

const char * EdModel::GetMassModelString(){
  if (m_model == 1) return "Breit-Wigner";
  else if (m_model == 2) return "Flat";
  else if (m_model == 3) return "m=mass";
  else return "Sorry: No mass model supported. Check your input file";

}

