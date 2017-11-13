#include "EdOutputBos.h"
using namespace std;

EdOutput::EdOutput(EdInput *inp, const char *fileout){
    char defaultname[255] = "output.root";

    if( fileout ){
	strcpy(fOutName, fileout);
    } else {
	strcpy(fOutName, defaultname);
    }
    printf("output file set as %s\n", fOutName);

    fTree = new TTree("T", "HG Monte Carlo");
    fNevt    = ((double) inp->GetNevt());
    n_part = inp->GetNpart();
    fnvertex = inp->GetNvertex();
    printf("Total weight of event will be = "); 
    for (int i=0; i< fnvertex; i++) {
      if (i==0) f1part[i] = 0;
      else {
	f1part[i] = f1part[i-1] + inp->GetNpvert(i-1);
	printf(" x ");
      }
      printf("w[%i]",f1part[i]);
    }
    printf("\n");
    int k=0;
    for (int i=0; i< fnvertex;i++) {
      for (int j=0; j<inp->GetNpvert(i) ; j++) {
	overt[k] = inp->GetOvert(i);
	k++;
      }
    }

    
    InitTree();

    return;
}

EdOutput::~EdOutput(){
    delete fTree;

    return;
}

void EdOutput::InitTree(){
    fTree->Branch("n_part", &n_part, "n_part/I");
    fTree->Branch("theta", theta, "theta[n_part]/D");
  
    fTree->Branch("phi", phi, "phi[n_part]/D");
    fTree->Branch("x", &x, "x/D");
    fTree->Branch("Ein_beam", &Ein, "Ein_beam/D");
    fTree->Branch("Ef", Ef, "Ef[n_part]/D");
 
    fTree->Branch("Q2", &Q2, "Q2/D");
    fTree->Branch("nu", &nu, "nu/D");

    fTree->Branch("W", &W, "W/D");
    fTree->Branch("t", &t_val, "t/D");
    fTree->Branch("y", &y, "y/D");

  
    fTree->Branch("Z_ion", &Z_ion, "Z_ion/I");
    fTree->Branch("N_ion", &N_ion, "N_ion/I");
    
    fTree->Branch("particle_id", particle_id, "particle_id[n_part]/I");
    fTree->Branch("charge", charge, "charge[n_part]/I");

    fTree->Branch("pf", pf, "pf[n_part]/D");
    fTree->Branch("px", px, "px[n_part]/D");
    fTree->Branch("py", py, "py[n_part]/D");
    fTree->Branch("pz", pz, "pz[n_part]/D");

    fTree->Branch("vx", &vx, "vx[n_part]/D");
    fTree->Branch("vy", &vy, "vy[n_part]/D");
    fTree->Branch("vz", &vz, "vz[n_part]/D");

    fTree->Branch("weight", &weight, "weight[n_part]/D");


    return;
}

void EdOutput::SetTheta( double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) theta[i] = 0; 
    else theta[i] = val[i];
  }
  n_part = tot;
}
void EdOutput::SetPhi(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) phi[i] = 0; 
    else phi[i] = val[i];
  }
}
void EdOutput::SetEf(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) Ef[i] = 0; 
    else Ef[i] = val[i];
  }

}
void EdOutput::Setpf(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) pf[i] = 0; 
    else pf[i] = val[i];
  }

}
void EdOutput::Setpx(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) px[i] = 0; 
    else px[i] = val[i];
  }

}
void EdOutput::Setpy(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) py[i] = 0; 
    else py[i] = val[i];
  }

}
void EdOutput::Setpz(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) pz[i] = 0; 
    else pz[i] = val[i];
  }

}
void EdOutput::Setparticle_id(int *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) particle_id[i] = 0; 
    else particle_id[i] = val[i];
  }

}
void EdOutput::Setcharge(int *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) charge[i] = 0; 
    else charge[i] = val[i];
  }

}
void EdOutput::Setvx(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) vx[i] = 0; 
    else vx[i] = val[i];
  }

}
void EdOutput::Setvy(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) vy[i] = 0; 
    else vy[i] = val[i];
  }

}
void EdOutput::Setvz(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) vz[i] = 0; 
    else vz[i] = val[i];
  }

}
void EdOutput::Setweight(double *val, int tot){
  for (int i =0; i<tot; i++) {
    if (std::isnan(val[i])) weight[i] = 0; 
    else weight[i] = val[i];
  }

}
void EdOutput::Settowrite(int *val, int tot){
  for (int i =0; i<tot; i++) {
    towrite[i] = val[i];
  }

}

void EdOutput::Write(){
    fTree->Fill();

    return;
}

void EdOutput::Close(){
    fOutfile = new TFile(fOutName, "RECREATE");
    fTree->Write();
    fOutfile->Close();
    delete fOutfile;

    return;
}

void  EdOutput::MakeFileLUND(){


  

  TString file(fOutName);
  file.ReplaceAll("root","lund"); 
  char outstring[200];

  
  Int_t nentries = (Int_t)fTree->GetEntries();
  int tot_part = 0;
  for (int j=0; j<n_part; j++) {
    tot_part = tot_part +towrite[j];
  }

  double vxcm,vycm,vzcm;
  double fweight;
  int active;
  int daughter;
  int at_daughter;
  
  std::ofstream OUT (file.Data());
  for (int i=0; i<nentries ; i++) {
    fTree->GetEntry(i);
    if(i % 1000000 == 0 ){
      printf("Analyzed %09d events of total %09d \n",i,nentries);
    }
    sprintf(outstring,"%i %i %i 0 0 %1.4e %1.4e %1.4e %1.4e %1.4e",n_part,(Z_ion+N_ion),Z_ion,x,y,W,Q2,nu);
    OUT << outstring << endl;
    fweight = 1.0;
    at_daughter = 0;
    for (int j=0; j<fnvertex; j++) {
      fweight = weight[f1part[j]] * fweight;
    }

    //    OUT << tot_part << " " << (Z_ion + N_ion)  << " " << Z_ion  << " " << "0"  << " " << "0" << " "  << x << " " << y  << " \t " << W  << " \t " << Q2  << " \t " << nu << endl;
    for (int j=0; j<n_part; j++) {
      if (towrite[j] == 1) active = 1;
      else active = 0;
      vxcm = vx[j]*100.0;
      vycm = vy[j]*100.0;
      vzcm = vz[j]*100.0;
      daughter = 0;
      if (active == 0) {
	at_daughter++;
	daughter = f1part[at_daughter] + 1;
      }
      
      sprintf(outstring,"  %i %i %i %i %i %i %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e",j+1,charge[j],active,particle_id[j],overt[j],daughter,px[j],py[j],pz[j],Ef[j],fweight,vxcm,vycm,vzcm); 
      OUT << outstring << endl;
    }
  }

  OUT.close();

}


void  EdOutput::MakeFileBOS(){

#ifdef CLAS6LIB
#define h HEAD->head[0]
#endif
  TString file(fOutName);
  file.ReplaceAll("root","bos"); 
#ifdef CLAS6LIB  
  clasHEAD_t *HEAD;
  clasMCEV_t *MCEV;
  clasMCTK_t *MCTK;
  clasMCVX_t *MCVX;
  clasMCHD_t *MCHD;
  clasTAGR_t *TAGR;

  
  int BosOutputUnitNo = 2; // unit to open for writing (!?!?!?!? FORTRAN)
  int maxFileLength = 2000000;  // Maximum length of the bos file (bos files can be maximum this size)
  int nFileWrite = 0;
  int icode;
  int mctk_array_n;
  char mess[1024];
  double fweight;
  double c_speed = 299792458.; // speed of light in meter per second

  Int_t nentries = (Int_t)fTree->GetEntries();
  int tot_part;

  printf("Output Bos file: %s\n",file.Data());

  sprintf(mess," OPEN BOSOUTPUT UNIT=%3d FILE=\"%s\" WRITE STATUS=NEW RECL=32760 SPLITMB=2047 FORM=BINARY", BosOutputUnitNo,file.Data());
  fparm_c(mess);
  initbos(); // c_bos_io format

  bankList(&bcs_, "C=","HEADMCEVMCTKMCVXMCHDTAGR");  // Write HEAD MCTK MCVX banks into the bos file

   for (int i=0; i<nentries ; i++) {
    fTree->GetEntry(i);
    if(i % 1000000 == 0 ){
       printf("Analyzed %09d events of total %09d \n",i,nentries);
     }
    tot_part = 0;
    // will need to modify the output array in order to take care that some particle can be in the final state or not
    for (int j=0; j<n_part; j++) {
      tot_part = tot_part +towrite[j];
    }

    //   // Filling the array for bcs_ from the TTree
    HEAD = (clasHEAD_t *) makeBank(&bcs_,"HEAD",0,8,1); // void *makeBank(BOSbank *bcs, char *bankname, int banknum, int ncol, int nrow)
    MCEV = (clasMCEV_t *) makeBank(&bcs_,"MCEV",0,2,1); 
    MCTK = (clasMCTK_t *) makeBank(&bcs_,"MCTK",0,11,tot_part); // void *makeBank(BOSbank *bcs, char *bankname, int banknum, int ncol, int nrow)  
    MCVX = (clasMCVX_t *) makeBank(&bcs_,"MCVX",0,5,tot_part); // void *makeBank(BOSbank *bcs, char *bankname, int banknum, int ncol, int nrow)
    MCHD = (clasMCHD_t *) makeBank(&bcs_,"MCHD",0,16,1); // void *makeBank(BOSbank *bcs, char *bankname, int banknum, int ncol, int nrow)
    TAGR = (clasTAGR_t *) makeBank(&bcs_,"TAGR",0,6,1); // void *makeBank(BOSbank *bcs, char *bankname, int banknum, int ncol, int nrow)

    //    gettimeofday(&tp, NULL);
    //    int32_t time_sec = time(NULL) -1400000000;
    //    printf("Time = %i \n",time_sec);
    time(&secs);
    struct tm y2k = {0};
    int seconds;

    y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
    y2k.tm_year = 70;  y2k.tm_mon = 0; y2k.tm_mday = 1;

    seconds = int(difftime(secs,mktime(&y2k)));

    h.version = 1;
    h.nrun = 10; // gsim run
    h.nevent = i+1; // number of event
    h.time = seconds; // time in seconds since Jan 1, 1970
    h.type = -2;
    h.roc = 0;
    h.evtclass = 7;
    h.trigbits = 1;
    mctk_array_n = 0;
    
    
    MCEV->mcev[0].i1=(int)(rand()*100000);
    MCEV->mcev[0].i2=(int)(rand()*100000); 
    
    fweight = 1.0;
    for (int j=0; j<fnvertex; j++) {
      fweight = weight[f1part[j]] * fweight;
    }
    MCHD->mchd[0].nrun = 10; // gsim run
    MCHD->mchd[0].nevent = i+1; // number of event
    MCHD->mchd[0].type = -2;
    MCHD->mchd[0].weight = fweight;
    MCHD->mchd[0].w = W;
    MCHD->mchd[0].q2 = Q2;
 
    TAGR->tagr[0].erg = Ein;
    TAGR->tagr[0].stat = 7;
    TAGR->tagr[0].ttag = -vz[0] / c_speed * pow(10,9);  // time at center of CLAS respect to the interaction time (ns) 
    TAGR->tagr[0].tpho = -vz[0] / c_speed * pow(10,9);
    
 

    for (int j=0; j<n_part; j++) {
      if (towrite[j] == 1) {

	MCVX->mcvx[mctk_array_n].x =vx[j]*100.0;
	MCVX->mcvx[mctk_array_n].y =vy[j]*100.0;
	MCVX->mcvx[mctk_array_n].z =vz[j]*100.0;
	
	MCTK->mctk[mctk_array_n].id = particle_id[j];
	MCTK->mctk[mctk_array_n].cx = px[j]/pf[j];
	MCTK->mctk[mctk_array_n].cy = py[j]/pf[j];
	MCTK->mctk[mctk_array_n].cz = pz[j]/pf[j];
	MCTK->mctk[mctk_array_n].pmom = pf[j];
	MCTK->mctk[mctk_array_n].mass = pow( (pow(Ef[j],2)-pow(pf[j],2)) ,0.5) ;
	MCTK->mctk[mctk_array_n].charge = charge[j];
	MCTK->mctk[mctk_array_n].beg_vtx = 1;
	MCTK->mctk[mctk_array_n].end_vtx = 0;
	MCTK->mctk[mctk_array_n].parent = 1;
	MCTK->mctk[mctk_array_n].flag = 0;
	mctk_array_n++;
	// One can still declare better the beg_vtx, end_vtx, parent
	
      }
    }


   // Writing into bos file
    //   printf("at event %i \n",i);
    icode = putBOS(&bcs_, BosOutputUnitNo, "C");
    if(!icode){
      fprintf(stdout,"ERROR - Trouble writing out BOS bank. \n");
    }

    dropAllBanks(&bcs_,"C");
    cleanBanks(&bcs_);

   }
   close_fpack_unit("BOSOUTPUT");
#endif
}

