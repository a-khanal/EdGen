#ifndef __EdInput_h
#define __EdInput_h

#include "TObject.h"
#include "TString.h"

struct inputdata {
    int    nevt;
    int    nprnt;
    int    model;
    int    mass_model;
  TString ifile;  //input spectrum for beam
  TString qf_file; //input fermi momentum distribution file 
  TString qf_fermi; //input fermi momentum distribution histogram
  int isqf; //is the target quasifree?
  int qftpdg;  //pdg number of qf target
  int qfspdg;  //pdg number of quasifree spectator
    int    npart; // total number of particles involved
    int    nvertex; // total number of vertexes
    int    npvert[10]; // total number of particle at one vertex
    int    overt[10]; // particle at origin of one vertex (need to be already existing)
    int    dvert[10][10]; // particle daughters at vertex 
    int    v_type[10]; // vertex type (1= ratio  of probability , 2= cross section
    double v_ratio[10]; // ration to be applied to vertex

    int    pid[100];
    int    beam_pid;
    double e_energy;
    double e_energy_min;
    double e_energy_max;
    double length;
    double lenx;
    double leny;
    double sigx;
    double sigy;
 
    double theta_min[100];
    double theta_max[100];

    int    tg_Z;
    int    tg_N;
    double tg_mass;
    double off_x;
    double off_y;
    double off_z;

    int    out_fmt; // 1=root 2=root+LUND 

};

class EdInput{ 
    public:
	 EdInput(const char *file);
	~EdInput();

	int    GetNevt(){ return fData.nevt; }
	int    GetNprnt(){ return fData.nprnt; }
	int    GetNpart(){ return fData.npart; }
	int    GetModel(){ return fData.model; }
	int    GetMassModel(){ return fData.mass_model; }
      TString  GetIfile(){ return fData.ifile; }
      TString  GetQFFile(){ return fData.qf_file; }
      TString  GetQFFermi(){ return fData.qf_fermi; }
      int IsQF(){return fData.isqf;}
      inputdata GetInData(){return fData;}
	int    GetBeamPID(){ return fData.beam_pid; }
	double GetSigBeamX(){ return fData.sigx; }
	double GetSigBeamY(){ return fData.sigy; }
	int    GetNvertex(){ return fData.nvertex; }
	int    GetOvert(int i){ return fData.overt[i]; }
	int    GetNpvert(int i){ return fData.npvert[i]; }
	int    GetDvert(int i, int j){ return fData.dvert[i][j]; }
	int    GetV_type(int i){ return fData.v_type[i]; }
	int    GetV_ratio(int i){ return fData.v_ratio[i]; }
	int    GetPid(int i){ return fData.pid[i]; }
	double Get_eEnergy(){ return fData.e_energy; }
	double GetEnergy_min(){ return fData.e_energy_min; }
	double GetEnergy_max(){ return fData.e_energy_max; }
	double Get_length(){ return fData.length; }
	int    Get_tg_Z(){ return fData.tg_Z; }
	int    Get_tg_N(){ return fData.tg_N; }
	double Get_tg_mass(){ return fData.tg_mass; }
	double Get_lenx(){ return fData.lenx; }
	double Get_leny(){ return fData.leny; }
	int    Get_fmt(){ return fData.out_fmt; }
	double Get_thetaMin(int i){ return fData.theta_min[i]; }
	double Get_thetaMax(int i){ return fData.theta_max[i]; }

	double GetTgtXoff(){ return fData.off_x; }
	double GetTgtYoff(){ return fData.off_y; }
	double GetTgtZoff(){ return fData.off_z; }
	


    private:
	inputdata fData;
	int numpart;
};
#endif//__EdInput_h
