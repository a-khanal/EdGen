#ifndef __EdModel_h
#define __EdModel_h

#include "EdInput.h"
#include "EdHisto.h"
#include <math.h>
#include <cstdlib>

#include <TVector3.h>
#include "TTree.h"
#include "TRandom2.h"

// HG electron beam

class EdModel{
    public:
	EdModel(EdInput *inp = 0);
	~EdModel();

	double GetLength() {return length; }
	double GetLx() {return len_x; }
	double GetLy() {return len_y; }
	int    Get_tgZ() {return tg_Z; }
	int    Get_tgN() {return tg_N; }
	double Get_tgMass() {return tg_mass; }
	double GetEnergy();
	const char * GetMassModelString();
	double GetBeamPID() {return beam_pid;}
	//	char* GetInFileName() {return ifile.data()};
	int    GetNpart() {return npart;}
	int    GetPhModel() {return ph_model;}
	int    GetMassModel() {return m_model;}
	TVector3 GetTgtOffset(){ return offset; }
	int    GetPid(int i){ return pid[i]; }
	int    GetNvertex(){ return nvertex; }
	int    GetOvert(int i){ return overt[i]; }
	int    GetNpvert(int i){ return npvert[i]; }
	int    GetDvert(int i, int j){ return dvert[i][j]; }
	int    GetV_type(int i){ return v_type[i]; }
	int    GetV_ratio(int i){ return v_ratio[i]; }
	double GetTheta_min(int i) {return theta_min[i];}
	double GetTheta_max(int i) {return theta_max[i];}
	Bool_t IsQF() {return fIsQF;}
	TH1F* GetFermi(){return fFermiMomentum;}
	EdInput* GetInput(){return fInp;}
	void SetRandom(TRandom2* thisRandom) {fRandom = thisRandom;};
   private:
	double length;
	double len_x;
	double len_y;
	double e_out_min;
	double e_out_max;
	int    tg_Z;
	int    tg_N;
	double tg_mass;
	int    ph_model;
	int    m_model;
	TString ifile;
	double energy;
	int    beam_pid;
	double theta_min[100];
	double theta_max[100];
	int    pid[100];
	int    npart;
	int    nvertex; // total number of vertexes
	int    npvert[10]; // total number of particle at one vertex
	int    overt[10]; // particle at origin of one vertex (need to be already existing)
	int    dvert[10][10]; // particle daughters at vertex 
	int    v_type[10]; // vertex type (1= ratio  of probability , 2= cross section
	double v_ratio[10]; // ration to be applied to vertex
	TVector3 offset;
	EdHisto   *H1_spec;
	TH1F* fFermiMomentum; //possible fermi momentum distribution for qf target
	Bool_t fIsQF;  //Check if QuasiFree model
	EdInput* fInp; //Saves having to duplicate getter functions....
	TRandom2* fRandom;
	int histo_Random_set;

};
#endif//__EdModel_h
