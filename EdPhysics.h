#ifndef __EdPhysics_h
#define __EdPhysics_h


#include "EdModel.h"
#include "EdOutput.h"
#include <iostream>
#include "TLorentzVector.h"
#include "TF2.h"
#include "TF1.h"
#include "TVector3.h"
#include "EdGenPhaseSpace.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TRandom2.h"



using namespace std;

class EdPhysics: public EdGenPhaseSpace{
    public:
	EdPhysics(EdModel *);
	~EdPhysics();

	void MakeEvent(EdOutput * , EdModel *);
	Double_t Generate_event(EdModel *, int);

	enum nucl {kProton, kNeutron};
    private:

	int n_part;
	int mass_model;
	int ph_model;
	int nvertex; // total number of vertexes
	int npvert[10]; // total number of particle at one vertex
	int overt[10];  // particle at origin of one vertex (need to be already existing)
	int    v_type[10]; // vertex type (1= ratio  of probability , 2= cross section
	double v_ratio[10]; // ration to be applied to vertex
	TDatabasePDG *pdg;
	TParticlePDG *part_pdg[MAX_PART];
	TLorentzVector *p4vector[MAX_PART][MAX_PART];
	TLorentzVector *p4vector_c;


	double theta[MAX_PART];
	double phi[MAX_PART];
	double Ef[MAX_PART];
	double pf[MAX_PART];
	double px[MAX_PART];
	double py[MAX_PART];
	double pz[MAX_PART];
	double x;
	double W;
	double y;
	double t_val;
	double Q2;
	double nu;
	double e_lab;
	double weight[MAX_PART];

	int particle_id[MAX_PART];
	int charge[MAX_PART];
	int Z_ion;
	int N_ion;
	double vx[MAX_PART];
	double vy[MAX_PART];
	double vz[MAX_PART];
	
	int towrite[MAX_PART];
	double masses[10][10];
	double val_mass[10][10];
	double val_mass_t1[10];
	double val_mass_t2[10];
	double width[10][10];
	int distr_mass[10][10];
	double max_mass[10][10];
	double mass_meson;
	double t_calc;
	TRandom2 *fRandom;
	TRandom2 *gRandom;
	TLorentzVector Wtg;
	TLorentzVector beam;
	TLorentzVector gammastar;
	TLorentzVector target;
	TLorentzVector spectator;
	TLorentzVector W4vector;
	TLorentzVector Q4vector;
	TLorentzVector t4vector;
	TVector3 vertex;
	double theta_min[MAX_PART];
	double theta_max[MAX_PART];
	double energy_min[MAX_PART];
	double energy_max[MAX_PART];
	int count_phase;
	double total_t;

	TVector3 Decay_vertex(TLorentzVector *Vp_4, int i, TVector3 vert);
	double GetBeamProfile( double sigma = 1.);
	int Gen_Phasespace(EdModel *model);
	int Gen_Mass(int i,EdModel *model);
	int Gen_Mass_t(EdModel *model);
	Double_t Calc_gamma(double t_gen);
	void QFTarget(EdModel *model);  //calculate quasi free target with fermi momentum

	double t_reaction(TLorentzVector *Vrecoil_tg_4 );
	double phi_meson();
	double theta_meson();
	double phi_3decay();
	double theta_3decay();
	double phi_2decay();
	double theta_2decay();
	 
};
#endif//__EdPhysics_h
