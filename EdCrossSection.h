#ifndef __EdCrossSection_h
#define __EdCrossSection_h

// This Code is an adaptation from the code that we got from Vincent Mathieu
// "As the publication time scale with CLAS is greater than the postdoc life-time, 
// I propose that experimentalists offer beers in exchange of theoretical postdoc work as an incentive for them." 
// You know now what’s the prize for using this code ;)

#include "EdModel.h"
#include "TGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>    	/* Standard Library of Complex Numbers */

// Structure with all the kinematical quantities
struct Kin
{
	double m1;				// Mass particle 1
	double m2;				// Mass particle 2
	double m3;				// Mass particle 3
	double m4;				// Mass particle 4
	double complex s;				// Mandelstam s variable
    double complex t;				// Mandelstam t variable
    double complex phi;				// phi = stu + ...
    double complex ks;		// photon  momentum in s-channel frame
    double complex qs;		// pion    momentum in s-channel frame
    double complex kt;		// photon  momentum in t-channel frame
    double complex pt;		// nucleon momentum in t-channel frame
    double complex E2s;		// initial nucleon energy in s-channel frame
    double complex E4s;		// final   nucleon energy in s-channel frame
    double complex E2t;		// initial nucleon energy in t-channel frame
    double complex E4t;		// final   nucleon energy in t-channel frame
    double complex coss;	// cos \theta_s
    double complex sins;	// sin \theta_s
    double complex cossh;	// cos \theta_s/2
    double complex sinsh;	// sin \theta_s/2
    double complex cost;	// cos \theta_t
    double complex sint;	// sin \theta_t
    double complex costh;	// cos \theta_t/2
    double complex sinth;	// sin \theta_t/2
};


class EdCrossSection{
    public:
	EdCrossSection(EdModel *);
	~EdCrossSection();

 /*
 ============================================================================
 Description of the subroutines:
 Pi0PhotCS_S
 	 differential cross section from s-channel amplitudes
 Pi0PhotAmpS
 	 s-channel helicity amplitudes with inputs = four vectors and helicities
 CGLN_Ai
 	 CGLN Scalar amplitudes Ai - the model has to be specified there
 CGLNA2F
 	 transformation between CGLN Fi (ouputs) and Ai (inputs)
 kin2to2
 	 four vectors of a + b --> c + d for given Ecm and theta_s
 kinematics
 	 compute all kinematical quantities (s- and t-channel)
 	 and store then in a struct Kin

 functions:
 ---------
 Cos2T
 	 return Mandelstam t from cos(theta_s)
 cgamma
 	 return gamma(z) or log( gamma(z) )
 lambda
 	 return a*a + b*b + c*c - 2*(a*b + b*c + c*a)
 snorm
 	 return p[0]*p[0] - ( p[1]*p[1] + p[2]*p[2] + p[3]*p[3] );
 ============================================================================
 */
	
	void kin2to2(double Ecm, double theta, double mass2[5], double pa[4],double pb[4],double pc[4], double pd[4]);
	double complex lambda(double complex a, double b, double c);
	double lambda_real(double a, double b, double c);
	double snorm(double p[4]);
	double complex cgamma(double complex z,int OPT);
	double complex Pi0PhotAmpS(double pa[4],double pb[4],double pc[4], int hel[4]);
	void Regge_CGLNAi(double complex s, double complex t, double complex CGLNA[5]);
	double Pi0PhotCS_S(double E,double theta);
	double Cos2T(double E, double z);
	void CGLNA2F(struct Kin var, double complex Ai[5], double complex Fi[5]);
	void kinematics(double complex s, double complex t, double mass2[], struct Kin *var);
	void MltPole2F(double s, double t, double complex MltPole[5][10] ,double complex CGLNF[5]);
	void MAID_Pi0(double s, double complex MltPole[5][10]);
	double Coef(double E);
	TGraph* Get_graph(){return fgraph;};


   private:
	double MP;
	double MPI;
	double EPS;
	double mass[5];
	TGraph *fgraph;

};
#endif//__EdCrossSection_h
