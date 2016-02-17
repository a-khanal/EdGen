#include "EdCrossSection.h"

/*
 * t cannot be 0 ; tmax < t < tmin < 0
 *  only photoproduction Q^2 = 0 = M1^2
 */


EdCrossSection::EdCrossSection(EdModel *model){

  // Masses in GeV
  MP= 0.938272046;
  MPI = 0.1349766;
  EPS = 0.0000;
  mass[0] = 0.0;
  mass[1] = 0.0;
  mass[2] = MP;
  mass[3] = MPI;
  mass[4] = MP;
  double plab;
  double E, z, theta;

  plab = 1; // Units = GeV
  E = sqrt( MP*MP + 2*MP*plab ); // Ecm = 3.2031 GeV corresponds to Plab = 9 GeV
  printf("Plab   E      Cos     t        Dsig/Dt    Dsig/DOmega\n");
  printf("GeV    GeV    -       GeV^2    mb/Gev^2     mb\n");
  for( z = -1; z <= 1; z = z + 0.1 ){
    printf("%.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n",plab, E, z,Cos2T(E,z),
	   Pi0PhotCS_S(E,acos(z))/1000, Pi0PhotCS_S(E,acos(z))/Coef(E)/1000 );
  }

}

// *********************************************************************************
void EdCrossSection::MAID_Pi0(double s, double complex MltPole[5][10]){
	/* The MAID 2007 solution computed for the reaction
	 *  gamma p --> pi0 p
	 *  Other reactions are related to isospin
	 *
	 * return double complex MltPole[5][10]
	 * MltPole[1][L] is E_{L+}
	 * MltPole[2][L] is E_{L-}
	 * MltPole[3][L] is M_{L+}
	 * MltPole[4][L] is M_{L-}
	 *
	 * MAXIMUM L is 9
	 */


	int L,i,k;
	// solution valid up to W = 2 GeV
	if(sqrt(s)>2.0){
		for(L=0;L<11;L++){
		for(i=0;i<5 ;i++){
			MltPole[i][L] = 0.0;
		}}
		return ;
	}

	int MAXBIN = 185 ;
	double wval[800]={0}; double wlab[800]={0};	double wcm[800]={0};

	double ReL0[10][1000]={{0},{0}}, ImL0[10][1000]={{0},{0}}; // first index is L ;  second is W
	double ReLP[10][1000]={{0},{0}}, ImLP[10][1000]={{0},{0}};
	double ReLM[10][1000]={{0},{0}}, ImLM[10][1000]={{0},{0}};
	double complex LPi0[10][1000]={{0},{0}};

	double complex Ep[10][186] = {{0},{0}}, Em[10][186] = {{0},{0}};
	double complex Mp[10][186] = {{0},{0}}, Mm[10][186] = {{0},{0}};
	double complex EpPi0[10] = {0}, EmPi0[10] = {0};
	double complex MpPi0[10] = {0}, MmPi0[10] = {0};

	FILE *pi0maidL[6];	// Load MAID multipoles
		pi0maidL[0] = fopen("MAID/pi0maidL0.txt","r");
		pi0maidL[1] = fopen("MAID/pi0maidL1.txt","r");
		pi0maidL[2] = fopen("MAID/pi0maidL2.txt","r");
		pi0maidL[3] = fopen("MAID/pi0maidL3.txt","r");
		pi0maidL[4] = fopen("MAID/pi0maidL4.txt","r");
		pi0maidL[5] = fopen("MAID/pi0maidL5.txt","r");

	/*    W          EL+(0)           EL+(+)           EL+(-)       w(lab)  w(cm)
 	 	 (MeV)      Re      Im       Re      Im       Re      Im      (MeV)   (MeV)*/

	for( i =1 ; i<4*MAXBIN+1 ; i++){
	for( L = 0 ; L < 6 ; L++){
		fscanf(pi0maidL[L],"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&wval[i],
				&ReL0[L][i],&ImL0[L][i],&ReLP[L][i],&ImLP[L][i],
				&ReLM[L][i],&ImLM[L][i],&wlab[i],&wcm[i]);

		// Neutral Pion = (0) + (+)
		// CHANGE HERE for other reaction
		LPi0[L][i] = ReL0[L][i] + I*ImL0[L][i] + ReLP[L][i] + I*ImLP[L][i] ;
	}}

	fclose(pi0maidL[0]);fclose(pi0maidL[1]);fclose(pi0maidL[2]);
	fclose(pi0maidL[3]);fclose(pi0maidL[4]);fclose(pi0maidL[5]);

	for( i=1; i<MAXBIN+1;i++){
	for(L = 0 ; L < 6 ; L++){
		k = i + 0*MAXBIN;  	Ep[L][i] = LPi0[L][k];
		k = i + 1*MAXBIN; 	Em[L][i] = LPi0[L][k];
		k = i + 2*MAXBIN;	Mp[L][i] = LPi0[L][k];
		k = i + 3*MAXBIN;	Mm[L][i] = LPi0[L][k];
	}}

	double w = sqrt(s);
	int ind = floor((w - 1.08)/0.005 + 1);
	double eps = (w-wval[ind]/1000)/0.005;

	for(L = 0 ; L < 10 ; L++){
		// Linear extrapolation between W bins
		EpPi0[L] = (1-eps)*Ep[L][ind] + eps*Ep[L][ind+1];
		EmPi0[L] = (1-eps)*Em[L][ind] + eps*Em[L][ind+1];
		MpPi0[L] = (1-eps)*Mp[L][ind] + eps*Mp[L][ind+1];
		MmPi0[L] = (1-eps)*Mm[L][ind] + eps*Mm[L][ind+1];

		// store multipoles for the output
		MltPole[1][L] = EpPi0[L];
		MltPole[2][L] = EmPi0[L];
		MltPole[3][L] = MpPi0[L];
		MltPole[4][L] = MmPi0[L];
	}
	return;
}

// *********************************************************************************
void EdCrossSection::MltPole2F(double s, double t, double complex MltPole[5][10] ,double complex CGLNF[5]){
	/*
	 * MltPole[1][L] is E_{L+}
	 * MltPole[2][L] is E_{L-}
	 * MltPole[3][L] is M_{L+}
	 * MltPole[4][L] is M_{L-}
	 *
	 * return CGLN Fi in mFm
	 */
	int L,k;

	// kinematics
	double w = sqrt(s);
	double som = 2*MP*MP + MPI*MPI;
	double u = -s - t + som;
	double kcm = (s - MP*MP)/(2.0*w);
	double Ex = (s + MP*MP - MPI*MPI)/(2*w);
	double qcm = sqrt(Ex*Ex - MP*MP);
	double z = ( s *(t-u) - MP*MP*(MPI*MPI - MP*MP))/(4*s*kcm*qcm);	// z = cos theta
	//double tmin = pow(MPI,4)/4/s - pow(qcm - kcm,2);
	//double tmax = tmin - 4*kcm*qcm;

	double pol1p[10] = {0}; // First derivative of Legendre polynomials
	double pol2p[10] = {0}; // second derivative of Legendre polynomials
	pol1p[1] = 1.0;
	pol1p[2] = 3*z;
	pol1p[3] = (-3 + 15*z*z)/2.0;
	pol1p[4] = (-60*z + 140*z*z*z)/8;
	pol1p[5] = (15 - 210*z*z + 315*z*z*z*z)/8;
	pol1p[6] = (210*z - 1260*z*z*z + 1386*pow(z,5))/16;
	pol1p[7] = (-35 + 945*z*z - 3465*pow(z,4) + 3003*pow(z,6))/16;
	pol1p[8] = (-2520*z + 27720*z*z*z - 72072*pow(z,5) + 51480*pow(z,7))/128;
	pol1p[9] = (315 - 13860*z*z + 90090*pow(z,4) - 180180*pow(z,6) + 109395*pow(z,8))/128;

	pol2p[2] = 3;
	pol2p[3] = 15*z;
	pol2p[4] = (-60+420*z*z)/8;
	pol2p[5] = (-420*z + 1260*z*z*z)/8.0;
	pol2p[6] = (210 - 3780*z*z + 6930*pow(z,4))/16.0;
	pol2p[7] = (1890*z - 13860*z*z*z + 18018*pow(z,5))/16.0;
	pol2p[8] = (-2520 + 83160*z*z - 360360*pow(z,4) + 360360*pow(z,6) )/128;
	pol2p[9] = ((-27720*z + 360360*z*z*z - 1081080*pow(z,5) + 875160*pow(z,7)))/128;

	CGLNF[1] = 0; CGLNF[2] = 0; CGLNF[3] = 0; CGLNF[4] = 0;
	for( L=0; L<10;L++){
		k = L-1; if(L==0){k = 0;}
		/*
		* MltPole[1][L] is E_{L+}
		* MltPole[2][L] is E_{L-}
	    * MltPole[3][L] is M_{L+}
		* MltPole[4][L] is M_{L-}
		*/
		CGLNF[1] = CGLNF[1] + (L*MltPole[3][L] + MltPole[1][L])*pol1p[L+1];
		CGLNF[1] = CGLNF[1] + ( (L+1)*MltPole[4][L] + MltPole[2][L] )*pol1p[k];
		CGLNF[2] = CGLNF[2] + ( (L+1)*MltPole[3][L] +L*MltPole[4][L])*pol1p[L];
		CGLNF[3] = CGLNF[3] + (MltPole[1][L] - MltPole[3][L])*pol2p[L+1];
		CGLNF[3] = CGLNF[3] + (MltPole[2][L] + MltPole[4][L])*pol2p[L-1];
		CGLNF[4] = CGLNF[4] + (MltPole[3][L] - MltPole[1][L] - MltPole[4][L] - MltPole[2][L])*pol2p[L];
	}

	double factor = 1.46188; // hbar*c
	for(k=1;k<5;k++){
	CGLNF[k] = CGLNF[k] * factor;
	}

	return;
}

// *********************************************************************************
double EdCrossSection::Pi0PhotCS_S(double E,double theta){
/*
 * compute the differential cross section from s-channel helicities
 * in micro barns/Gev^2
 */
	double pa[4], pb[4], pc[4], pd[4];
	int hel[4] = {2,1,0,1};						// hel = {1,+,0,+} (x2)
	double complex  res[4] ;
	double sig;
	kin2to2(E, theta, mass, pa, pb, pc, pd);

	res[0] = Pi0PhotAmpS(pa, pb, pc, hel);		// hel = {1,+,0,+} (x2)

	hel[1] = -1;								// hel = {1,-,0,+} (x2)
	res[1] = Pi0PhotAmpS(pa, pb, pc, hel);

	hel[3] = -1;								// hel = {1,-,0,-} (x2)
	res[2] = Pi0PhotAmpS(pa, pb, pc, hel);

	hel[1] = +1;								// hel = {1,+,0,-} (x2)
	res[3] = Pi0PhotAmpS(pa, pb, pc, hel);

	sig = pow(cabs(res[0]),2) + pow(cabs(res[1]),2) + pow(cabs(res[2]),2) + pow(cabs(res[3]),2);
	sig = sig/32/M_PI/(E*E-MP*MP)/(E*E-MP*MP) * 389.3;
	return sig;
}

// *********************************************************************************
double complex EdCrossSection::Pi0PhotAmpS(double pa[4],double pb[4],double pc[4], int hel[4]){
	/* Regge amplitudes for gamma + p --> pi0  + p
	 * See arXiv:1505.02321 for the formulas
	 * Amplitudes in the S-CHANNEL !
	 * Helicities are defined in the c.o.m. of (gamma,p)
	 * Inputs:
	 * 		pa, pb, pc 					: four vectors in the c.o.m. frame.
	 * 		hel={mua, mub, muc, mud} 	: 2 x helicities of a,b,c,d
	 * Outputs:
	 * 		ampl						: the amplitude ; complex number
	 */
	double complex amp = 0.0;
	double complex Ai[5] = {0}, FiReg[5] = {0}, FiMAID[5] = {0}, Fi[5] = {0};
	double pab[4], pbc[4], pca[4], pd[4];
	struct Kin var;
	double mass2[5]={0};
	double complex MltPole[5][10];
	int i;
	// Only valid for Q^2 = 0 --> photon helicity should be +1 or -1
	if ( abs(hel[0]) != 2 ) return 0.0;

	for(i=0;i<4;i++){							// compute the Mandelstam invariant
		pab[i] = pa[i] + pb[i];						// from four-vectors
		pbc[i] = pb[i] - pc[i];						//
		pca[i] = pc[i] - pa[i];						//
		pd[i]  = pa[i] + pb[i] - pc[i];				//
	}
	var.s = snorm( pab ) + I*EPS;
	var.t = snorm( pca ) - I*EPS;
	// photon should be real
	mass2[1]   = 0.0;
	mass2[2]  = sqrt( snorm( pb ) );
	mass2[3]  = sqrt( snorm( pc ) );
	mass2[4]  = sqrt( snorm( pd ) );

	kinematics(var.s,var.t, mass2, &var);	// fill the structure var with all kin. quantities

	// Compute the CGLN Fi for a given model
	// CHANGE HERE TO CHANGE THE MODEL
	// Regge for Ecm > 2 GeV
	if(creal(var.s)>4){
		Regge_CGLNAi(var.s,var.t,Ai);				// call the Regge model
		CGLNA2F(var, Ai, FiReg);				// Regge model: Convert CGLN Ai to CGLN Fi
	}
	// MAID for Ecm < 2 GeV
	if(creal(var.s)<4){
		MAID_Pi0(creal(var.s),  MltPole);		// MAID multipoles for pi0 photoproduction
		MltPole2F(creal(var.s), creal(var.t), MltPole,FiMAID);	// MAID CGLN Fi for pi0 photoproduction
	}

	// Combine the two models
	for( i=0;i<5;i++){
		Fi[i] = FiMAID[i] + FiReg[i];
	}

	// Test nucleon helicities hel[1] and hel[3]
	if ( hel[1] == 1 && hel[3] == 1 ) {
		// hel = {1,+,0,+}
		amp = sqrt(2) * var.sinsh * (Fi[2] + Fi[1] )
			  + 1/sqrt(2) * var.sins * var.cossh * (Fi[3] + Fi[4]);
	}
	else if  ( hel[1] == -1 && hel[3] == -1 ){
		// hel = {1,-,0,-}
		amp = -1/sqrt(2) * var.sins * var.cossh * (Fi[3] + Fi[4]);
	}
	else if ( hel[1] == -1 && hel[3] == 1 ){
		// hel = {1,-,0,+}
		amp = sqrt(2) * var.cossh * (Fi[2] - Fi[1] )
			  + 1/sqrt(2) * var.sins * var.sinsh * (Fi[3] - Fi[4]);
	}
	else if ( hel[1] == 1 && hel[3] == -1 ){
		// hel = {1,+,0,-}
		amp = 1/sqrt(2) * var.sins * var.sinsh * (Fi[3] - Fi[4]);
	}
	else return 0.0 ;

	// if negative photon helicity, nucleon flip amplitude changes sign
	if (hel[0] == -2 && hel[1] != hel[3] ) amp = -1*amp;

	// there is a factor 8*PI*W between Hi and helicity amplitudes
	return amp * 8 * M_PI * csqrt(var.s);
}

// *********************************************************************************
void EdCrossSection::Regge_CGLNAi(double complex s, double complex t, double complex CGLNA[5]){
	/*
	 * Model for gamma p --> pi0 p from arXiv:1505.02321
	 * Vincent Mathieu June 2015
	 */

	// Model parameters
	double alpV[3] = { 0.442457, 1.099910, 0.0};		// vector trajectory
	double alpC[3] = { 0.461811, 0.166543, 0.0};		// vector cut trajectory
	double alpA[3] = {-0.193332, 1.021300, 0.0};		// axial-vector trajectory
	double g1 = -3.8873, g4 = 10.1403, g1c = 1.76187, g4c = 3.58089, g2 = 8.887; // residues

	double complex Rv, Rc, Ra;
	double complex avec, acut, aaxi;
	// trajectories:
	avec = alpV[0] + t*alpV[1] + t*t*alpV[2];
	acut = alpC[0] + t*alpC[1] + t*t*alpC[2];
	aaxi = alpA[0] + t*alpA[1] + t*t*alpA[2];

	// Regge factors:
	Rv = cgamma( 1.0 - 1.0*avec, 0)/2 * ( 1.0-1.0*cexp(-1.0*I*M_PI*avec) ) * cpow(s,avec-1);
	Rc = cgamma( 1.0 - 1.0*acut, 0)/2 * ( 1.0-1.0*cexp(-1.0*I*M_PI*acut) ) * cpow(s,acut-1);
	Ra = cgamma( 1.0 - 1.0*aaxi, 0)/2 * ( 1.0-1.0*cexp(-1.0*I*M_PI*aaxi) ) * cpow(s,aaxi-1);
	Rc = Rc / clog(s);
	// IF ONLY VECTOR POLE
	//Rc = 0; Ra = 0;

	// Scalar amplitudes:
	CGLNA[1] = t* ( g1*Rv + g1c*Rc);
	CGLNA[2] = g2 * Ra - ( g1*Rv + g1c*Rc);
	CGLNA[3] = 0;
	CGLNA[4] = g4*Rv + g4c*Rc;

	return ;

}

// *********************************************************************************
void EdCrossSection::CGLNA2F(struct Kin var, double complex Ai[5], double complex Fi[5]){
	/*
	 * Compute CGLN Fi(s,t) from CGLN Ai(s,t)
	 */
	double complex w;
	double complex E2p, E2m, E4p, E4m;

	w = csqrt(var.s);
	E2p = var.E2s + var.m2 ;
	E2m = var.E2s - var.m2 ;
	E4p = var.E4s + var.m4 ;
	E4m = var.E4s - var.m4 ;

	Fi[1] = (w - var.m2) * Ai[1] + (var.m3*var.m3 - var.t)/2 * (Ai[3] - Ai[4]);
	Fi[1] = Fi[1]  + ( w - var.m2 ) * ( w - var.m4 ) * Ai[4];
	Fi[1] = Fi[1] * csqrt( E2p * E4p ) / ( 8 * M_PI * w);

	Fi[2] = -1.0*(w + var.m2) * Ai[1] + (var.m3*var.m3 - var.t)/2 * (Ai[3] - Ai[4]);
	Fi[2] = Fi[2]  + ( w + var.m2 ) * ( w + var.m4 ) * Ai[4];
	Fi[2] = Fi[2] * csqrt( E2m * E4m ) / ( 8 * M_PI * w);

	Fi[3] = (w + var.m2) * (  (w - var.m2)*Ai[2] + Ai[3] - Ai[4] );
	Fi[3] = Fi[3] * csqrt( E2m * E4p) * var.qs  / ( 8 * M_PI * w);

	Fi[4] = (w - var.m2) * ( -1.0*(w + var.m2)*Ai[2] + Ai[3] - Ai[4] );
	Fi[4] = Fi[4] * csqrt( E2p * E4m) * var.qs  / ( 8 * M_PI * w);

	return;
}

// *********************************************************************************
void EdCrossSection::kin2to2(double Ecm, double theta, double mass2[5], double pa[4],double pb[4],double pc[4], double pd[4]){
/*
 * Kinematics of the 2-to-2 reaction, a + b --> c + d
 * Inputs:
 * 		Ecm center of mass energy
 * 		cos cosine of the scattering angle in the center of mass
 * 		mass2 = {ma,mb,mc,md} vector with the masses of external particles
 * 	Outputs:
 * 		pa, pb, pc, pd are the momenta of the particles
 */
	double Ea, Eb, Ec, Ed; 			// Energies of the particles
	double pi, pf;					// initial and final breakup momenta
	double ma, mb, mc, md;			// masses of the particles

	ma = mass2[1]; mb = mass2[2]; mc = mass2[3]; md = mass2[4];

	// Check that inputs are valid
	if( ( ma<0 ) || ( mb<0 ) || ( mc<0 ) || ( md<0 )  ) {
		printf("\n*** Wrong masses in kin2to2 ! *** \n\n");
		return ;
	}
	if( ( Ecm < ma + mb ) || ( Ecm < mc + md )  ) {
		printf("\n*** Wrong total energy in kin2to2 ! *** \n\n");
		return ;
	}

	Ea = ( Ecm*Ecm + ma*ma - mb*mb )/2/Ecm;
	Eb = ( Ecm*Ecm - ma*ma + mb*mb )/2/Ecm;
	Ec = ( Ecm*Ecm + mc*mc - md*md )/2/Ecm;
	Ed = ( Ecm*Ecm - mc*mc + md*md )/2/Ecm;

	pi = sqrt(Ea*Ea - ma*ma);
	pa[0] = Ea; pa[1] = 0; pa[2] = 0; pa[3] = +pi;
	pb[0] = Eb; pb[1] = 0; pb[2] = 0; pb[3] = -pi;

	pf = sqrt(Ec*Ec - mc*mc);
	pc[0] = Ec; pc[1] = +pf*sin(theta); pc[2] = 0; pc[3] = +pf*cos(theta);
	pd[0] = Ed; pd[1] = -pf*sin(theta); pd[2] = 0; pd[3] = -pf*cos(theta);

	return;
}

// *********************************************************************************
void EdCrossSection::kinematics(double complex s, double complex t, double mass2[5], struct Kin *var)
{
	double m12, m22, m32, m42; 	// masses squared
	double complex t0, t1, u ;

	var->s = s;
	var->t = t;

	var->m1 = mass2[1];
	var->m2 = mass2[2];
	var->m3 = mass2[3];
	var->m4 = mass2[4];

	m12 = mass2[1] * mass2[1] ;
	m22 = mass2[2] * mass2[2] ;
	m32 = mass2[3] * mass2[3] ;
	m42 = mass2[4] * mass2[4] ;

	var->ks = csqrt( lambda(s, m12, m22) / 4 / s );
	var->qs = csqrt( lambda(s, m32, m42) / 4 / s );
	var->kt = csqrt( lambda(t, m12, m32) / 4 / t );
	var->pt = csqrt( lambda(t, m22, m42) / 4 / t );

	// nuclei energies in s- and t-channel frames
	var->E2s = ( s + m22 - m12 ) / 2 / csqrt(s);
	var->E4s = ( s + m42 - m32 ) / 2 / csqrt(s);
	var->E2t = ( t + m22 - m42 ) / 2 / csqrt(t);
	var->E4t = ( t + m42 - m22 ) / 2 / csqrt(t);

	// t1 = tmin ; t0 = tmax ; t0 < t < t1
	t1 = cpow(m12 - m32 - m22 + m42,2)/(4*s) - (var->ks - var->qs) * (var->ks - var->qs);
	t0 = t1 - 4 * var->ks *var->qs ;
	u = -1.0* s - 1.0*t + m12 + m22 + m32 + m42 ;	// Mandelstam s variable
	var->phi = s * (t - t1) * (t0 - t) ;	// Kibble function

	var->coss = 1 + (t - t1)/(2 * var->qs * var->ks) ;
	//var->cost = (t*(s-u) + (m12-m32) * (m22-m42) )/(4 * t * var->pt * var->kt) ;
	var->cost = (t*(s-u) + (m12-m32) * (m22-m42) )/(csqrt(lambda(var->t,m12,m32) * lambda(var->t,m22,m42) ) ) ;
	var->sins = csqrt( var->phi / s ) / ( 2 * var->qs * var->ks) ;
	var->sint = csqrt( var->phi / t ) / ( 2 * var->kt * var->pt) ;
	var->cossh = csqrt( (1 + var->coss ) / 2 );
	var->sinsh = csqrt( (1 - var->coss ) / 2 );
	var->costh = csqrt( (1 + var->cost ) / 2 );
	var->sinth = csqrt( (1 - var->cost ) / 2 );

	return ;
}

// *********************************************************************************
double EdCrossSection::Cos2T(double E, double z){
	/*
	 * Compute t from Ecm = W and cos = z
	 */
	double t;
	double s = E*E;
	double ma, mb, mc, md;			// masses of the particles
	ma = mass[1]; mb = mass[2]; mc = mass[3]; md = mass[4];

	double qi, qf;
	qi = sqrt(lambda_real(s, ma*ma, mb*mb) /4/s);
	qf = sqrt(lambda_real(s, mc*mc, md*md) /4/s);
	double som = ma*ma + mb*mb + mc*mc + md*md;

	t = 2*qi*qf*z + som/2 - s/2 - (ma*ma - mb*mb) * (mc*mc - md*md) /2/s ;
	return t;
}

// *********************************************************************************
double EdCrossSection::Coef(double E){
	/*
	 * dSig/dt = Coef dSig/dOmega
	 * Coef = Pi/qi/qf
	 */
	double s = E*E;
	double ma, mb, mc, md;			// masses of the particles
	ma = mass[1]; mb = mass[2]; mc = mass[3]; md = mass[4];

	double qi, qf;
	qi = sqrt(lambda_real(s, ma*ma, mb*mb) /4/s);
	qf = sqrt(lambda_real(s, mc*mc, md*md) /4/s);
	return M_PI/qi/qf;
}

// *********************************************************************************
double complex EdCrossSection::lambda(double complex a, double b, double c){
// triangle function
	return a*a + b*b + c*c - 2*(a*b + b*c + c*a);
}

// *********************************************************************************
double EdCrossSection::lambda_real(double a, double b, double c){
// triangle function
	return a*a + b*b + c*c - 2*(a*b + b*c + c*a);
}


// *********************************************************************************
double EdCrossSection::snorm(double p[5]){
// Norm squared of a quadri-vector ; p[0] is the energy ; p[1-3] are x,y,z components
	return p[0]*p[0] - ( p[1]*p[1] + p[2]*p[2] + p[3]*p[3] );
}

// *********************************************************************************
double complex EdCrossSection::cgamma(double complex z,int OPT)
{
    double complex g, infini= 1e308+ 0*I; // z0,z1
    double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
    double na,t,x1,y1,sr,si;
    int j,k;

    static double a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};

    x = creal(z);
    y = cimag(z);
    if (x > 171) return infini;
    if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
        return infini;
    else if (x < 0.0) {
        x1 = x;
        y1 = y;
        x = -x;
        y = -y;
    }
    x0 = x;
    if (x <= 7.0) {
        na = (int)(7.0-x);
        x0 = x+na;
    }
    q1 = sqrt(x0*x0+y*y);
    th = atan(y/x0);
    gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
    gi = th*(x0-0.5)+y*log(q1)-y;
    for (k=0;k<10;k++){
        t = pow(q1,-1.0-2.0*k);
        gr += (a[k]*t*cos((2.0*k+1.0)*th));
        gi -= (a[k]*t*sin((2.0*k+1.0)*th));
    }
    if (x <= 7.0) {
        gr1 = 0.0;
        gi1 = 0.0;
        for (j=0;j<na;j++) {
            gr1 += (0.5*log((x+j)*(x+j)+y*y));
            gi1 += atan(y/(x+j));
        }
        gr -= gr1;
        gi -= gi1;
    }
    if (x1 <= 0.0) {
        q1 = sqrt(x*x+y*y);
        th1 = atan(y/x);
        sr = -sin(M_PI*x)*cosh(M_PI*y);
        si = -cos(M_PI*x)*sinh(M_PI*y);
        q2 = sqrt(sr*sr+si*si);
        th2 = atan(si/sr);
        if (sr < 0.0) th2 += M_PI;
        gr = log(M_PI/(q1*q2))-gr;
        gi = -th1-th2-gi;
        x = x1;
        y = y1;
    }
    if (OPT == 0) {
        g0 = exp(gr);
        gr = g0*cos(gi);
        gi = g0*sin(gi);
    }
    g = gr + I*gi;
    return g;
}


