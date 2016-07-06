#include "EdPhysics.h"


#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "TRandom.h"




EdPhysics::EdPhysics(EdModel *model){
    pdg = new TDatabasePDG();
    pdg->ReadPDGTable("eg_pdg_table.txt");
    fRandom = new TRandom2(0);
    SetRandom(fRandom);
    model->SetRandom(fRandom);
    //   gRandom->Delete();
    gRandom = fRandom;
    printf("Seed number %d\n",fRandom->GetSeed());
    target.SetPxPyPzE(0.0, 0.0, 0.0, model->Get_tgMass()); // target at rest
    //    if(model->GetInData().IsQF())spectator.SetXYZM(0,0,0,pdg->GetParticle(model->GetInData().qfspdg);
    mass_model = model->GetMassModel();
    n_part = model->GetNpart();
    nvertex = model->GetNvertex();
    part_pdg[n_part] = pdg->GetParticle(model->GetBeamPID()); // Beam particle stored in part_pdg[n_part]
    double masses2[n_part];
    double width2[n_part];
    for (int i=0; i<n_part; i++) {
      towrite[i] = 1;
      particle_id[i] = model->GetPid(i);
      theta_min[i] = model->GetTheta_min(i);
      theta_max[i] = model->GetTheta_max(i);
      part_pdg[i] = pdg->GetParticle(particle_id[i]); 
      charge[i] = part_pdg[i]->Charge()/3; // Charge is in unit of |e|/3
      masses2[i] = part_pdg[i]->Mass();
      width2[i] = part_pdg[i]->Width();
      if (width2[i] > 0.001) printf("Particle n.%i \t pid=%i \t mass=%.3e GeV width=%.3e Stable(0/1)=%i: Mass will be generated as %s ; theta_min=%.3e theta_max=%.3e \n",i+1,particle_id[i],masses2[i],width2[i],part_pdg[i]->Stable(),model->GetMassModelString(),theta_min[i],theta_max[i]);
      else printf("Particle n.%i \t pid=%i \t mass=%.3e GeV width=%.3e Stable(0/1)=%i; theta_min=%.3e theta_max=%.3e \n",i+1,particle_id[i],masses2[i],width2[i],part_pdg[i]->Stable(),theta_min[i],theta_max[i]);
      
    }
    nvertex = model->GetNvertex();
    int atpart = 0;
    for (int i=0; i<nvertex; i++) {
      overt[i] = model->GetOvert(i);
      if (overt[i] > 0) towrite[overt[i]-1] = 0; // Setting not to write particle already decayed in the LUND or BOS file
      npvert[i] = model->GetNpvert(i);
      v_type[i] = model->GetV_type(i);
      v_ratio[i] = model->GetV_ratio(i);
      if (overt[i] == 0) printf("Vertex n. %i, Origin (Beam + Tg) --> ",i+1);
      else printf("Vertex n. %i, Origin (pid=%i %.3e GeV,  Lifetime=%.3e) --> ",i+1,particle_id[overt[i]-1],part_pdg[overt[i]-1]->Mass(),part_pdg[overt[i]-1]->Lifetime());
      for(int j=0; j<npvert[i] ; j++) {
	p4vector[i][j] = new TLorentzVector();
	masses[i][j] = part_pdg[atpart]->Mass();
	width[i][j] = part_pdg[atpart]->Width();
	//	if (width[i][j] <= 0.001) val_mass[i][j] = masses[i][j];
	val_mass[i][j] = masses[i][j];
	max_mass[i][j] = 0.;
	printf("(pid=%i %.3e GeV) ", particle_id[atpart], part_pdg[atpart]->Mass());
	atpart++;
	if (j==(npvert[i]-1)) printf(" \n");
	else printf(" + ");
      }
      for(int j=0; j<npvert[i] ; j++) {
	for(int k=0; k<npvert[i] ; k++) {
	  if(width[i][k] <= 0.001 && j != k) max_mass[i][j] = max_mass[i][j] - masses[i][k] ;
	}
      }

    }
    if(model->IsQF()){//Also want to write out spectator for quasifree (add at end =npart)
      cout<<"Quasi free so store specator information "<<endl;
      towrite[n_part] = 1;
      particle_id[n_part] = model->GetInput()->GetInData().qfspdg;
      theta_min[n_part] = model->GetTheta_min(0);
      theta_max[n_part] = model->GetTheta_max(TMath::Pi());
      part_pdg[n_part] = pdg->GetParticle(particle_id[n_part]); 
      charge[n_part] = part_pdg[n_part]->Charge()/3; // Charge is in unit of |e|/3
      masses2[n_part] = part_pdg[n_part]->Mass();
      width2[n_part] = part_pdg[n_part]->Width();
    }
    
 

 
 
    return;
}

EdPhysics::~EdPhysics(){
    return;
}

void EdPhysics::MakeEvent(EdOutput *out , EdModel *model){
  // target info

  // Energy of the PHOTON BEAM is sometimes below the threshold of the pi0 mass. The random generator is constant in generation, so the 26-th event is the one with energy below the threshold. This could be also the problem with the multiple particle vertex one, since with multiple particle I require more energy in the final state. I should try to put a single beam with really high energy and check if ken hicks simulation works.

    // double  e_lab = model->GetEnergy();
    // out->SetEin(e_lab);
    // beam.SetPxPyPzE(0.0, 0.0,e_lab,e_lab);
    double tglx = model->GetLx();
    double tgly = model->GetLy();
    double tglength = model->GetLength();

    TVector3 tgtoff = model->GetTgtOffset();


  

    nucl n;
    double A = ((double) (model->Get_tgZ()+model->Get_tgN()));
    double prot_prob = ((double) model->Get_tgZ())/A;
    // Determine which type of nucleon we hit
    if( fRandom->Uniform() < prot_prob ){
	n = kProton;
    } else {
	n = kNeutron;
    }    

    
    //    for (int i=0; i<(n_part+1) ; i++) p4vector[i]->SetPxPyPzE(0.,0.,0.,0.);
    W4vector.SetXYZT(0.,0.,0.,0.);
    Q4vector.SetXYZT(0.,0.,0.,0.);

    double pos_x = GetBeamProfile(tglx);
    double pos_y = GetBeamProfile(tgly);
    double pos_z = fRandom->Uniform(-0.5*tglength,0.5*tglength);
    vertex.SetXYZ(pos_x,pos_y,pos_z);
    vertex = vertex + tgtoff;
    //    printf("vertex at X=%.3e Y=%.3e Z=%.3e \n",pos_x,pos_y,pos_z);
    int test_gen = 0;
    count_phase = 0;
    while (test_gen < nvertex ) test_gen = Gen_Phasespace(model);

    if(model->IsQF())n_part++; //dirty hack to get it to write spectator!

    out->SetEin(e_lab);
    out->SetTheta(theta,n_part);
    out->SetPhi(phi,n_part);
    out->SetEf(Ef,n_part);
    out->Setpf(pf,n_part);
    out->Setpx(px,n_part);
    out->Setpy(py,n_part);
    out->Setpz(pz,n_part);
    out->Setparticle_id(particle_id,n_part);
    out->Setcharge(charge,n_part);
    out->Setweight(weight,n_part);
    out->Settowrite(towrite,n_part);
    
    // Here I need to fix the energy, because I do not know the kinematics at now and I generate randomly the momentum.
    
  
    Z_ion = model->Get_tgZ(); 
    N_ion = model->Get_tgN();
    W = W4vector.M();
    Q2= -Q4vector.M2(); 
    nu= Q4vector.Dot(target)/target.M(); 
    x = Q2/(2*target.M()*nu);
    y = Q4vector.Dot(target)/beam.Dot(target);

    out->SetZ_ion(Z_ion);
    out->SetN_ion(N_ion);
    out->Setx(x);
    out->Sety(y);
    out->SetW(W);
    out->SetQ2(Q2);
    out->Setnu(nu);
    out->Setvx(vx,n_part);
    out->Setvy(vy,n_part);
    out->Setvz(vz,n_part);
 
    //    printf("written vertex at X=%.3e Y=%.3e Z=%.3e \n",vx[0],vy[0],vz[0]);

   if(model->IsQF())n_part--; //dirty hack to get it to write spectator!

    
    out->Write();

    return;
}



double EdPhysics::GetBeamProfile( double sigma){
  return fRandom->Gaus(0.,sigma);
}



TVector3 EdPhysics::Decay_vertex(TLorentzVector *Vp_4, int i, TVector3 vert) {
  TVector3 b_3 ; // beta to boost the LAB frame for going in the pi0 rest frame 
  b_3 = Vp_4->BoostVector(); // return (px/E,py/E,pz/E) (is all in GeV)
  double lifetime = part_pdg[i]->Lifetime();
  TLorentzVector test(0.,0.,0.,TMath::C()*lifetime);
  TVector3 result;

  test.Boost(b_3);
  if (test.Rho() < 0.0001) result = vert;  // Delta vertex for t=lifetime different less than 0.1mm
  else {
    double toptime = lifetime * 20; // exp(-20) = 2.0e-9
    //define vertex of the decayed particles...
    TF1 *fr = new TF1("fr","exp(-x/([0]))",0,toptime) ; // 8.4e-17 is the mean lifetime of the pi0
    fr->SetParameter(0,lifetime);
    double time = fr->GetRandom(0.,toptime);
    TLorentzVector move(0.,0.,0.,TMath::C()*time); // displacement for the creation of the two gammas in the pi0 rest frame (ready to boost (c*t) )
    move.Boost(b_3); // displacement for the creation of the two gammas in the LAB frame
    result.SetX( vert.X() + move.X() );
    result.SetY( vert.Y() + move.Y() );
    result.SetZ( vert.Z() + move.Z() );
    delete fr;
  }
  return result;

}



int EdPhysics::Gen_Mass(int i,EdModel *model) {
  // need to put all values of val_mass[i][j]with this function. The return integer is in case the generation is correct (could be a loop with while ( output < npvert[i] ) In this way I can generate all masses according to the value here generated. Also the order of generation, considering the limit should be random.
  double prob[10];
  fRandom->RndmArray(npvert[i],prob);
  int good_gen = 1;
  int k;
  double total_gen = 0.;
  //  for (int j=0; j<MAX_PART; j++) {
    //  }
  //  printf("Energy = %f \n",e_lab);
  if (overt[i] == 0) { // (Origin Beam + Tg)
    e_lab = model->GetEnergy();
    beam.SetPxPyPzE(0.0, 0.0,e_lab,e_lab);
    //    printf("Set energy from beam \n");
    if(!model->IsQF()) //standard target
      Wtg = beam + target;
    else{  //qf target
      QFTarget(model);
      Wtg = beam + target;
      //cout<<"W "<<Wtg.M()<<" "<<e_lab<<endl;
    }
  }
  else {
    //    TLorentzVector *p4vector_calc = p4vector[i][0]; 
    p4vector_c = new TLorentzVector(*p4vector[i][0]); 
    
    Wtg = *p4vector_c;
    // printf("Vertex %i  particle n. %i mass%.3e \n",i,overt[i]-1,Wtg.M());
  }
  //  printf("Mass at vertex %i part %i = %.3e \n",i,overt[i]-1,Wtg.M());

  for (int j=0; j<npvert[i] ; j++) {
    k = (int)TMath::LocMin(npvert[i],prob);
    val_mass[i][k] = masses[i][k];
    if (width[i][k] > 0.001) {
      val_mass[i][k] = 0.;
      prob[k] = prob[k] +1;
      if (mass_model==1 && (Wtg.M()+max_mass[i][k] - total_gen) > 0.001 ) {
	while (val_mass[i][k] <= 0.001) val_mass[i][k] = fRandom->BreitWigner(masses[i][k],width[i][k]); // Sometimes the random value is outside the limits 
	if ( val_mass[i][k] > (Wtg.M()+max_mass[i][k] - total_gen))  good_gen = 0;
      }
      else if (mass_model==1 && (Wtg.M()+max_mass[i][k] - total_gen) < 0.001 ) good_gen = 0;
      else if (mass_model==2 && (Wtg.M()+max_mass[i][k] - total_gen) > 0.001 ) {
	val_mass[i][k] = fRandom->Uniform(0.001,Wtg.M()+max_mass[i][k] - total_gen); // Sometimes the random val_massm is outside the limits ?!??!?!
      }
      else if (mass_model==2 && (Wtg.M()+max_mass[i][k] - total_gen) < 0.001 ) good_gen = 0; 
      else if (mass_model==3 && (Wtg.M()+max_mass[i][k] - total_gen) > masses[i][k] ) val_mass[i][k] = masses[i][k] ;
      else if (mass_model==3 && (Wtg.M()+max_mass[i][k] - total_gen) > masses[i][k] ) good_gen = 0;
      else if (mass_model==4 && (Wtg.M()+max_mass[i][k] - total_gen) > 0.001 ) {
	while (val_mass[i][k] <= 0.001 || val_mass[i][k] > (Wtg.M()+max_mass[i][k] - total_gen)) val_mass[i][k] = fRandom->BreitWigner(masses[i][k],width[i][k]); // Sometimes the random value is outside the limits 
	       }

      else if (mass_model < 1 || mass_model > 3 ) printf("Mass model %i not allowed: Please check your input file \n",mass_model);
      else good_gen = 0;
      total_gen = total_gen + val_mass[i][k] ; 
    }
  }  // Take away from the mass the stable particle
  //  printf("good_gen = %d \n",good_gen);
  return good_gen; 
}


int EdPhysics::Gen_Phasespace(EdModel *model){

      //    if (valid_event>0) printf("valid events =%i but nvertex=%i",valid_event,nvertex);
  // TLorentzVector *p4vector[n_part+1];
  double weight2;
  double total_mass;
  int atpart = 0;
  int valid_event = 0;
  int good_mass = 0; 
  int failed_event = 0;
  double good_weight = 0.;
  valid_event = 0;
  double tglx = model->GetLx();
  double tgly = model->GetLy();
  double tglength = model->GetLength();

  TVector3 tgtoff = model->GetTgtOffset();

  double pos_x = GetBeamProfile(tglx);
  double pos_y = GetBeamProfile(tgly);
  double pos_z = fRandom->Uniform(-0.5*tglength,0.5*tglength);
  vertex.SetXYZ(pos_x,pos_y,pos_z);
  vertex = vertex + tgtoff;


  for (int i=0; i<nvertex; i++) {
    good_mass=0;
    while (good_mass == 0) good_mass = Gen_Mass(i,model); 
    total_mass = 0.;
    for (int j=0; j<npvert[i]; j++) {
      // val_mass[i][j] = -1.;
      // if (width[i][j]>0.001) {
      // 	while (val_mass[i][j]< 0.001) val_mass[i][j] = Gen_Mass(i,j); // If width > 1MeV, generate the event witha Breit-Wigner probability 
      // }
      // else val_mass[i][j] = masses[i][j];
      total_mass = total_mass + val_mass[i][j];
      //      printf("mass vertex %i particle %i total=%.3e mass=%.3e max_mass%.3e \n",i,j,total_mass,val_mass[i][j],max_mass[i][j]);
    }
    // if (Wtg.M() < total_mass) good_mass = Gen_Mass(i);
    //    printf("mass generated Wtg=%.3e total=%.3e good_mass=%i \n",Wtg.M(),total_mass,good_mass);      
    if (Wtg.M() > total_mass) { // mass check at each vertex
      //   printf("mass generated\n");

      SetDecay(Wtg, npvert[i], val_mass[i]);
      weight2 = Generate();
      good_weight = fRandom->Uniform(1.0);
      if (good_weight <= weight2) {
	weight2 = 1.0;
	valid_event++;
      }
      //      printf("weight= %.3e\n",weight2);
      //   printf("event generated\n");
      for (int j=0; j<npvert[i]; j++) {
	p4vector[i][j+1] = GetDecay(j);
	//	cout << "Particle n." << atpart << " Mass=" << p4vector[atpart]->M() << endl; 
	for (int k=i+1; k<nvertex; k++) {
	  //	  printf ("i=%i k=%i \n",i,k);
	  // if (overt[k] == atpart) p4vector[k][0] = new TLorentzVector(p4vector[i][j+1]->Px(),p4vector[i][j+1]->Py(),p4vector[i][j+1]->Pz(),p4vector[i][j+1]->E());
	  if ( (overt[k]-1) == atpart) p4vector[k][0] = new TLorentzVector(*p4vector[i][j+1]);

	}
	theta[atpart] = p4vector[i][j+1]->Theta();
	phi[atpart] = p4vector[i][j+1]->Phi();
	Ef[atpart] = p4vector[i][j+1]->E();
	pf[atpart] = p4vector[i][j+1]->Rho();
	px[atpart] = p4vector[i][j+1]->Px();
	py[atpart] = p4vector[i][j+1]->Py();
	pz[atpart] = p4vector[i][j+1]->Pz();
	if (atpart<npvert[0] && atpart>0) W4vector += *p4vector[i][j+1]; // I am assuming that the first particle is the scattered beam
	if (atpart == 0)      Q4vector= beam - *p4vector[0][j+1];
	weight[atpart] = weight2;
	//	else weight[atpart] =  weight[overt[i]-1] * weight2;
	if (overt[i] ==0) {
	  vx[atpart] = vertex.X();
	  vy[atpart] = vertex.Y();
	  vz[atpart] = vertex.Z();
	  //	  printf("part %i written vertex at X=%.3e Y=%.3e Z=%.3e \n",atpart,vx[atpart],vy[atpart],vz[atpart]);

	}
	else {
	  if (part_pdg[overt[i]-1]->Width() == 0.0) {
	    //	    printf("Origin particle %i at vertex %i is stable??? vertexes of daughters particles as mother \n", particle_id[overt[i]-1],i); 
	    vx[atpart] = vx[overt[i]-1];
	    vy[atpart] = vy[overt[i]-1];
	    vz[atpart] = vz[overt[i]-1];
	  }
	  else if (part_pdg[overt[i]-1]->Width() > 0.0 && j==0 && i!=0){
	    // cout<<"Going to do decay vertex "<< particle_id[atpart]<<endl;
	    vertex.SetXYZ(vx[overt[i]-1],vy[overt[i]-1],vz[overt[i]-1]);
	    vertex = Decay_vertex(p4vector[i][0],(overt[i]-1),vertex);
	    vx[atpart] = vertex.X() ;
	    vy[atpart] = vertex.Y();
	    vz[atpart] = vertex.Z();
	  }
	  else if (part_pdg[overt[i]-1]->Width() > 0.0 && j>0){
	    vx[atpart] = vertex.X() ;
	    vy[atpart] = vertex.Y();
	    vz[atpart] = vertex.Z();
	  }
	  
	}
	atpart++;
      }
    }
  }
  if(model->GetInput()->IsQF()){
    theta[n_part] = spectator.Theta();
    phi[n_part] = spectator.Phi();
    Ef[n_part] = spectator.E();
    pf[n_part] = spectator.Rho();
    px[n_part] = spectator.Px();
    py[n_part] = spectator.Py();
    pz[n_part] = spectator.Pz();
    vx[n_part] = vertex.X();
    vy[n_part] = vertex.Y();
    vz[n_part] = vertex.Z();
  }
  
  count_phase++;
  if ( (count_phase%100000) == 0) printf("Generated %d events without passing your angle cuts. Could be you want to check your limits\n",count_phase);
  // theta_v_min = TMath::Pi();
  // theta_v_max = 0.;

  // For selecting which particle is going to be created first I will need to use creation time if the space of creation time corresponds to the size where there will still be interaction between the two packets 

  for (int i=0; i<n_part; i++) {
    //    if (towrite[i] == 1) {  // do the angle cut only for particles in the output tht will hit the detector
    if (theta[i] < theta_min[i]) valid_event--; 
    if (theta[i] > theta_max[i]) valid_event--;
      //    }
  }

  return valid_event;


}
void EdPhysics::QFTarget(EdModel *model){
  Float_t ptar, Etar, Espec, costhtar, thtar, phtar, pxtar, pytar, pztar,smass;
  //Get Fermi momentum sampled from distribution given in input file
  ptar = model->GetFermi()->GetRandom()/1000.;
  //Get random flat cos theta
  costhtar   = fRandom->Uniform( -1., 1. );
  thtar      = TMath::ACos( costhtar );
  //Get random flat phi
  phtar      = fRandom->Uniform( -TMath::Pi(), TMath::Pi() );
  //calculate momentum components
  pxtar      = ptar * TMath::Sin( thtar ) * TMath::Cos( phtar );
  pytar      = ptar * TMath::Sin( thtar ) * TMath::Sin( phtar );
  pztar      = ptar * TMath::Cos( thtar );
  
  // Force spectator on mass shell
  //get its mass from input data and PDG database
  smass=pdg->GetParticle(model->GetInput()->GetInData().qfspdg)->Mass();
  Espec = TMath::Sqrt(ptar*ptar + smass * smass);
  //set spectator 4 momentum (opposite p to quasi target)
  spectator.SetXYZT(-pxtar,-pytar,-pztar, Espec);
  
  //calculate quasi target energy from energy conservation
  Etar  = model->Get_tgMass() - Espec;
  target.SetXYZT(pxtar,pytar,pztar, Etar);
  
}

// double EdPhysics::t_reaction(TLorentzVector *Vrecoil_tg_4 ) {
//   TLorentzVector tg_v4;
//   tg_v4.SetPxPyPzE(Vrecoil_tg_4->Px(),Vrecoil_tg_4->Py(),Vrecoil_tg_4->Pz(),Vrecoil_tg_4->E());
//   TLorentzVector t_v4;
//   t_v4 = target - tg_v4; // Operation allowed just between TLorentzVector and not TLorentzVector*
//   return t_v4.M2();

// }


// Double_t TPart_ident_ct::theta_pip_mc() {

//   p_theta_pip_mc = 0;

//   TVector3 V3_e, V3_e2;
//   TVector3 V3_1 , V3_2 , V3_3, V3_4;

//   V3_e.SetXYZ(0,0,Ebeam); // e beam momentum
//   V3_e2.SetXYZ(px_mc(11),py_mc(11),pz_mc(11)); // Scattered electron  
//   V3_1.SetXYZ(px_mc(-211),py_mc(-211),pz_mc(-211)); // pi- 3-momentum
//   V3_2.SetXYZ(px_mc(211),py_mc(211),pz_mc(211)); // pi+ 3-momentum
//   V3_4 = ( V3_e - V3_e2 ) - ( V3_1 + V3_2) ; // outgoing proton
//      // 3-momentum = - t =
//      // gamma* - rho0

//   TLorentzVector V4_1(V3_1,sqrt(V3_1.Mag2() + pow(0.13957,2))); // pi- quadrimomentum
//   TLorentzVector V4_2(V3_2,sqrt(V3_2.Mag2() + pow(0.13957,2))); // pi+ quadrimomentum
//   TLorentzVector V4_4(V3_4,sqrt(V3_4.Mag2() + pow(0.938272,2))); // outgoing proton quadrimomentum


//   TLorentzVector V4_3 ; // rho momentum
//   V4_3 = V4_1 + V4_2; 
//   V3_3 = V3_1 + V3_2;

//   TVector3 b_3 ; // beta to boost the LAB frame for going in the rho0
//  // rest frame 

//   b_3 = V4_3.BoostVector(); // return (px/E,py/E,pz/E) (is all in GeV)

//   b_3 = - b_3; // the Boost function in TLorentzVector is an anti-lorentz transform 

//   V4_4.Boost(b_3); // proton in the rho0 frame

//   V4_2.Boost(b_3); // pi+ in the rho0 frame

//   V3_4.SetXYZ(V4_4.Px(),V4_4.Py(),V4_4.Pz()) ; // proton momentum in
//        // the rho0 frame
//   V3_2.SetXYZ(V4_2.Px(),V4_2.Py(),V4_2.Pz()) ; // pi+ momentum in
//        // the rho0 frame
  
//   p_theta_pip_mc = TMath::Pi() - V3_2.Angle(V3_4) ;  // angle between
//      // the opposite 
//      // of the
//      // direction of
//      // the proton and
//      // the pi+
  
//   return p_theta_pip_mc;
  
  
  


// }

// //-----------------------------------------
// double EdPhysics::phi_lep_had(TLorentzVector *Vrecoil_tg_4, TLorentzVector *Velectron_4, TLorentzVector *Vmeson_4, double e_lab) {

// //   // 11 = electron , -211 = pi- , 211 = pi+

//   double phi_lep_had_v = 0;

//   TVector3 V3_1 , V3_2 ,V3_3 ,V3_6;
//   V3_1 = Velectron_4->Vect(); // e' momentum
//   V3_2.SetXYZ(0,0,e_lab); // e momentum
//   V3_3 = V3_2 - V3_1;   // Virtual photon momentum vector
//   V3_6 = Vmeson_4->Vect(); // rho0 momentum


//   TVector3 n_lep, n_had;

//   n_lep = V3_2.Cross(V3_1);

//   n_had = V3_3.Cross(V3_6);

  
//   phi_lep_had_v = n_lep.Angle(n_had);

//   return phi_lep_had_v;

// }


// double EdPhysics::phi_pip_mc() {

//    p_phi_pip_mc = 0;
  
//   TVector3 V3_e, V3_e2, V3_g;
//   TVector3 V3_1 , V3_2 , V3_3, V3_4;

//   V3_e.SetXYZ(0,0,Ebeam); // e beam momentum
//   V3_e2.SetXYZ(px_mc(11),py_mc(11),pz_mc(11)); // Scattered electron  
//   V3_1.SetXYZ(px_mc(-211),py_mc(-211),pz_mc(-211)); // pi- 3-momentum
//   V3_2.SetXYZ(px_mc(211),py_mc(211),pz_mc(211)); // pi+ 3-momentum
//   V3_4 = ( V3_e - V3_e2 ) - ( V3_1 + V3_2) ; // outgoing proton
//      // 3-momentum = - t =
//      // gamma* - rho0
//   V3_g = V3_e - V3_e2;  // gamma* momentum

//   TLorentzVector V4_1(V3_1,sqrt(V3_1.Mag2() + pow(0.13957,2))); // pi- quadrimomentum
//   TLorentzVector V4_2(V3_2,sqrt(V3_2.Mag2() + pow(0.13957,2))); // pi+ quadrimomentum
//   TLorentzVector V4_4(V3_4,sqrt(V3_4.Mag2() + pow(0.938272,2))); // outgoing proton quadrimomentum


//   TLorentzVector V4_3 ; // rho momentum
//   V4_3 = V4_1 + V4_2; 
//   V3_3 = V3_1 + V3_2;


//   TVector3 n_had, n_z, n_x , opip;

//   n_had = V3_g.Cross(V3_3);  // p_gamma* X p_rho    : is parallel to
//      //                       the y axis and
//      //                       is perpendicular
//      //                       to the boost
//                              // it is let unchanged from the boost

//   n_had = n_had.Unit();      // creatin the unit vector in the
//      // y direction 

//   TVector3 b_3 ; // beta to boost the LAB frame for going in the rho0
//  // rest frame 

//   b_3 = V4_3.BoostVector(); // return (px/E,py/E,pz/E) (is all in GeV)

//   b_3 = - b_3;

//   V4_4.Boost(b_3); // proton in the rho0 frame

//   V4_2.Boost(b_3); // pi+ in the rho0 frame

//   V3_4.SetXYZ(V4_4.Px(),V4_4.Py(),V4_4.Pz()) ; // proton momentum in
//        // the rho0 frame
//   V3_2.SetXYZ(V4_2.Px(),V4_2.Py(),V4_2.Pz()) ; // pi+ momentum in
//        // the rho0 frame
  
//   n_z = - V3_4.Unit();  // z axis is opposite to the outgoing proton
//   n_x = n_had.Cross(n_z); // x axis = (y axis) X (z axis)  

//   opip = n_z.Cross(V3_2);

//   Double_t sen_phi_ev = asin( - n_x.Dot(opip) / opip.Mag());

//   p_phi_pip_mc = acos(n_had.Dot(opip) / opip.Mag()); // this angle is
//      // defined from 0
//      // to 2 Pi
//   if (sen_phi_ev < 0 ) {
//     p_phi_pip_mc = 2 * TMath::Pi() - p_phi_pip_mc ;  // if sen (phi)<0
//      // phi= 2Pi - phi
     
//   }
   
  
//   return p_phi_pip_mc;
  
  


// }
