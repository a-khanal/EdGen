EdGen Event Generator
==================

This code is for a common Event Generator for the HASPECT collaboration.
It uses the ROOT (root.cern.ch) PhaseSpace generator has basis and 
See http://www2.ph.ed.ac.uk/~lzana/Documents/zana_haspect_genova2014.pdf for a quick presentation on this code.
This version respect to the common version it supports the BOS output format (for CLAS6 analysis) if the system has it.

Prerequisites
-------------
* CERN ROOT  (tested at now with version ROOT 5.34/09
* git 
* cmake
* Tested on ifarm.jlab.org on June 5 2014. Needs CLAS software (and OLD CERN LIBS) 
* Environment variables on ifarm.jlab.org

setenv CERNLIB /apps/cernlib/x86_64_rhel4/2005/lib

setenv CVSROOT /group/clas/clas_cvs

setup cernlib tcl

source /group/clas/builds/environment.csh

Install
-------
This version is running correctly on ifarm.jlab.org. The version of cmake there is too old, so, just use /work/halla/parity/disk2/zana/Cmake/cmake-2.8.8/bin/cmake rather than the default cmake.
* cd EdGen (go to the EdGen directory)
* mkdir build
* cd build
* if you want to install BOS output support you will need to have correctly setup in your environment variables CLAS6LIB CLAS6INC and CERNLIB (the environment show before has been tested)
* if you don't want BOS support, just unset those environment variables (CLAS6LIB CLAS6INC and CERNLIB). If you will request a BOS output the code will not create any bos file (strangely enough)
* cmake ../ ( at jlab /work/halla/parity/disk2/zana/Cmake/cmake-2.8.8/bin/cmake ../ )
* make 

Running
-------
* cd EdGen/output (go to the output directory)
* A table with particle properties is in this directory (from PDG) eg_pdg_table.txt :Modify it if you need it. An example of how to add particles to this table are written at the end of the file. At now the information used are just the mass, the charge of the particle and the lifetime (if the particle has vertex in your reaction diagram for this MonteCarlo. 
* A template is in input.dat (input_test.dat it is just for my developments at now), and modify the file to fit your reaction
* In input.dat one can specify the input spetrum (for example for CLAS photon beam). The input spectrum file is written in a txt files format with raw that represent for each bin E_min, E_max, Counts (Does not need to be normalize, the code is going to normalize it if need it). The bin size do not need to be the same for each bin.  An example is written in the output directory as energy.txt 
* ./EdGen -h will give you the options
* ./EdGen -i input.dat -o output.root        (Change input.dat or output.root (need to be a *.root) to your desired input and output file format 
  
Input file
----------
* nevt:    20000;                  NUMBER OF EVENTS TO GENERATE
* nprint:  1000;                   NUMBER OF EVENTS TO HAVE A PRINTOUT (FOR DEBUGGING)
* model:   2;      		 MODELS AVAILABLE (SEE BOTTOM FOR DIFFERENT OPTIONS)
* M_mode:  1;          MASS MODELS AVAILABLE (SEE BOTTOM FOR DIFFERENT OPTIONS) 
* ifile:	 energy.txt; 		 INPUT FILE SPECTRUM FOR BEAM (NEEDED FOR OPTION MODEL = 2) 
* beam:    22;			 BEAM PARTICLE ID
* en:	 11.0    GeV;		 BEAM ENERGY (NEEDED FOR OPTION MODEL = 1)
* tg_Z:    1;	 		 TARGET Z (NOT WORKING YET)
* tg_N:    1;			 TARGET N (NOT WORKING YET)
* tg_mass  1.876  GeV;           TARGET MASS
* length:	 40	cm;		 LENGTH TARGET
* ras_x:	 0.2	cm;		 BEAM PROFILE (GAUSSIAN SIGMA IN THE X DIRECTION)
* ras_y:	 0.2	cm;		 BEAM PROFILE (GAUSSIAN SIGMA IN THE Y DIRECTION)
* offset:  0.12,0.14,1.1 cm;	 OFFSET TARGET (X,Y,Z)
* npart:   5;	       		 NUMBER OF PARTICLE INVOLVED IN THE INTERACTION (EXCLUDING BEAM AND TARGET)
* pid:     11,2212,113,211,-211;	 PARTICLE ID OF THE PARTICLE SPECIFIED WITH npart
* theta_min:   2.5,4.0,5.0,4.0,5.0 deg;		 THETA CUT FOR SINGLE PARTICLE (FROM 'pid:' flag) (AT NOW IS AN HARD CUT ON THE SIMULATED DATA)
* theta_max:   180.0,180.0,180.0,180.0,180.0 deg;		 THETA CUT FOR SINGLE PARTICLE (FROM 'pid:' flag) (AT NOW IS AN HARD CUT ON THE SIMULATED DATA)
* nvertex: 2;			 NUMBER OF VERTEXES IN THE INTERACTION
* vertex:  0,3;			 ORIGIN OF THE VERTEX (0 STANDS FOR BEAM+TARGET),PARTICLE GOING OUT OF VERTEX (READ IN SEQUENCE FROM pid
* v_type:  1,1.0;		 TYPE OF VERTEX (AT NOW JUST OPTION 1, IN THE FUTURE DIFFERENT CROSS SECTION CAN BE APPLIED AT EACH VERTEX)
* vertex:  3,2;			 ORIGIN OF VERTEX (3 IN THIS CASE IS PARTICLE WITH pid=113), PARTICLE GOING OUT OF VERTEX
* v_type:  1,1.0;		 SAME AS BEFORE
* output:  2;			 OUTPUT FORMAT (SEE BELOW FOR OPTIONS)
* END

Models
-------
* 1 Phase Space Single Energy (for example e-)
* 2 Phase Space Energy Spectrum (for example gamma)
* 3 Cross Section (sorry, not yet)
* 4 Amplitudes (sorry, not yet) 
* 5 Data Points (sorry, not yet)

Mass Models
-------
(Mass model: if width of particle>1MeV, one can generate the mass according to different distributions)
* 1 Breit-Wigner (Mass and Width are automatically read from output/eg_pdg_table.txt)
* 2 Flat (Generated flat in mass in the allowed range)
* 3 Just the mass at the center of the distribution
* 4 Gaussian and more to come


output
-------
* 1  ROOT only
* 2  ROOT + LUND
* 3  ROOT + BOS

Examples
-------
* 3 particles in a vertex, Dalitz plots are generate using the weight (for particles that decay with more than 2 particles in a vertex). The weight is an array of all the particles and describe the weight at creation (the weight is given to the decayed particles) <br />
** Create generated output file: ./EdGen -i input_test2.dat <br />
** Analyze the output (with TProof) of the generated file (files analysis.C , analysis.h, run_analysis.C): root run_analysis.C <br />
* 2 particles per vertex, but 3 vertex <br />
** Create generated output file: ./EdGen -i input_test.dat <br />
** Analyze the output (with TProof) of the generated file (files newAnalysis.C , newAnalysis.h, run_newAnalysis.C): root run_newAnalysis.C <br />
* See other examples of input file with input.dat (default), etc.
