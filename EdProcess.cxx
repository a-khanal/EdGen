#include "EdProcess.h"

EdProcess::EdProcess(const char *file, char *file2){

    finp   = new EdInput(file);
    fout  = new EdOutput(finp,file2);
    fmodel = new EdModel(finp);
    fphy  = new EdPhysics(fmodel);
   
 
    return;
}

void EdProcess::Run(){
    int nevt = finp->GetNevt();

    printf("nevt = %d\n", nevt);

    int evt;
    int nprnt = finp->GetNprnt();
    TString fmtst[4];    
    fmtst[0] = "ROOT";
    fmtst[1] = "ROOT and LUND";
    fmtst[2] = "ROOT and BOS";
    fmtst[3] = "ROOT and A2";
   
    int fmt = finp->Get_fmt();
 
    for( evt = 0; evt < nevt; evt++ ){
      if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
      fphy->MakeEvent(fout, fmodel);
    }
    if (fmt == 2 ) { // also LUND format 
      printf("Creating LUND output \n") ;
      fout->MakeFileLUND();
    }
    if (fmt == 3 ) { // also LUND format 
      printf("Creating BOS output \n") ;
      fout->MakeFileBOS();
    }
    if (fmt == 4 ) { // also LUND format 
      printf("Creating A2 output \n") ;
      fout->MakeFileA2();
    }


    fout->Close();


    

    return;
}

EdProcess::~EdProcess(){
    return;
}
