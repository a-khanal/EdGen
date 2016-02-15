#include "EdProcess.h"

EdProcess::EdProcess(const char *file, char *file2){

    finp   = new EdInput(file);
    fout  = new EdOutput(finp,file2);
    fmodel = new EdModel(finp);
    fphy  = new EdPhysics(fmodel);
    fcross = new EdCrossSection(fmodel);
 
    return;
}

void EdProcess::Run(){
    int nevt = finp->GetNevt();

    printf("nevt = %d\n", nevt);

    int evt;
    int nprnt = finp->GetNprnt();
    TString fmtst[3];    
    fmtst[0] = "ROOT";
    fmtst[1] = "ROOT and LUND";
    fmtst[2] = "ROOT and BOS";
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


    fout->Close();


    

    return;
}

EdProcess::~EdProcess(){
    return;
}
