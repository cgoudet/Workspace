#include "Workspace/Workspace.h"
#include <iostream>

int main () {


  Workspace ws;
  //ws.Configure( "/afs/in2p3.fr/home/c/cgoudet/private/Couplings/Workspace/python/StatChallenge011.boost" );
  ws.Configure( "/afs/in2p3.fr/home/c/cgoudet/private/Couplings/Workspace/python/StatChallenge011_pdf.boost" );
  ws.CreateWS();
  return 0;
}
