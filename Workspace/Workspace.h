#ifndef WORKSPACE_H
#define WORKSPACE_H

#include "Workspace/Category.h"
#include "RooCategory.h"

class Workspace {

 public :

  Workspace();
  ~Workspace();
  Workspace( string name );

  void Configure( string configFileName );
  void CreateWS();

 private :

  RooWorkspace *m_workspace;
  vector<Category*> m_categories;
  string m_configFileName;
  vector<string> m_categoriesNames;
  vector<string> m_processes;
  map<string, int> m_sDef;
  string m_systFileName;
  RooCategory *m_category;
  RooDataSet* addGhosts(RooDataSet* orig, const RooArgSet *observables );
  map<string, RooArgSet*> m_mapSet;
  string m_name;
  bool m_debug;

  void readConstraintFile();
  //  void makeAsimovData( RooRealVar* mH, ModelConfig* mcInWs, bool doConditional, RooWorkspace* combWs, RooAbsPdf* combPdf, RooDataSet* combData, bool b_only);
};
#endif
