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
  void SetDebug( int debug ) { m_debug=debug; }
 private :

  RooWorkspace *m_workspace;
  vector<Category*> m_categories;
  string m_configFileName;
  vector<string> m_categoriesNames;
  vector<string> m_processes;
  //  map<string, int> m_sDef;
  string m_systFileName;
  RooCategory *m_category;
  RooDataSet* addGhosts(RooDataSet* orig, const RooArgSet *observables );
  map<string, RooArgSet*> m_mapSet;
  string m_name;
  int m_debug;

  RooDataSet* MakeAsimovData();

  /**
     \brief Extract all constraint pdf
   */
  void UnfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter);
};
#endif
