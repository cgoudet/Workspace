#ifndef WORKSPACE_H
#define WORKSPACE_H

#include "Workspace/Category.h"
#include "RooCategory.h"

class Workspace {

 public :

  Workspace();
  ~Workspace();
  Workspace( std::string name );

  void Configure( std::string configFileName );
  void CreateWS();
  void SetDebug( int debug ) { m_debug=debug; }
 private :

  RooWorkspace *m_workspace;
  std::vector<Category*> m_categories;
  std::string m_configFileName;
  std::vector<std::string> m_categoriesNames;
  std::vector<std::string> m_processes;
  //  std::map<std::string, int> m_sDef;
  std::string m_systFileName;
  RooCategory *m_category;
  RooDataSet* addGhosts(RooDataSet* orig, const RooArgSet *observables );
  std::map<std::string, RooArgSet*> m_mapSet;
  std::string m_name;
  int m_debug;

  RooDataSet* MakeAsimovData();

  /**
     \brief Extract all constraint pdf
   */
  void UnfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter);
};
#endif
