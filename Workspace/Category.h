#ifndef CATEGORY_H
#define CATEGORY_H
#include "PlotFunctions/Arbre.h"
#include <string>
#include <map>

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include <vector>
#include "RooWorkspace.h"
#include "RooDataSet.h"


enum ConstraintModel { NO_CONSTRAINT, GAUSS_CONSTRAINT, LOGNORM_CONSTRAINT, ASYM_CONSTRAINT };
enum BackgroundModel { BKG_MODEL_EXPO, BKG_MODEL_EXPO_POL2, BKG_MODEL_BERN2, BKG_MODEL_BERN3, BKG_MODEL_BERN4 };

class Category {

 public :
  Category();
  Category( std::string name );
  ~Category();

  void SetProcesses( std::vector<std::string> *processes );
  void SetCategoriesNames( std::vector<std::string> *categoriesNames ) { m_categoriesNames = categoriesNames;}
  RooWorkspace *GetWorkspace() { return m_workspace; }

  std::string GetName() { return m_name; }
  void SetSystFileName( std::string systFileName ) { m_systFileName = systFileName; }

  void LoadParameters( std::string configFileName );
  void CreateWS();
  void SetDebug( int debug ) { m_debug=debug; }

 private : 
  /**
     Name of the Category, related to the label of the roocategory.
   */
  std::string m_name;

  std::string m_correlatedVar;
  /**
     Map containing the the roorealvar necessary for signal parametrization
   */
  std::map<std::string, RooRealVar*> m_mapVar;
  std::map<std::string, RooFormulaVar*> m_mapFormula;
  std::map<std::string, RooAbsPdf*> m_mapPdf;
  std::map<std::string, RooArgSet*> m_mapSet;

  std::map<std::string, int> m_sDef;
  std::vector<std::string> *m_processes;

  RooRealVar *GetCurrentSyst( int constraint, std::string NPName, double upVal, double downVal=0 );
  RooRealVar* defineSystematic_Gauss(TString name, double sigma_value, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, std::string &channel_correlated_np, RooArgSet  *allConstraints, TString process="common", double sigmaRightBifurGauss=0);
  RooRealVar* defineSystematic_LogNorm(TString name, double sigma_value, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, std::string &channel_correlated_np, RooArgSet  *allConstraints, TString process="common", double sigmaRightBifurGauss=0);
  RooRealVar* defineSystematic_asymmetric(TString name, double sigma_value_up, double sigma_value_down, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, std::string &channel_correlated_np, RooArgSet  *allConstraints, TString process="common", double sigmaRightBifurGauss=0);

  /**
     m_signalModel
     0 : signal per process
     1 : common signal sum of processes
   */

  RooDataSet *m_dataset;
  std::string m_dataCut;

  void CreateBackgroundModel();
  void CreateSpurious();
  void ReadNuisanceParameters();
  void ReadNuisanceParametersXML();
  void GetData();
  void DefineSet( std::string set );
  void GetPdfFromWS( const ChrisLib::Arbre &arbre, std::map<std::string, std::stringstream> &editStr );
  void GetYieldFromWS( const ChrisLib::Arbre &arbre );
  void ChangeVars( const ChrisLib::Arbre &arbre );
  std::string m_systFileName;
  std::string m_outName;
  RooCategory* m_category;
  RooRealVar *GetNuisanceParameter(TString name, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, std::string &channel_correlated_np, RooArgSet  *allConstraints, double sigmaRightBifurGauss );
  RooWorkspace *m_workspace;
  int m_debug;

  void SelectInputWorkspace( const std::string &fileName );

  TFile *m_readInputFile;
  RooWorkspace *m_readInputWorkspace;
  void SignalFromPdf();
  std::map<std::string, std::string> m_mapPdfInfo;
  std::map<std::string, double> m_changeVar;
  std::vector<std::string> *m_categoriesNames;
  void ReadConstraintFile();

  ChrisLib::Arbre m_catProperties;
};
#endif
