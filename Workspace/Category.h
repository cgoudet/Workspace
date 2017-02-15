#ifndef CATEGORY_H
#define CATEGORY_H
#include "PlotFunctions/Arbre.h"
#include <string>
using std::string;
#include <map>
using std::map;
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include <vector>
using std::vector;
#include "RooWorkspace.h"
#include "RooDataSet.h"


enum ConstraintModel { NO_CONSTRAINT, GAUSS_CONSTRAINT, LOGNORM_CONSTRAINT, ASYM_CONSTRAINT };
enum BackgroundModel { BKG_MODEL_EXPO, BKG_MODEL_EXPO_POL2, BKG_MODEL_BERN2, BKG_MODEL_BERN3, BKG_MODEL_BERN4 };

class Category {

 public :
  Category();
  Category( string name );
  ~Category();

  void SetProcesses( vector<string> *processes );
  void SetCategoriesNames( vector<string> *categoriesNames ) { m_categoriesNames = categoriesNames;}
  RooWorkspace *GetWorkspace() { return m_workspace; }

  string GetName() { return m_name; }
  void SetSystFileName( string systFileName ) { m_systFileName = systFileName; }

  void LoadParameters( string configFileName );
  void CreateWS();
  void SetDebug( int debug ) { m_debug=debug; }

 private : 
  /**
     Name of the Category, related to the label of the roocategory.
   */
  string m_name;

  string m_correlatedVar;
  /**
     Map containing the the roorealvar necessary for signal parametrization
   */
  map<string, RooRealVar*> m_mapVar;
  map<string, RooFormulaVar*> m_mapFormula;
  map<string, RooAbsPdf*> m_mapPdf;
  map<string, RooArgSet*> m_mapSet;

  map<string, int> m_sDef;
  vector<string> *m_processes;

  RooRealVar *GetCurrentSyst( int constraint, string NPName, double upVal, double downVal=0 );
  RooRealVar* defineSystematic_Gauss(TString name, double sigma_value, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, TString process="common", double sigmaRightBifurGauss=0);
  RooRealVar* defineSystematic_LogNorm(TString name, double sigma_value, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, TString process="common", double sigmaRightBifurGauss=0);
  RooRealVar* defineSystematic_asymmetric(TString name, double sigma_value_up, double sigma_value_down, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, TString process="common", double sigmaRightBifurGauss=0);

  /**
     m_signalModel
     0 : signal per process
     1 : common signal sum of processes
   */

  RooDataSet *m_dataset;
  string m_dataCut;

  void CreateBackgroundModel();
  void CreateSpurious();
  void ReadNuisanceParameters();
  void ReadNuisanceParametersXML();
  void GetData();
  void DefineSet( string set );
  void GetPdfFromWS( const ChrisLib::Arbre &arbre, std::map<std::string, std::stringstream> &editStr );

  string m_systFileName;
  string m_outName;
  RooCategory* m_category;
  RooRealVar *GetNuisanceParameter(TString name, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, double sigmaRightBifurGauss );
  RooWorkspace *m_workspace;
  int m_debug;

  void SelectInputWorkspace( string fileName );

  TFile *m_readInputFile;
  RooWorkspace *m_readInputWorkspace;
  void SignalFromPdf();
  map<string, string> m_mapPdfInfo;
  map<string, double> m_changeVar;
  vector<string> *m_categoriesNames;
  void ReadConstraintFile();

  ChrisLib::Arbre m_catProperties;
};
#endif
