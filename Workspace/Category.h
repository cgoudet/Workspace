#ifndef CATEGORY_H
#define CATEGORY_H

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

  void SetSDef( map<string, int> *sDef ) { m_sDef = sDef; }
  void SetProcesses( vector<string> *processes );
  void SetCategoriesNames( vector<string> *categoriesNames ) { m_categoriesNames = categoriesNames;}
  RooWorkspace *GetWorkspace() { return m_workspace; }

  string GetName() { return m_name; }
  void SetSystFileName( string systFileName ) { m_systFileName = systFileName; }

  void LoadParameters( string configFileName );
  void CreateWS();

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

  map<string, int> *m_sDef;
  vector<string> *m_processes;

  RooRealVar *GetCurrentSyst( int constraint, string NPName, double upVal, double downVal=0 );
  RooRealVar* defineSystematic_Gauss(TString name, double sigma_value, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, TString process="common", double sigmaRightBifurGauss=0);
  RooRealVar* defineSystematic_LogNorm(TString name, double sigma_value, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, TString process="common", double sigmaRightBifurGauss=0);
  RooRealVar* defineSystematic_asymmetric(TString name, double sigma_value_up, double sigma_value_down, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, TString process="common", double sigmaRightBifurGauss=0);

  /**
     m_signalModel %10
     0 : signal per process
     1 : common signal sum of processes
   */
  unsigned int m_signalModel;
  unsigned int m_signalInput;

  RooDataSet *m_dataset;
  string m_dataFileName;
  string m_dataCut;

  void CreateBackgroundModel();
  void CreateSignalModel();
  void CreateSpurious();
  void ReadNuisanceParameters();
  void ReadNuisanceParametersXML();
  void GetData();
  void DefineSet( string set );

  vector<string> m_coef;
  vector<string> m_form;
  vector<string> m_param;
  string m_systFileName;
  string m_outName;
  RooCategory* m_category;
  RooRealVar *GetNuisanceParameter(TString name, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, double sigmaRightBifurGauss );
  RooWorkspace *m_workspace;
  bool m_debug;

  void SelectInputWorkspace( vector<string> &infos );

  TFile *m_readInputFile;
  RooWorkspace *m_readInputWorkspace;
  void SignalFromParameters();
  void SignalFromPdf();
  map<string, string> m_mapPdfInfo;
  map<string, double> m_changeVar;
  vector<string> *m_categoriesNames;
};
#endif
