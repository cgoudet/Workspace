#include "Workspace/Category.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "RooExponential.h"
#include "RooGaussian.h"
#include <iostream>
using std::cout;
using std::endl;
#include "PlotFunctions/SideFunctions.h"
#include "PlotFunctions/DrawPlot.h"
#include "RooProduct.h"
#include "RooBernstein.h"
#include "RooAddition.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooCBShape.h"
#include "RooProdPdf.h"
#include "TFile.h"
#include "TROOT.h"
#include "TObjString.h"
#include "TTree.h"
#include "RooBifurGauss.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "TIterator.h"
#include "TXMLEngine.h"
#include "TSystem.h"
#include "TDOMParser.h"
#include "TXMLDocument.h"
#include "TXMLNode.h"
#include "TXMLAttr.h"
#include "TIterator.h"
using std::ifstream;
using std::stringstream;
using namespace RooStats;
using namespace RooFit;


Category::Category() : m_name( "inclusive" ), m_signalModel(0), m_signalInput(0), m_debug(1), m_catProperties()
{

  m_dataset = 0;
  m_processes = 0;

  m_mapVar["lumi"] = new RooRealVar( "lumi_dum", "lumi_dum", 3.21296+10.0638 );
  m_mapVar["lumi"]->setConstant(1);
  m_mapVar["invMass"] = new RooRealVar ("invariant_mass","invariant_mass",126.5, 110.,160.); 
  m_mapVar["mHcomb"] = new RooRealVar("mHcomb","mHcomb",125.09, 110, 160); // reference is mH = 125 GeV
  m_mapFormula["mHRen"] = new RooFormulaVar("mHRen","mHRen","(@0-100)/100.", RooArgList(*m_mapVar["mHcomb"])); 
  m_mapVar["mu"] = new RooRealVar( "mu", "mu", 1 );
  m_mapVar["mu_BR_yy"] = new RooRealVar( "mu_BR_yy", "mu_BR_yy", 1 );
  m_correlatedVar = "mHRen,mu,mu_BR_yy,one,zero";
  //  m_correlatedVar += ",XS13_ggH_yy,XS13_VBF_yy,XS13_VH_yy,XS13_WH_yy,XS13_ZH_yy,XS13_ttH_yy,XS13_bbH_yy,XS13_tHjb_yy,XS13_tWH_yy";
  m_coef = { "a", "b", "c", "d" };
  m_form = { "CB", "GA", "Var" };
  m_param = { "mean", "sigma", "alpha", "yield" };

  vector<string> setsToDefine = { "pdfProc", "pdfToAdd", "yieldsToAdd", "yieldsSpurious", "observables", "parametersOfInterest", "modelParameters" };
  for ( auto vSet : setsToDefine )   m_mapSet[vSet] = new RooArgSet(vSet.c_str());

  m_mapSet["observables"]->add( *m_mapVar["invMass"] );
  m_mapSet["parametersOfInterest"]->add( *m_mapVar["mu"] );
  m_mapSet["parametersOfInterest"]->add( *m_mapVar["mu_BR_yy"] );
  m_mapSet["parametersOfInterest"]->add( *m_mapVar["mHcomb"] );

  m_workspace = 0;
  m_readInputFile=0;
  m_readInputWorkspace=0;

}

Category::Category( string name ) : Category()
{
  m_name = name;
  //  m_mapSet["parametersOfInterest"]->Print();
  string dumName = "ws_" + m_name;
  m_workspace = new RooWorkspace( dumName.c_str(), dumName.c_str() );


}

Category::~Category() {

  for ( auto vVar = m_mapVar.begin(); vVar != m_mapVar.end(); vVar++ )
    if ( vVar->second ) delete vVar->second;
  for ( auto vVar = m_mapFormula.begin(); vVar != m_mapFormula.end(); vVar++ )
    if ( vVar->second ) delete vVar->second;
  for ( auto vVar = m_mapPdf.begin(); vVar != m_mapPdf.end(); vVar++ )
    if ( vVar->second ) delete vVar->second;
  for ( auto vVar = m_mapSet.begin(); vVar != m_mapSet.end(); vVar++ )
    if ( vVar->second ) delete vVar->second;

  if ( m_workspace ) delete m_workspace;
}

//=========================================
void Category::LoadParameters( string configFileName ) {
  cout << "LoadParameters" << endl;
  //Reas xml configuration file
  cout << configFileName << endl;
  Arbre wsProperties = Arbre::ParseXML( configFileName );
  //  wsProperties.Dump();
  //  exit(0);
  vector<string> vectNodeNames{ "category", "CreateWorkspace" };
  vector< map< string, string > > vectOptions;
  for ( unsigned int i=0; i<vectNodeNames.size(); i++ ) vectOptions.push_back(map< string, string >());
  vectOptions.front()["Name"] = m_name;

  vector<Arbre> vectNodes;
  Arbre::GetArbresPath( wsProperties, vectNodes, vectNodeNames,  vectOptions );
  if ( vectNodes.size() == 1 ) m_catProperties = vectNodes.front();
  else  { cout << "No or too many nodes for the given category : " << m_name << " " << vectNodes.size() << endl; exit(0); }
  //  m_catProperties.Dump();

  vectNodes.clear();
  vectNodeNames = { "correlatedVar", "category" };
  Arbre::GetArbresPath( m_catProperties, vectNodes, vectNodeNames );
  for ( auto vArbre : vectNodes ) m_correlatedVar += "," + vArbre.GetAttribute( "text" );

}
// //=========================================
// void Category::LoadParameters( string configFileName ) {
//   boost::property_tree::ptree pt;
//   boost::property_tree::ini_parser::read_ini(configFileName, pt);
  
//   m_systFileName=pt.get<string>( m_name + ".systFileName" );
//   if ( m_systFileName.find( ".xml" ) == string::npos )  readConstraintFile();

//   vector<string> inputParamInfo( 2, "" );
//   string pdfInfoName;
//   vector<string> processes = *m_processes;
//   processes.push_back( "all" );
//   string name;
//   for ( auto vProc : processes ) {    
//     m_mapPdfInfo["invMass"] = pt.get<string>( m_name + ".invMass", "" );
//     m_mapPdfInfo["mHcomb"] = pt.get<string>( m_name + ".mHcomb", "" );
//     //Loading pdf and formulas directly
//     name = "signal_" + vProc;
//     m_mapPdfInfo[name] = pt.get<string>( m_name + "." + name, "" );
//     if ( m_mapPdfInfo[name] != "" ) m_signalInput=2;
    
//     name = "yield_" + vProc;
//     m_mapPdfInfo[name] = pt.get<string>( m_name + "." + name, "" );
//     if ( m_mapPdfInfo[name] != "" ) m_signalInput=2;
    
//   }//end process

//   m_signalModel = pt.get<unsigned int>( m_name + ".signalModel", 0 );

//   m_dataFileName = pt.get<string>( m_name + ".dataFileName" );
//   cout << "m_dataFileName : " << m_dataFileName << endl;
//   m_dataCut = pt.get<string>( m_name + ".dataCut", "" );
//   m_mapPdfInfo["dataWeight"] = pt.get<string>( m_name + ".dataWeight", "" );
//   cout << "end LoadingParameters" << endl;
//   exit(0);
//   if ( !m_signalInput ) {
//     cout << "no signal input defined" << endl;
//     exit(0);
//   }

//   //Get change of variables
//   vector<string> changeVarName;
//   vector<double> changeVarVal;
//   string dumString = pt.get<string>( m_name + ".changeVarName","" );
//   ParseVector( dumString, changeVarName );
//   dumString=pt.get<string>( m_name + ".changeVarVal","" );
//   ParseVector( dumString, changeVarVal );

//   if ( changeVarVal.size() != changeVarName.size() ) { cout << "chnageVar sizes do not match : Name=" << changeVarName.size() << " Val=" << changeVarVal.size() << endl; exit(0); }
//   for ( unsigned int iName=0; iName< changeVarVal.size(); iName++ )
//     m_changeVar[changeVarName[iName]]= changeVarVal[iName];

// }


//=========================================
void Category::ReadNuisanceParameters() {
  cout << "ReadNuisanceParameter" << endl;
  m_mapSet["systematicValues"] = new RooArgSet();
  for( auto iter : m_sDef ) {
    string name = iter.first;
    m_mapVar[name] = new RooRealVar(name.c_str(),name.c_str(),0); // initialize to 0                  
    //    if ( m_debug )  cout << iter.first << " " << iter.second << " " << m_mapVar[name]->getVal() << endl;                                         
    m_mapSet["systematicValues"]->add(*m_mapVar[name]);
  }
  m_mapSet["systematicValues"]->readFromFile(m_systFileName.c_str(),0,"Common_2015");
  //  m_mapSet["systematicValues"]->readFromFile(m_systFileName,0, TString::Format( "Common_%s", Year.Data() )); // read the common block                     
  m_mapSet["systematicValues"]->readFromFile(m_systFileName.c_str(),0,m_name.c_str()); // read values corresponding to channelname section only.                     

  vector<string> processes( *m_processes );
  processes.push_back( "common" );
  for ( auto vProc = processes.begin(); vProc != processes.end(); vProc++ ) {
    m_mapSet["systematic_yield_"+*vProc] = new RooArgSet();
    m_mapSet["systematic_mass_"+*vProc] = new RooArgSet();
    m_mapSet["systematic_sigma_"+*vProc] = new RooArgSet();
  }
  
  for( auto iter : m_sDef ) {

    TString fullName = iter.first;
    if (fullName == "bkg_model") continue;

    double current_value =  ((RooRealVar*)m_mapSet["systematicValues"]->find(fullName))->getVal();
    // do not add systematics at 0 (surcharge the workspace without valid reason)
    if (current_value==0) continue; 

    double current_err_lo =  ((RooRealVar*)m_mapSet["systematicValues"]->find(iter.first.c_str()))->getMin();
    double current_err_hi =  ((RooRealVar*)m_mapSet["systematicValues"]->find(iter.first.c_str()))->getMax();


    if ( fullName == "dumVar" )  ((RooRealVar*)m_mapSet["systematicValues"]->find(fullName))->Print();
    cout << "systFullName : " << fullName << " " << current_value << endl;   

    //Impose all np to be correlated between categories (same name)    
    bool containsCategory= TString( fullName ).Contains( m_name );
    cout << "containsCategory : " << containsCategory << endl;
    fullName.ReplaceAll( m_name, "" );
    for ( auto vCatName : *m_categoriesNames ) fullName.ReplaceAll( vCatName, "" );


    bool containsPROCESS = false;
    cout << "containsProcess : " << containsPROCESS << endl;
    TString process = "common";
    vector<string> processes = *m_processes;
    processes.push_back( "WH" );
    processes.push_back( "ZH" );
    for ( auto vProc : processes ) {
      if ( !fullName.Contains( vProc ) ) continue;
      containsPROCESS = true;
      process = vProc;
      break;
    }
    fullName.ReplaceAll( process, "" );
    if ( process == "WH" || process =="ZH" ) process = "VH";
    //this line imposes all processes to be correlated
    //Its a shortcut for a non needed possilitu of definesystematic
    string processForName = string(process);
    cout << "process : " << process << endl;

    //This imposes all NP from different years to be correlated
    bool containsYear = fullName.Contains( "2015" )  || fullName.Contains( "2016" );
    cout << "contains Year : " << containsYear << endl;
    fullName.ReplaceAll( "2015", "" );
    fullName.ReplaceAll( "2016", "" );

    CleanName( fullName );
    //Do correlation model
    TString NPName = fullName;
    cout << "NPName : " << NPName << endl;
    if ( NPName == "ATLAS_MET" ) {
      if ( process == "VH" || process == "WH" || process == "ZH" ) NPName+="_VH";
      else NPName+="_nonVH";
    }
    else if ( NPName.Contains("BIAS") )  NPName+="_"+m_name;
    else if ( NPName == "ATLAS_pdf_gg" || NPName== "ATLAS_pdf_qq" ) NPName = iter.first;
    else if ( NPName.Contains("pdf_acc_gg") && process=="ttH" ) NPName+="_"+process;
    else if ( NPName.Contains("pdf_acc") ) NPName+="_"+process+"_"+m_name;
    else if ( NPName.Contains("QCDscale") ) {
      NPName = iter.first;
      NPName.ReplaceAll( m_name, "" );
      if ( process == "VH" ) NPName.ReplaceAll( "WH", "VH").ReplaceAll( "ZH", "VH" );
    }
    else if ( NPName.Contains("syst_Shape_Bkg") ) {
      NPName=iter.first;
      iter.second = LOGNORM_CONSTRAINT;
    }
    else if ( NPName.Contains("syst_Yield_Bkg") ) {
      NPName=iter.first;
      iter.second = LOGNORM_CONSTRAINT;
    }
    else if ( NPName == "ATLAS_TRIGGER_HLT_g35_loose_g25_loose" ) iter.second = LOGNORM_CONSTRAINT;
    else if ( NPName == "ATLAS_LUMI" ) iter.second = LOGNORM_CONSTRAINT;
    else if ( NPName == "DeltaPhi_jj" ) iter.second = LOGNORM_CONSTRAINT;
    else if ( NPName == "EtaStar" ) iter.second = LOGNORM_CONSTRAINT;

    CleanName( NPName );

    //This int is the functional form of the constraint
    int current_constraint = iter.second;
    RooRealVar *current_syst=0;
    if ( current_constraint == ASYM_CONSTRAINT ) current_syst = GetCurrentSyst( current_constraint, string(NPName), current_err_hi, current_err_lo );
    else current_syst = GetCurrentSyst(  current_constraint, string(NPName), current_value );
    cout << "current constraint : " << current_constraint << endl;

    if (NPName.Contains("spurious") || NPName.Contains("BIAS") ) m_mapVar["spurious"]  = current_syst; // the systematics value is the spurious signal
    else if ( NPName.Contains("MSS" ) ) m_mapSet["systematic_mass_"+string(process)]->add(*current_syst);
    else if ( NPName.Contains("MRES" ) ) m_mapSet["systematic_sigma_"+string(process)]->add(*current_syst);
    else if ( NPName.Contains("lhcMass") ) m_mapSet["systematic_mass_"+string(process)]->add(*current_syst);
    else { // means that type == YIELD
      cout << "mapSet : " << m_mapSet[string("systematic_yield_"+process)] << endl;
      m_mapSet[string("systematic_yield_"+process)]->add(*current_syst);    
      cout << "adding" << endl;
    }

  } // end loop on systematics 

  cout << "end nuisanceParameters" << endl;
}//end ReadNuisanceParameter


RooRealVar* Category::defineSystematic_Gauss(TString name, double sigma_value, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, TString process, double sigmaRightBifurGauss) //, bool useBifurGauss=false) 
{
  
  RooRealVar *nui = GetNuisanceParameter(name, nuisance_parameters, global_parameters, constraints_pdf_list, channel_correlated_np, allConstraints, sigmaRightBifurGauss);
  
  TString suffix;
  if (process=="common") suffix = name;
  else  suffix = name+"_"+process;
  // the name of the systematics needs to take into account its process dependance, so add a suffix to it. E.g. JES_ggH
  
  RooRealVar *value = new RooRealVar("value_"+suffix, "value_"+suffix, sigma_value);  
  RooProduct *prod = new RooProduct("prod_"+suffix, "prod_"+suffix, RooArgSet(*nui,*value));  // use the nuisance parameter here
  double central_val = 1.;
  if (name.Contains("spurious") || name.Contains("BIAS") )   central_val = 0.;
  TString oneName = central_val ? "one" : "zero";
  RooRealVar *one = new RooRealVar(oneName, oneName, central_val);
  // cout << "oneName : " << one->GetName() << endl;
  // cout << channel_correlated_np << endl;
  // exit(0);
  RooAddition *systematic = new RooAddition("systematic_"+suffix, "systematic_"+suffix, RooArgSet(*one, *prod));
  return (RooRealVar*) systematic;
}



  /***************************************************************************************************/
  /***************************************************************************************************/
  /***************************************************************************************************/



  RooRealVar* Category::defineSystematic_LogNorm(TString name, double sigma_value, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, TString process, double sigmaRightBifurGauss) //, bool useBifurGauss=false) 
  {
    RooRealVar *nui = GetNuisanceParameter(name, nuisance_parameters, global_parameters, constraints_pdf_list, channel_correlated_np, allConstraints, sigmaRightBifurGauss);

    TString suffix;
    if (process=="common")
      suffix = name;
    else 
      suffix = name+"_"+process;// the name of the systematics needs to take into account its process dependance, so add a suffix to it. E.g. JES_ggH

    RooRealVar *logNorm = new RooRealVar("logNorm_"+suffix, "logNorm_"+suffix, sqrt( log( 1+pow(sigma_value,2))));
    RooExponential *expTerm = new RooExponential("expTerm_"+suffix, "expTerm_"+suffix, *nui, *logNorm);
    RooRealVar *nom_nui = new RooRealVar("nom_nui_"+suffix, "nom_nui_"+suffix, 1);
    //  RooProduct *systematic = new RooProduct("systematic_"+suffix, "systematic_"+suffix, RooArgSet(*nom_nui, *expTerm));
    RooProduct *systematic = new RooProduct("systematic_"+suffix, "systematic_"+suffix, RooArgSet(*nom_nui, *expTerm));

    return (RooRealVar*) systematic;
  }



  /***************************************************************************************************/
  /***************************************************************************************************/
  /***************************************************************************************************/
  RooRealVar* Category::defineSystematic_asymmetric(TString name, double sigma_value_up, double sigma_value_down, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, TString process, double sigmaRightBifurGauss) //, bool useBifurGauss=false)
  {
    TString suffix;
    if (process=="common") 
      suffix = name;
    else 
      suffix = name+"_"+process;// the name of the systematics needs to take into account its process dependance, so add a suffix to it. E.g. JES_ggH

    vector<double> vals_up, vals_down; 
    double pdf_mean = 1.; 
    vals_up.push_back(1+(sigma_value_up/100.)); 
    vals_down.push_back(1/(1-(sigma_value_down/100.))); 
    RooRealVar* nui = GetNuisanceParameter(name, nuisance_parameters, global_parameters, constraints_pdf_list, channel_correlated_np, allConstraints, sigmaRightBifurGauss);
    RooArgList fix;  
    fix.add(*nui);

    HistFactory::FlexibleInterpVar *special_value = new HistFactory::FlexibleInterpVar("asymParam_"+suffix,"asymParam_"+suffix,fix,pdf_mean,vals_down,vals_up); // value that will be asymmetric
    // special_value->Print();
    // see http://www.usatlas.bnl.gov/~fisyak/star/root/roofit/histfactory/src/FlexibleInterpVar.cxx
    //  special_value->setAllInterpCodes(1);
    //  special_value->setAllInterpCodes(1);
    //       special_value->setAllInterpCodes(4);
    RooRealVar *one = new RooRealVar("one", "one", 1);
    RooProduct *systematic = new RooProduct("systematic_"+suffix, "systematic_"+suffix, RooArgSet(*one, *special_value)); // "1+" contained in vals_up and vals_down !
    //  systematic->Print();

    return (RooRealVar*) systematic;

  }

//#################################################
void Category::CreateBackgroundModel() {
  cout << "define the background model : " << m_name << endl;

  vector<Arbre> vectNodes;
  vector<string> vectPath = { "bkg", "category" };
  Arbre::GetArbresPath( m_catProperties, vectNodes, vectPath );
  if ( vectNodes.size() != 1 ) { cout << "too many or too few bkg descriptions : " << vectNodes.size() << endl;exit(0); }
  int bkg_model = BKG_MODEL_EXPO;
  if ( vectNodes.front().GetAttribute( "form" ) == "expPol2" ) bkg_model = BKG_MODEL_EXPO_POL2;

  // ((RooRealVar*)m_mapSet["systematicValues"]->find("bkg_model"))->getVal();
  cout << "invMass : " << m_mapVar["invMass"] << endl;
  RooFormulaVar *x = new RooFormulaVar("x","x", "(@0-100.)/100.", *m_mapVar["invMass"]);  
  x->Print();
  m_mapVar["invMass"]->Print();
  switch (bkg_model ) {
  case  BKG_MODEL_EXPO : {
    RooRealVar *slope = new RooRealVar("slope", "slope", -2e-2, -1, 0.);
    m_mapSet["nuisanceParameters"]->add(*slope);
    RooExponential *exp_pdf = new RooExponential ("bkgExp", "bkgExp", *m_mapVar["invMass"], *slope);
    m_mapPdf["bkg"] = exp_pdf;
    break;
  }
  case  BKG_MODEL_EXPO_POL2 : {
    RooRealVar *p1 = new RooRealVar("p0", "p1", 0, -10, 10);
    RooRealVar *p2 = new RooRealVar("p1", "p2", 0, -10, 10);
    RooArgSet *coefficients = new RooArgSet(*p1,*p2);
    m_mapSet["modelParameters"]->add(*coefficients);
    m_mapSet["nuisanceParameters"]->add(*coefficients);

    RooGenericPdf *exp_pol2 = new RooGenericPdf("bkgExpPol2", "bkgExpPol2", "exp(@1*@0+@2*@0*@0)", RooArgList( *x, *p1, *p2 ) ); // quite long, WHY???
    m_mapPdf["bkg"] = exp_pol2;
    break;
  }
  case  BKG_MODEL_BERN4 :   {
    RooRealVar *a0 = new RooRealVar("a0", "a0", 1);
    RooRealVar *a1 = new RooRealVar("a1", "a1", 0, 0, 20);
    RooRealVar *a2 = new RooRealVar("a2", "a2", 0, 0, 30);
    RooRealVar *a3 = new RooRealVar("a3", "a3", 0, 0, 30);
    RooRealVar *a4 = new RooRealVar("a4", "a4", 0, 0, 10);
    RooArgSet *coefficients = new RooArgSet(*a0, *a1, *a2, *a3, *a4);
    m_mapSet["modelParameters"]->add(*coefficients);
    m_mapSet["nuisanceParameters"]->add(*coefficients);
    RooBernstein *bern4 = new RooBernstein("bkgBern4", "bkgBern4", *x, *coefficients);
    m_mapPdf["bkg"] = bern4;
    break;
  }
  case  BKG_MODEL_BERN3 :   {
    RooRealVar *a0 = new RooRealVar("a0", "a0", 1, -10, 10);
    RooRealVar *a1 = new RooRealVar("a1", "a1", 1, -10, 10);
    RooRealVar *a2 = new RooRealVar("a2", "a2", 1, -10, 10);
    RooRealVar *a3 = new RooRealVar("a3", "a3", 1);
    RooArgList *coefficients = new RooArgList(*a0, *a1, *a2, *a3);


    m_mapSet["modelParameters"]->add(*coefficients);
    m_mapSet["nuisanceParameters"]->add(*coefficients);

    m_mapPdf["bkg"] = new RooBernstein("bkgBern3", "bkgBern3", *m_mapVar["invMass"], *coefficients);
    //    m_mapPdf["bkg"] = new RooBernstein("bkgBern3", "bkgBern3", *x, *coefficients);
    break;
  }
  }

  RooRealVar *nbkg = new RooRealVar( "nbkg", "nbkg", m_dataset ? m_dataset->sumEntries() : 1 , 0, 1e9 );
  nbkg->setConstant(0);
  m_mapPdf["bkg"]->SetNameTitle( "bkg", "bkg" );
  m_mapSet["yieldsToAdd"]->add( *nbkg );
  m_mapSet["pdfToAdd"]->add( *m_mapPdf["bkg"] );
  m_mapSet["modelParameters"]->add(*nbkg);
  m_mapSet["nuisanceParameters"]->add(*nbkg);
  if ( m_dataset && m_dataset->sumEntries() > 3  ) m_mapPdf["bkg"]->fitTo( *m_dataset, SumW2Error(kFALSE) );

  DrawPlot( m_mapVar["invMass"], {m_dataset, m_mapPdf["bkg"] }, "/sps/atlas/c/cgoudet/Plots/HgamBkg_"+m_name, {"nComparedEvents=55"} );
  //  exit(0);
  m_mapVar["invMass"]->Print();
  cout << "DefineBackground done" << endl;
}


//##################################
void Category::CreateSignalModel() {

  // if ( m_signalInput == 1 ) SignalFromParameters();
  // else if ( m_signalInput == 2 ) 
  SignalFromPdf();
  m_mapSet["yieldsSpurious"]->Print();
  m_mapSet["yieldsToAdd"]->Print();

  switch ( m_signalModel ) {
  // case 1 :
  //   for ( unsigned int iProc = 0; iProc < m_processes->size(); iProc++ ) {
  //     string name = "signalSumPdf_" + (*m_processes)[iProc];
  //     m_mapPdf[name] = new RooAddPdf( name.c_str(), name.c_str(), *m_mapSet["pdfProc"], *m_mapSet["yieldsSpurious"] );
  //     if ( m_mapSet["yieldsToAdd"]->find( string( "yieldFactor_" + (*m_processes)[iProc] ).c_str() ) ) m_mapSet["pdfToAdd"]->add( *m_mapPdf[name] );
  //     }
  //   break;
  default :
    //    cout << "pdfProc : " << m_mapSet["pdfProc"] << endl;
    m_mapSet["pdfToAdd"]->add( *m_mapSet["pdfProc"] );
    break;
  }
  //  m_mapSet["pdfToAdd"]->Print();
}

//===================================
void Category::CreateSpurious() {
  if ( m_debug ) cout << "Category::CreateSpurious()" << endl;
  if ( !m_mapVar["spurious"] ) return;

  if ( m_debug ) {
    m_mapSet["pdfProc"]->Print();
    m_mapSet["yieldsSpurious"]->Print();
  }
  RooAddPdf *spuriousPdf = new RooAddPdf( "spuriousPdf", "spuriousPdf", *m_mapSet["pdfProc"], *m_mapSet["yieldsSpurious"] );
  // cout << spuriousPdf << endl;
  m_mapSet["pdfToAdd"]->add( *spuriousPdf );
  m_mapSet["yieldsToAdd"]->add( *m_mapVar["spurious"] );
  //  cout << "spurious done" << endl;
}

//====================================
void Category::CreateWS() {

  GetData();

  //Define datasets of nuisance parameters for ReadNuisanceParameters
  m_mapSet["nuisanceParameters"] = new RooArgSet();
  m_mapSet["globalObservables"] = new RooArgSet();
  m_mapSet["constraintPdf"] = new RooArgSet();
  m_mapSet["allConstraint"] = new RooArgSet();

  string systFileName = m_catProperties.GetAttribute( "systFileName" );
  if ( systFileName.find( ".xml" ) != string::npos ) ReadNuisanceParametersXML();
  else ReadNuisanceParameters();
  
  CreateSignalModel();
  m_mapVar["invMass"]->setRange( 105, 160);

  //If only the pdf all is available (or yield all)
  if ( ( m_mapSet["yieldsToAdd"]->getSize()==1 && m_mapSet["pdfToAdd"]->getSize()!=1  )
       || ( m_mapSet["yieldsToAdd"]->getSize()!=m_mapSet["pdfToAdd"]->getSize() && m_mapSet["pdfToAdd"]->getSize()>1  )
       ) { cout << "ToAdd's sizes do not match" << endl; exit(0); }
  if ( m_mapSet["pdfToAdd"]->getSize()==1 && m_mapSet["yieldsToAdd"]->getSize()>1 ) {
    RooAddition *sumYieldCateg = new RooAddition( "sumYieldCateg", "sumYieldCateg",  *m_mapSet["yieldsToAdd"] );
    //    delete m_mapSet["yieldsToAdd"];
    m_mapSet["yieldsToAdd"] = new RooArgSet( "yieldsToAdd" );
    m_mapSet["yieldsToAdd"]->add( *sumYieldCateg );
    m_mapSet["yieldsSpurious"] = new RooArgSet( "yieldsSpurious" );
    m_mapSet["yieldsSpurious"]->add( *sumYieldCateg );
  }


  CreateSpurious();
  CreateBackgroundModel();

  RooWorkspace *dumWS = 0;
  if ( m_workspace )  {
    dumWS = m_workspace;
    dumWS->SetName( "dumWS" );
    string dumName = "ws_" + m_name;
    m_workspace = new RooWorkspace( dumName.c_str(), dumName.c_str() );
    m_workspace->importClassCode();
  }

  m_mapSet["pdfToAdd"]->Print();
  m_mapPdf["modelSB"] =  new RooAddPdf( "modelSB", "modelSB", *m_mapSet["pdfToAdd"], *m_mapSet["yieldsToAdd"] );
  m_workspace->import( *m_mapPdf["modelSB"], RecycleConflictNodes(), RenameAllVariablesExcept( m_name.c_str(), m_correlatedVar.c_str() ), RenameAllNodes( m_name.c_str()), Silence() );
  cout << "imported SB" << endl;
  // m_workspace->Print();
  // exit(0);
  RooArgSet prodPdf;
  prodPdf.add( *m_workspace->pdf( string( string(m_mapPdf["modelSB"]->GetName()) + "_" + m_name ).c_str() ) );
  cout << "pdfConstraint : " << m_mapSet["constraintPdf"] << endl;
  TIterator* iter = m_mapSet["constraintPdf"]->createIterator();
  RooAbsPdf* parg;
  cout << "importing constraint" << endl;
  while((parg=(RooAbsPdf*)iter->Next()) ) {
    TString name = parg->GetName();
    //    cout << name << endl;
    m_workspace->import( *parg, Silence() );
    prodPdf.add( *m_workspace->pdf( name ) );
  }

  string modelName = "model_" + m_name;
  m_mapPdf["model"] = new RooProdPdf( modelName.c_str(), modelName.c_str(), prodPdf );
  m_workspace->import( *m_mapPdf["model"], RecycleConflictNodes(), Silence() );
  cout << "importing model" << endl;
  cout << m_workspace->var( "mHcomb" ) << endl;
  m_workspace->var( "mHcomb" )->setConstant(1);
  m_workspace->var( "mHcomb" )->setVal(125.09);
  cout << m_workspace->var( "mu" ) << endl;
  m_workspace->var( "mu" )->setConstant(0);
  m_workspace->var( "mu" )->setVal(1);

  for ( auto vMap : m_changeVar ) {
    RooRealVar *changeVar = (RooRealVar*) m_workspace->var( vMap.first.c_str() );
    if ( changeVar ) changeVar->setVal( vMap.second );
  }

  //m_mapPdf["model"]->fitTo( *m_dataset, Strategy(1), SumW2Error(1) );
  DrawPlot( m_mapVar["invMass"], { m_dataset, m_mapPdf["model"] }, "/sps/atlas/c/cgoudet/Plots/HgamModel_"+m_name );

  m_workspace->import( *m_dataset );
  cout << "seting sets" << endl;
  vector<string> sets = { "nuisanceParameters", "globalObservables", "observables", "parametersOfInterest", "modelParameters" };
  for ( auto set : sets ) DefineSet( set );

  m_workspace->importClassCode();
  cout << "end CreateWS" << endl;

  cout << "fitting cat only" << endl;
  // m_mapPdf["model"]->fitTo( *m_dataset, RooFit::SumW2Error(kFALSE) );
  // cout << "fitted BR : " << m_workspace->var( "ATLAS_BR_gamgam" )->getVal() << " " << m_workspace->var( "ATLAS_BR_gamgam" )->getError() << endl;

}

//=========================================
void Category::GetData() {
  //  m_catProperties.Dump();

  vector<string> vectNodeNames;
  vector<Arbre> vectNodes;
  Arbre::GetArbresPath( m_catProperties, vectNodes, { "dataFile", "data", "category" } );
  cout << "nNodes : " << vectNodes.size() << endl;
  vectNodes.front().Dump();

  for ( auto vDataArbre : vectNodes ) {
    RooDataSet *dumDataset = 0;
    string inFileName = vDataArbre.GetAttribute( "inFileName" );
    if ( TString(inFileName).Contains( ".txt" ) ) {
      RooDataSet *newData = RooDataSet::read(m_dataFileName.c_str(), *m_mapSet["observables"]);
      dumDataset = newData;
      continue;
    }
    
    //Get the root file
    TFile *inFile = new TFile( inFileName.c_str() );
    if ( !inFile ) { cout << inFileName << " does not exist." << endl; exit(0);}
    
    //Get the TTree or the workspace
    string datasetName = vDataArbre.GetAttribute( "datasetName" );
    TTree *inTree = (TTree*) inFile->Get( datasetName.c_str() );
    if ( inTree ) {
      delete inTree; inTree=0;
      delete inFile; inFile=0;
      continue;
    }
    
    cout << datasetName  << " does not refer to TTree. Trying RooDataSet in workspace." << endl;
    RooWorkspace *inWS = (RooWorkspace*) inFile->Get( FindDefaultTree( inFile, "RooWorkspace" ).c_str() );
    if ( !inWS ) { cout << "TTree and Workspace failed." << endl; exit(0); }
    
    dumDataset = (RooDataSet*) inWS->data( datasetName.c_str() );
    if ( !dumDataset ) { cout << "dataset failed." << endl; }
    
    string invMassName = vDataArbre.GetAttribute( "varName" );
    if ( invMassName == "" ) { cout << "invariant mass parameter missiong from config file." << endl; exit(0); }
    m_mapVar["invMass"]->SetName( invMassName.c_str() );
    m_correlatedVar += "," + string( m_mapVar["invMass"]->GetName() );
    m_mapSet["observables"]->add( *m_mapVar["invMass"] );
    
    string weightVarName = vDataArbre.GetAttribute( "weightName" );
    if ( weightVarName != "" && 0 ) {
      m_mapVar["dataWeight"] = inWS->var(  weightVarName.c_str() );
      cout << weightVarName << " " << m_mapVar["dataWeight"] << endl;
      m_mapVar["dataWeight"]->SetName( "dataWeight" );
      m_mapSet["observables"]->add( *m_mapVar["dataWeight"] );
    }

    
    datasetName = "obsData_" + m_name;
    string catVarName = vDataArbre.GetAttribute( "catName" );
    string catIndex = vDataArbre.GetAttribute( "catIndex" );
    if ( catVarName != "" && catIndex!="") {
      RooDataSet *dumDataSet=dumDataset;
      RooCategory *eventCateg = (RooCategory*) inWS->cat( catVarName.c_str() );
      m_mapSet["observables"]->add( *eventCateg );
      dumDataset = new RooDataSet( datasetName.c_str(), datasetName.c_str(),  
				  dumDataSet, 
				  *m_mapSet["observables"],
				  catIndex.c_str() );
      }
    
    if ( !m_dataset ) m_dataset = dumDataset;
    else m_dataset->append( *dumDataset );

    // delete inWS; inWS=0;
    // delete inFile; inFile=0;
  }//end for dataFile

  m_dataset->SetName( ("obsData_"+m_name).c_str());
  //  exit(0);
  cout << "end GetData" << endl;
}

//================
void Category::DefineSet( string set ) {
  RooArgSet dumSet( *m_mapSet[set] );
  delete m_mapSet[set];
  m_mapSet[set] = new RooArgSet( set.c_str() );

  TIterator* iter_observable = dumSet.createIterator();
  RooRealVar* parg_observable ;
  while((parg_observable=(RooRealVar*)iter_observable->Next()) ) {
    TString name_of_observable = parg_observable->GetName()+TString("_")+m_name;
    if( m_workspace->obj(name_of_observable) ) m_mapSet[set] -> add( *(RooRealVar*)m_workspace->obj(name_of_observable) );
    else if ( m_workspace->obj(parg_observable->GetName() ) ) m_mapSet[set]-> add( *(RooRealVar*)m_workspace->obj(parg_observable->GetName()) );
  }
    
  m_workspace->defineSet(set.c_str(), *m_mapSet[set], kTRUE);
}

 //===============================
//==============================
RooRealVar *Category::GetNuisanceParameter(TString name, RooArgSet *nuisance_parameters, RooArgSet *global_parameters, RooArgSet *constraints_pdf_list, string &channel_correlated_np, RooArgSet  *allConstraints, double sigmaRightBifurGauss ) //bool useBifurGauss=false) 2
  {
    //    if ( m_debug ) 
    cout << "Category::GetNuisanceParameter " << name << endl;
    RooRealVar *nui = ((RooRealVar*)nuisance_parameters->find(name));
    if (!nui) { // if the np does not exist yet in this channel
      RooRealVar* glob_nui = new RooRealVar("glob_nui_"+name, "glob_nui_"+name, 0, -5, 5);
      glob_nui->setConstant();
      nui = new RooRealVar(name, name, 0, -5, 5 );
      nuisance_parameters->add(*nui);
      global_parameters->add(*glob_nui);

      // (Re)create the constraint
      RooRealVar *sigma_gauss_constraint = new RooRealVar("sigma_"+name, "sigma_"+name, 1);
      channel_correlated_np += "," + string(nui->GetName()) + "," + string(glob_nui->GetName()) +"," +string(sigma_gauss_constraint->GetName());
      RooAbsPdf *constraint;
      if (sigmaRightBifurGauss>0) {
	RooRealVar *sigma_gauss_constraint_right = new RooRealVar("sigma_right_"+name, "sigma_right_"+name, sigmaRightBifurGauss);
	constraint = new RooBifurGauss("constraint_"+name, "constraint_"+name, *nui, *glob_nui, *sigma_gauss_constraint,  *sigma_gauss_constraint_right);
      }
      else {
	constraint = new RooGaussian("constraint_"+name, "constraint_"+name, *nui, *glob_nui, *sigma_gauss_constraint);
      }
      //      if ( name.Contains( "spurious" )  ) constraint->Print();
      if (! allConstraints->find(constraint->GetName())) { // if constraint not yet considered, should be the case since the NP does not exist yet
	constraints_pdf_list->add(*constraint); // add it to the list of constraints to be applied to this channel
	allConstraints->add(*constraint); // add it to the list of all constraints (for all channels)not to apply it twice
      }
      m_correlatedVar += "," + string(constraint->GetName());
    }

    return nui;
  }


//=========================================
void Category::SetProcesses( vector<string> *processes ) {
  m_processes = processes;
  for ( auto vProc : *m_processes ) {
    string muName = "mu_XS_"+ vProc;
    m_mapVar[muName] = new RooRealVar( muName.c_str(), muName.c_str(), 1 );
    m_mapVar[muName]->setConstant(1);
    m_mapSet["parametersOfInterest"]->add( *m_mapVar[muName] );
    m_correlatedVar += "," + muName;
  }
}

//================================================
void Category::SelectInputWorkspace( string fileName ) {

  if ( m_readInputFile && fileName == m_readInputFile->GetName() ) return;
  if ( m_readInputFile ) delete m_readInputFile; m_readInputFile=0; 
  if ( m_readInputWorkspace ) delete m_readInputWorkspace; m_readInputFile=0;
  m_readInputFile = new TFile( fileName.c_str() );
  if ( !m_readInputFile ) {
    cout << fileName.front() << " does not exists" << endl;
    exit(0);
  }
  m_readInputWorkspace = (RooWorkspace*) m_readInputFile->Get( FindDefaultTree( m_readInputFile, "RooWorkspace" ).c_str() );
}

//========================================
void Category::SignalFromParameters() {
  for ( auto vProcess = m_processes->begin(); vProcess != m_processes->end(); vProcess++ ) {
    //Stores the gaussian and CB pdf of the process in this category   
    RooArgSet signalForms;
    for ( auto vForm : m_form ) {
      map<string, RooProduct *> product;
      for ( auto vPar : m_param ) {
	if ( ( vPar == "yield"  && vForm != "Var" )
	     || ( vPar != "yield"  && vForm == "Var" )
	     || (  vPar == "alpha" && vForm!= "CB" )
	     ) continue;
	//defined formulas for cental values and width of signal parametrization
	RooArgList varSet;
	varSet.add( *m_mapVar["mHcomb"] );
	varSet.add( *m_mapFormula["mHRen"] );

	for ( auto vCoef = m_coef.begin(); vCoef != m_coef.end(); vCoef++ ) {
	  string dumString = string(TString::Format( "%s_%s%s_%s", vCoef->c_str(), vPar.c_str(), vForm.c_str(), vProcess->c_str() ));
	  if ( m_mapVar[dumString] ) varSet.add( *m_mapVar[dumString] );

	}//end vCoef

	if ( varSet.getSize() == 2 ) continue;
	string formulaName = string(TString::Format( "%s%s_%s", vPar.c_str(), vForm.c_str(), vProcess->c_str() ));
	//Define the formula for the parameters
	//Polynome is defined according to varSet size
	string formula = "", base = "";
	for ( int iPlot=2; iPlot < varSet.getSize(); iPlot++ ) {
	  if ( formula != "" ) formula+="+";
	  formula+=TString::Format( "@%d%s%s", iPlot, base=="" ? "" : "*", base.c_str() );
	  base+=( base=="" ? "" : "*" ) + string("@1");
	}
	if ( vPar == "mean" ) formula += "+@0";
	if ( vPar == "yield" ) formula = string(TString::Format( "%2.2f*max(0.,%s)", m_mapVar["lumi"]->getVal(), formula.c_str() ) );
	m_mapFormula[formulaName] = new RooFormulaVar(formulaName.c_str(), formulaName.c_str(), formula.c_str(), varSet );
	//Contains mu/sigma and all systematics related
	RooArgSet factors;
	factors.add( *m_mapFormula[formulaName] );

	if ( vPar == "sigma" ) {
	  factors.add( *m_mapSet["systematic_sigma_" + *vProcess] );
	  factors.add( *m_mapSet["systematic_sigma_common"] );
	}
	else if ( vPar == "mean" ) {
	  factors.add( *m_mapSet["systematic_mass_" + *vProcess] );
	  factors.add( *m_mapSet["systematic_mass_common"] );
	}
	formulaName = "prod_" + formulaName;

	product[vPar] = new RooProduct( formulaName.c_str(), formulaName.c_str(), factors );
      }//end vPar

      if ( ( vForm == "GA" || vForm == "CB" ) && ( !product["mean"] || !product["sigma"] ) ) continue;
      string pdfName = string( TString::Format( "%s_%s", vForm.c_str(), vProcess->c_str() ) );

      if ( vForm == "CB" ) {
	string alphaCBName = string(TString::Format("alphaCB_%s", vProcess->c_str() ) );
	string nCBName = string(TString::Format("nCB_%s", vProcess->c_str()) );
	if ( !m_mapFormula[alphaCBName] || !m_mapVar[nCBName] ) continue;
	m_mapPdf[pdfName] =  new RooCBShape( pdfName.c_str(), pdfName.c_str(), *m_mapVar["invMass"],*product["mean"],*product["sigma"],*m_mapFormula[alphaCBName],*m_mapVar[nCBName] );
      }
      else if ( vForm == "GA" ) m_mapPdf[pdfName] = new RooGaussian( pdfName.c_str(), pdfName.c_str(), *m_mapVar["invMass"], *product["mean"], *product["sigma"] );
      cout << vForm << endl;
      if ( vForm == "GA" || vForm == "CB" )  signalForms.add( *m_mapPdf[pdfName] );
      else {
	if ( signalForms.getSize() != 2 || ! product["yield"] ) continue;
	//Create a rooArglist with all variables to multiply yields with
	RooArgList yieldsFactors;
	yieldsFactors.add( *product["yield"] );
	yieldsFactors.add( *m_mapVar["mu_XS_" + *vProcess ] );//muXS
	yieldsFactors.add( *m_mapVar["mu"] );//globalMu
	yieldsFactors.add( *m_mapVar["mu_BR_yy"] );//muBR
	yieldsFactors.add( *m_mapSet["systematic_yield_" + *vProcess] );
	yieldsFactors.add( *m_mapSet["systematic_yield_common"] );
	string yieldFactorName = string(TString::Format( "yieldFactor_%s", vProcess->c_str()));
	RooProduct *yieldFactorProd = new RooProduct( yieldFactorName.c_str(), yieldFactorName.c_str(), yieldsFactors );
	m_mapSet["yieldsToAdd"]->add( *yieldFactorProd );
	m_mapSet["yieldsSpurious"]->add( *product["yield"] );
      }//end if yield

    }//end vForm


    string signalName = string( TString::Format( "signal_%s", vProcess->c_str() ) );
    RooArgList CBFraction = RooArgList( *m_mapVar[string( TString::Format( "fCB_%s", vProcess->c_str()))] );
    RooAddPdf *signalPdf = new RooAddPdf( signalName.c_str(), signalName.c_str(), signalForms, CBFraction );
    if ( m_mapSet["yieldsToAdd"]->getSize() == m_mapSet["pdfProc"]->getSize()+1 ) m_mapSet["pdfProc"]->add( *signalPdf );
  }//end vProcess
}


//=========================================
void Category::SignalFromPdf() { 
  if ( m_debug ) cout << "Category::SignalFromPdf()" <<endl;

  vector<string> processes  = *m_processes;
  processes.insert( processes.begin(), "all" );

  //####### Import all necessary variables from different workspace into a tmporary workspace
  //m_workspace = new RooWorkspace( "m_workspace", "m_workspace" );
  
  //Get alll the Arbres nodes wich refer to  pdf and yield
  map<string, stringstream> editStr; 
  string editedPdfName = "signal";
  

  vector<Arbre> vectNodes;
  vector<string> vectPath = { "", "category" };
  vector<string> tmpNodes = { "pdf", "yield" };
  for ( auto vNodeName : tmpNodes ) {
    vectPath.front() = vNodeName;
    Arbre::GetArbresPath( m_catProperties, vectNodes, vectPath );
  }
  
  for ( auto vArbre : vectNodes ) {
    SelectInputWorkspace( vArbre.GetAttribute( "inFileName" ) );
    string proc = vArbre.GetAttribute( "process" );
    RooAbsReal *absReal = (RooAbsReal*) m_readInputWorkspace->obj( vArbre.GetAttribute( "inVarName" ).c_str() );
    //    cout << vArbre.GetAttribute( "inName" ) << " " << absReal << endl;
    if ( absReal ) m_workspace->import( *absReal, RecycleConflictNodes(), RenameAllVariablesExcept( proc.c_str(), m_correlatedVar.c_str() ), RenameAllNodes( proc.c_str() ), Silence() );

    RooAbsPdf *pdf = (RooAbsPdf*) m_readInputWorkspace->pdf( vArbre.GetAttribute( "inVarName" ).c_str() );
    if ( pdf ) {
      m_workspace->import( *pdf, RecycleConflictNodes(), RenameAllVariablesExcept( proc.c_str(), m_correlatedVar.c_str() ), RenameAllNodes( proc.c_str() ), Silence() );
      editStr[proc] << "EDIT::" << editedPdfName + "_" + proc << "(" << pdf->GetName() << "_" << proc;
    }

  }

  //First loop over ched variables to change value and names
  vectPath.front() = "changeVar";
  vectNodes.clear();
  Arbre::GetArbresPath( m_catProperties, vectNodes, vectPath );
  for ( auto vArbre : vectNodes ) {
    string varName = vArbre.GetAttribute( "inName" );
    RooRealVar *var = m_workspace->var( varName.c_str() );
    vector<string> dumVect  ={ "" };
    if ( !var ) dumVect = processes;
    
    for ( auto vProc : dumVect ) {
      string tmpName = varName + "_" + vProc;
      if ( vProc != "" ) var = m_workspace->var( tmpName.c_str() );
      if ( !var ) continue;
      tmpName = vArbre.GetAttribute( "outName" );
      if ( vProc!="" ) tmpName += "_" + vProc;
      var->SetName( tmpName.c_str() );
      if ( vArbre.GetAttribute( "scale" ) != "" ) ((RooRealVar*) var)->setVal( var->getVal() * stod(vArbre.GetAttribute( "scale" )) );
      else if ( vArbre.GetAttribute( "outVal" ) != "" ) ((RooRealVar*) var)->setVal( stod(vArbre.GetAttribute( "outVal" )) );
    }

  }

  //SEcond loop to change name and parametrization of functions
    map<string,RooArgSet> mapSet;
    mapSet["mean"].add( *m_mapSet["systematic_mass_common"] );
    mapSet["sigma"].add( *m_mapSet["systematic_sigma_common"] );
    for ( auto vProc  : processes ) {
      if ( vProc != "tWH" && vProc != "bbH" && vProc != "tHjb" ) {
	mapSet["yield"].add( *m_mapVar["mu"] );//globalMu
	mapSet["yield"].add( *m_mapVar["mu_BR_yy"] );//muBR
	mapSet["yield"].add( *m_mapSet["systematic_yield_common"] );
      }
      if ( vProc != "all"  ) {
	mapSet["mean"].add(*m_mapSet["systematic_mass_" + vProc] );
	mapSet["sigma"].add( *m_mapSet["systematic_sigma_" + vProc] );
	mapSet["yield"].add( *m_mapSet["systematic_yield_" + vProc] );
	mapSet["yield"].add( *m_mapVar["mu_XS_" + vProc ] );//muXS
      }
    }



  for ( auto vArbre : vectNodes ) {
    string varName = vArbre.GetAttribute( "inName" );
    RooAbsReal *funct = m_workspace->function( varName.c_str() );
    vector<string> dumVect  ={ "" };
    if ( !funct ) dumVect = processes;
    
    for ( auto vProc : dumVect ) {
      string tmpName = varName + "_" + vProc;
      if ( vProc != "" ) funct = m_workspace->function( tmpName.c_str() );
      if ( !funct ) continue;
      string outName = vArbre.GetAttribute( "outName" );
      if ( vProc!="" ) outName += "_" + vProc;

      if ( vArbre.GetAttribute("systNP") == "" ) funct->SetName( outName.c_str() );
      else {
	if ( vArbre.GetAttribute("systNP") == "mean" ) {
	  string depName = "dependentMean";
	  if ( vProc != "" ) depName += "_" + vProc;
	  RooRealVar *mass  = m_workspace->var( "mHcomb" );
	  if ( !mass ) mass = new RooRealVar( "mHcomb", "mHcomb", 125.09 );
	  RooFormulaVar *depMean = new RooFormulaVar( depName.c_str(), depName.c_str(), "@0+(@1-125)", RooArgSet( *funct, *mass ) );
	  funct = depMean;
	}
	RooProduct *prod = new RooProduct( outName.c_str(), outName.c_str(), RooArgSet( *funct, mapSet[vArbre.GetAttribute( "systNP" )] ) );
	m_workspace->import( *prod, Silence(), RecycleConflictNodes() );

	string tmpEdit = "," + tmpName  + "=" + outName;
	if ( vProc == "" ) for ( map<string,stringstream>::iterator it = editStr.begin(); it != editStr.end(); ++it ) it->second << tmpEdit;
	else if ( editStr.find(vProc) != editStr.end() ) editStr[vProc] << tmpEdit;
      }
    } 
  }

  //Close the editStr 
  for ( map<string,stringstream>::iterator it = editStr.begin(); it != editStr.end(); ++it ) {
    it->second << ")";
    cout << it->second.str() << endl;
    m_workspace->factory(it->second.str().c_str());      
  }

  //Get from the workspace the pdf and yields and add them 
  for ( auto vProc : processes ) {
    string tmpName = "signal_" + vProc;
    RooAbsPdf *pdf = m_workspace->pdf( tmpName.c_str() );
    if ( pdf ) {
      //      m_workspace->import( *pdf, RecycleConflictNodes(), Silence() );
      m_mapSet["pdfProc"]->add( *pdf );
    }

    tmpName = "yieldSignal_"+vProc;
    RooAbsReal *yield = m_workspace->function( tmpName.c_str() );
    if ( yield ) {
      //      m_workspace->import( *yield, RecycleConflictNodes(), Silence() );
      //  yield = (RooProduct*) m_workspace->function( tmpName.c_str() );
      m_mapSet["yieldsToAdd"]->add( *yield );
      m_mapSet["yieldsSpurious"]->add( *yield );
    }
  }
  m_workspace->Print();
  //  delete tmpWS; tmpWS=0;



  // vector<string> inputParamInfo;
  // vector<string> varToEdit = { "meanCB", "sigmaCB" };
  // //  m_correlatedVar += "," + m_mapPdfInfo["mHcomb"] + "," + m_mapPdfInfo["invMass"];

  // for ( auto vProc : processes ) {

  //   RooWorkspace *dumWS = new RooWorkspace( "dumWS", "dumWS" );    
  //   string name = "signal_" + vProc;



  //   //Recovering pdf info in this process
  //   vector<Arbre> vectNodes;
  //   vector<string> vectPath = { "pdf", "category" };
  //   vector<map<string,string>> vectOpt;
  //   for ( unsigned int i = 0 ; i<vectPath.size(); ++i ) vectOpt.push_back( map<string,string>() );
  //   vectOpt.front()["process"] = vProc;
  //   Arbre::GetArbresPath( m_catProperties, vectNodes, vectPath , vectOpt );

  //   if ( vectNodes.size() == 1 ) {
  //     Arbre pdfInfo = vectNodes[0];
  //     //Get the workspace for the current process

      
  //     //Retrieve the raw signal pdf
  //     RooAbsPdf *tmpPdf = m_readInputWorkspace->pdf( pdfInfo.GetAttribute( "inVarName" ).c_str() );
  //     if ( !tmpPdf ) {
  // 	cout << "pdf not found : " << pdfInfo.GetAttribute( "inVarName" ) << endl;
  // 	continue;
  //     }
  //     dumWS->import( *tmpPdf, Silence() );

  //     //Give the m-yy variables the proper properties
  //     //if ( !dumWS->var( m_mapPdfInfo["invMass"].c_str() ) ) { cout << m_mapPdfInfo["invMass"].c_str() << " does not exist in dumWS" << endl; exit(0);}
  //     // dumWS->var( m_mapPdfInfo["invMass"].c_str() )->setRange(105,160);
  //     // dumWS->var( m_mapPdfInfo["invMass"].c_str() )->SetName( m_mapVar["invMass"]->GetName() );
      
  //     //start the editing line for the new pdf
  //     stringstream editStr; 
  //     string editedPdfName = "signal";
  //     editStr << "EDIT::" << editedPdfName << "(" << tmpPdf->GetName();

      
  //     //Retrieve the mass from the input workspace
  //     //RooRealVar *mH = m_readInputWorkspace->var(m_mapPdfInfo["mHcomb"].c_str() );
  //     // if ( !mH ) { cout << "mH not found in input workspace with name : " << m_mapPdfInfo["mHcomb"] << endl; exit(0); }

  //     vectPath = { "changeVar", "category" };
  //     vectNodes.clear();
  //     Arbre::GetArbresPath( m_catProperties, vectNodes, vectPath );

  //     //Iterate through node to import all variables which will be chnaged into the temporary workspace
  //     for ( auto vArbre : vectNodes ) {
  // 	RooAbsReal *obj =(RooAbsReal*) m_readInputWorkspace->obj( vArbre.GetAttribute( "inName" ).c_str() );
  // 	if ( obj ) dumWS->import( *obj, Silence(), RecycleConflictNodes() );
  //     }

  //     //Iterate a second time to change the names of roorealvar
  //     for ( auto vArbre : vectNodes ) {
  //     	map<string,string> mapAttr = vArbre.GetAttributes();
  // 	RooRealVar *var = dumWS->var( mapAttr["inName"].c_str() );
  // 	if ( !var ) continue;

  // 	if ( mapAttr.find("scale") != mapAttr.end() ) var->setVal( var->getVal() * std::stod(mapAttr["scale"]) );
  // 	else if ( mapAttr.find("outVal") != mapAttr.end() ) var->setVal( std::stod(mapAttr["outVal"] ) );
  // 	if ( mapAttr["outName"] != "" ) var->SetName( mapAttr["outName"].c_str() );
  //     }

  //     //Run again other the nodes to perform the changes of functions
  //     for ( auto vArbre : vectNodes ) {

  //     	map<string,string> mapAttr = vArbre.GetAttributes();

  // 	RooAbsReal *funct = dumWS->function( mapAttr["inName"].c_str() );
  // 	if ( !funct ) continue;

  // 	if ( mapAttr["outName"].find("meanCB") != string::npos ) { 
  // 	  RooRealVar *mass = dumWS->var( "mHcomb" );
  // 	  string dumString = "dependentMean_" + vProc;
  // 	  RooFormulaVar *dependentMass = new RooFormulaVar( dumString.c_str(), "@0+(@1-125)", RooArgSet( *funct, *mass ) );
  // 	  dependentMass->Print();
  // 	  RooProduct *prod = new RooProduct( mapAttr["outName"].c_str(), mapAttr["outName"].c_str(), RooArgSet( *dependentMass, mapSet["mean"] ) );
  // 	  dumWS->import( *prod, Silence(), RecycleConflictNodes() );
  // 	}
  // 	else {
  // 	  if ( mapAttr["systNP"] == "" ) funct->SetName( mapAttr["outName"].c_str() );
  // 	  else {
  // 	    RooProduct *prod = new RooProduct( mapAttr["outName"].c_str(), mapAttr["outName"].c_str(), RooArgSet( *funct, mapSet[mapAttr["systNP"]] ) );
  // 	    editStr << "," << funct->GetName() << "=" << prod->GetName();
  // 	    dumWS->import( *prod, Silence(),RecycleConflictNodes() );
  // 	  }
  // 	}
  //     }

  //     //Edit the pdf with new names and included systematics
  //     if ( m_debug ) cout << "Editing pdf" << endl;
  //     editStr << ")";
  //     cout << editStr << endl;
  //     dumWS->factory(editStr.str().c_str());      
  //     tmpPdf = dumWS->pdf( editedPdfName.c_str() );
      
  //     //import the final pdf for the process in the final workspace
  //     editedPdfName += "_" + vProc;

  //     m_workspace->import( *tmpPdf, RecycleConflictNodes(), RenameAllVariablesExcept( vProc.c_str(), m_correlatedVar.c_str() ), RenameAllNodes( vProc.c_str() ), Silence() );

  //     m_mapSet["pdfProc"]->add( *m_workspace->pdf( editedPdfName.c_str() ) );
  //     if ( m_debug ) cout << editedPdfName << " imported" << endl;
  //   }//end if pdf for process exist
  //   //----------------

  //   //Start getting the yields
  //   vectPath = { "yield", "category" };
  //   vectNodes.clear();
  //   Arbre::GetArbresPath( m_catProperties, vectNodes, vectPath , vectOpt );

  //   name = "yield_" + vProc;
  //   if ( m_debug ) cout << "starting yields : " << vProc << " " << name << " " << vectNodes.size() << endl;
  //   if ( vectNodes.size() > 1 ) cout << "Too many yield node with the process : " << vProc << endl;
  //   else if ( vectNodes.size() == 1 ) {
  //     Arbre yieldArbre = vectNodes.front();
  //     map<string,string> mapAttrYield = yieldArbre.GetAttributes();


  //     RooAbsReal *yieldLumi = 0;
      
  //     if ( mapAttrYield["inFileName"].find( ".txt" ) != string::npos ) {
  // 	//read the constant yield from txt file
  // 	fstream stream;
  // 	stream.open( inputParamInfo.front().c_str(), fstream::in );
  // 	int category;
  // 	double yield, mass;
  // 	while ( stream >> mass >> category >> yield ) {
  // 	  if ( category != mapAttrYield["inCatIndex"] ) continue;
  // 	  yieldLumi = new RooRealVar( "yieldPerFb", "yieldPerFb", yield );
  // 	  break;
  // 	}

  //     }
  //     else yieldLumi = m_readInputWorkspace->function( mapAttrYield["inVarName"].c_str() );
      
  //     if ( !yieldLumi ) {
  // 	cout << " No input yield for " << m_name << " " << vProc << endl;
  // 	continue;
  //     }
      
  //     name = "yield";
  //     yieldLumi->SetName( name.c_str() );
  //     mapSet["yield"].add( *yieldLumi );

  //     name = "yieldFactor";
  //     RooProduct *yield = new RooProduct( name.c_str(), name.c_str(), mapSet["yield"] );
  //     name = name + "_" + vProc;
  //     cout << "importYields" << endl;
  //     cout << m_correlatedVar << endl;
  //     m_workspace->Print();
  //     m_workspace->import( *yield, RenameAllVariablesExcept( vProc.c_str(), m_correlatedVar.c_str() ), RenameAllNodes( vProc.c_str() ), RecycleConflictNodes(), Silence() );

  //     if ( m_name == "ttHlep" && vProc == "ttH" ) {
  // 	cout << "process : " << vProc << endl;
  // 	cout << m_correlatedVar << endl;
  // 	//	exit(0);
  //     }

  //     yield = (RooProduct*) m_workspace->function( name.c_str() );
  //     m_mapSet["yieldsToAdd"]->add( *yield );
  //     m_mapSet["yieldsSpurious"]->add( *yield );

  //   }
    // delete dumWS; dumWS=0;
    // }//end vProc


}
 
//=============================================================
void Category::ReadNuisanceParametersXML() {
  if ( m_debug ) cout << "ReadNuisanceParameterXML" << endl;

  //Define variable which will contain the spurious signal
  m_mapVar["spurious"] = 0;

  TDOMParser xmlparser;
  //Check if the xml file is ok                                                                                                                                                                      
  xmlparser.ParseFile( m_systFileName.c_str() );
  TXMLDocument* xmldoc = xmlparser.GetXMLDocument();
  TXMLNode *rootNode  = xmldoc->GetRootNode();
  TXMLNode *systNode = rootNode->GetChildren();

  //Loop over systematics
  while ( systNode!=0 ) {
    map<string, string> mapAttr;
    //Affects all the attributes to a map
    TList *systAttr = systNode->GetAttributes();
    if(systAttr!=0) {
      TIterator *it = systAttr->MakeIterator();
      for ( auto attr = (TXMLAttr*) it->Next(); attr!=0; attr=(TXMLAttr*)it->Next() ) {
    	mapAttr[attr->GetName()] = attr->GetValue();
      }
    }
    TXMLNode *systEffectNode = systNode->GetChildren();
    systNode = systNode->GetNextNode();

    if ( mapAttr["Name"] == "" || !systAttr ) continue;
    if ( mapAttr["Name" ] == "bkg_model" ) continue;
    if ( m_debug ) cout << "systName : " << mapAttr["Name"] << endl;
    //Deal with empty attributes
    if ( mapAttr["correlation"]=="" ) mapAttr["correlation"] = "All";

    int constraint = GAUSS_CONSTRAINT;
    if (  mapAttr["constraint"]=="LogNorm" ) constraint = LOGNORM_CONSTRAINT;
    else if ( mapAttr["constraint"]=="Asym" ) constraint = ASYM_CONSTRAINT;

    //Loop over all effect in all categories and processes of the current systematic

    while ( systEffectNode != 0 ) {
      map<string,string> mapSystEffect;
      TList *systEffectAttr = systEffectNode->GetAttributes();
      if(systEffectAttr!=0) {
	TIterator *it = systEffectAttr->MakeIterator();
	for ( auto attr = (TXMLAttr*) it->Next(); attr!=0; attr=(TXMLAttr*)it->Next() ) {
	  cout << attr->GetName() << " " << attr->GetValue() << endl;
	  mapSystEffect[attr->GetName()] = attr->GetValue();
	}
      }
      systEffectNode = systEffectNode->GetNextNode();	
      if ( !systEffectAttr ) continue;
      if ( mapSystEffect["upVal"] == 0 ) continue;
      //If the systematic effect isnt common or in the category do nothing
      if ( mapSystEffect["category"] != m_name && mapSystEffect["category"]!="common" ) continue;
      if ( m_debug ) cout << "category/Process : " << mapSystEffect["category"] << " " << mapSystEffect["process"] << endl;

	  
      //Dealing with empty parameters
      if ( mapAttr["process"]=="" ) mapAttr["process"] = "All";
      if ( mapAttr["varName"]=="" ) mapAttr["varName"] = "yield";
      
      //Create the name of the nuisance parameter
      string NPName = mapAttr["Name"];
      if ( (mapAttr["correlation"] == "None" || mapAttr["correlation"] == "Category") && mapSystEffect["process"]!="All" ) NPName += "_" + mapSystEffect["process"];
      if ( mapAttr["correlation"] == "None" || mapAttr["correlation"] == "Process" )  NPName += "_" + m_name;
      cout << NPName << endl;
	  
      RooRealVar *currentSyst = GetCurrentSyst( constraint, NPName, std::stod( mapSystEffect["upVal"] ) , std::stod( mapSystEffect["downVal"] ) );
      
      if (mapAttr["Name"].find("spurious") != string::npos || mapAttr["Name"].find("BIAS") != string::npos ) m_mapVar["spurious"]  = currentSyst; // the systematics value is the spurious signal
      else {
	string setName = "systematic_" + mapSystEffect["varName"] + "_" + mapSystEffect["process"];
	cout << "setName : " << setName << endl;
	if ( !m_mapSet[setName] )  m_mapSet[setName] = new RooArgSet();
	m_mapSet[setName]->add( *currentSyst );
	m_mapSet[setName]->Print();
      }
      
    }//end systEffectNode
  }//end systNode

  exit(0);

}//end ReadNuisanceParameter


//===========================================================
RooRealVar *Category::GetCurrentSyst( int constraint, string NPName, double upVal, double downVal ) {

  string processForName = "common";

  RooRealVar *currentSyst = 0;
  switch (constraint ) {
  case GAUSS_CONSTRAINT :
    currentSyst  = defineSystematic_Gauss(NPName, upVal,
					   m_mapSet["nuisanceParameters"],
					   m_mapSet["globalObservables"],
					   m_mapSet["constraintPdf"],
					   m_correlatedVar,
					   m_mapSet["allConstraint"],
					   processForName);
    break;
  case LOGNORM_CONSTRAINT :
    currentSyst  = defineSystematic_LogNorm(NPName, upVal,
  					     m_mapSet["nuisanceParameters"],
  					     m_mapSet["globalObservables"],
  					     m_mapSet["constraintPdf"],
  					     m_correlatedVar,
  					     m_mapSet["allConstraint"],
					     processForName );
    break;
  case ASYM_CONSTRAINT : 
    currentSyst  = defineSystematic_asymmetric(NPName, upVal, downVal,
						m_mapSet["nuisanceParameters"],
						m_mapSet["globalObservables"],
						m_mapSet["constraintPdf"],
						m_correlatedVar,
						m_mapSet["allConstraint"],
						processForName );
    break;
  default : //if no constraint
    break;
  }//end switch on constraint

  return currentSyst;
}

//=======================================================
void Category::readConstraintFile()
{
  ifstream current_file(m_systFileName.c_str());
  do 
    {
      string tmpString; 
      getline(current_file, tmpString);
      if(tmpString.size() == 0) continue;
      if ( tmpString.find( "#" ) != string::npos ) continue;

      int defConstraint = -1;
      TString tmp(tmpString);
      TObjArray tmpAr = *(tmp.Tokenize(" "));

      if ( tmpAr.GetEntries() == 2 ) {
	TString tmpStrDefConst= ((TObjString*) tmpAr.At(1))->GetString();
	if(tmpStrDefConst == "NO_CONSTRAINT") defConstraint = NO_CONSTRAINT;
	else if(tmpStrDefConst == "GAUSS_CONSTRAINT") defConstraint = GAUSS_CONSTRAINT;
	else if(tmpStrDefConst == "LOGNORM_CONSTRAINT") defConstraint = LOGNORM_CONSTRAINT;
	else if(tmpStrDefConst == "ASYM_CONSTRAINT") defConstraint = ASYM_CONSTRAINT;
	else  defConstraint = GAUSS_CONSTRAINT;
      }
      else {
	if ( TString(tmpString).Contains("-100 L(" ) ) defConstraint = ASYM_CONSTRAINT;
	else defConstraint = GAUSS_CONSTRAINT;
	//	  cout << tmpString << " " << defConstraint << endl;
      }//end else 

      m_sDef[string(((TObjString*) tmpAr.First())->GetString())] = defConstraint;
    } while(!current_file.eof());
}
