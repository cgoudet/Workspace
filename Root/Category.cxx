#include "Workspace/Category.h"
#include "PlotFunctions/SideFunctions.h"
#include "PlotFunctions/SideFunctionsTpp.h"
#include "PlotFunctions/DrawPlot.h"
#include "PlotFunctions/Foncteurs.h"
using namespace ChrisLib;

#include "RooExponential.h"
#include "RooGaussian.h"
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
using namespace RooStats;
using namespace RooFit;

//#include <boost/algorithm/string.hpp>

#include <iostream>
using std::cout;
using std::endl;
#include <exception>
using std::runtime_error;
using std::ifstream;
using std::stringstream;
using std::to_string;

#include <algorithm>
using std::copy;
#include <iterator>
using std::ostream_iterator;
using std::distance;

Category::Category() : m_name( "inclusive" ), m_debug(1), m_catProperties()
{

  m_dataset = 0;
  m_processes = 0;

  m_mapVar["lumi"] = new RooRealVar( "lumi_dum", "lumi_dum", 3.21296+10.0638 );
  m_mapVar["lumi"]->setConstant(1);
  m_mapVar["invMass"] = new RooRealVar ("m_yy","m_yy",126.5, 110.,160.); 
  m_mapVar["invMass"]->setRange( 105, 160);
  m_mapVar["mHcomb"] = new RooRealVar("mHcomb","mHcomb",125.09, 110, 160); // reference is mH = 125 GeV
  m_mapFormula["mHRen"] = new RooFormulaVar("mHRen","mHRen","(@0-100)/100.", RooArgList(*m_mapVar["mHcomb"])); 
  m_mapVar["mu"] = new RooRealVar( "mu", "mu", 1, -10, 10 );
  m_mapVar["mu_BR_yy"] = new RooRealVar( "mu_BR_yy", "mu_BR_yy", 1 );
  m_correlatedVar = "mHRen,mu,mu_BR_yy,one,zero";

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

//#################################################
Category::Category( string name ) : Category()
{
  m_name = name;
  string dumName = "ws_" + m_name;
  m_workspace = new RooWorkspace( dumName.c_str(), dumName.c_str() );
}

//#################################################
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
  if ( m_debug )  cout << "Category::LoadParameters" << endl;
  //Read xml configuration file
  Arbre wsProperties = Arbre::ParseXML( configFileName );

  vector<string> vectNodeNames{ "category", "CreateWorkspace" };
  vector< map< string, string > > vectOptions;
  for ( unsigned int i=0; i<vectNodeNames.size(); i++ ) vectOptions.push_back(map< string, string >());
  vectOptions.front()["Name"] = m_name;
  vector<Arbre> vectNodes;
  Arbre::GetArbresPath( wsProperties, vectNodes, vectNodeNames,  vectOptions );
  if ( vectNodes.size() == 1 ) m_catProperties = vectNodes.front();
  else throw runtime_error( "Category::LoadParameters : No or too many nodes for the given category : " + m_name + " " + to_string(vectNodes.size()) );

  vectNodes.clear();
  vectNodeNames = { "correlatedVar", "category" };
  Arbre::GetArbresPath( m_catProperties, vectNodes, vectNodeNames );
  for ( auto vArbre : vectNodes ) m_correlatedVar += "," + vArbre.GetAttribute( "text" );
  if ( m_debug ) cout << "Category::LoadParameters done\n";
}

//=========================================
void Category::ReadNuisanceParameters() {
  if ( m_debug ) cout << "Category::ReadNuisanceParameter\n";

  ReadConstraintFile();

  m_mapSet["systematicValues"] = new RooArgSet();
  for( auto iter : m_sDef ) {
    string name = iter.first;
    m_mapVar[name] = new RooRealVar(name.c_str(),name.c_str(),0); // initialize to 0                  
    m_mapSet["systematicValues"]->add(*m_mapVar[name]);
  }

  string systFileName = m_catProperties.GetAttribute( "systFileName" );;
  if ( m_debug ) cout << "systFileName : " << systFileName << endl;
  m_mapSet["systematicValues"]->readFromFile(systFileName.c_str(),0,"Common_2015");
  m_mapSet["systematicValues"]->readFromFile(systFileName.c_str(),0,m_name.c_str()); // read values corresponding to channelname section only.                     

  
  for( auto iter : m_sDef ) {

    string fullName = iter.first;
    if (fullName == "bkg_model") continue;

    double current_value = static_cast<RooRealVar*>(m_mapSet["systematicValues"]->find(fullName.c_str()))->getVal();
    // do not add systematics at 0 (surcharge the workspace without valid reason)
    if (current_value==0) continue; 

    double current_err_lo =  static_cast<RooRealVar*>(m_mapSet["systematicValues"]->find(iter.first.c_str()))->getMin();
    double current_err_hi =  static_cast<RooRealVar*>(m_mapSet["systematicValues"]->find(iter.first.c_str()))->getMax();

    cout << "systFullName : " << fullName << " " << current_value << " " << current_err_lo << " " << current_err_hi << endl;   

    //Impose all np to be correlated between categories (same name)    
    bool containsCategory= fullName.find(m_name)!=string::npos;
    if ( m_debug ) cout << "containsCategory : " << containsCategory << endl;
    if ( containsCategory ) fullName = ReplaceString( "_"+m_name)( fullName);
    for ( auto vCatName : *m_categoriesNames ) fullName = ReplaceString( "_"+vCatName)(fullName);


    bool containsPROCESS = false;
    TString process = "common";
    vector<string> processes = *m_processes;
    processes.push_back( "WH" );
    processes.push_back( "ZH" );
    for ( auto vProc : processes ) {
      if ( fullName.find(vProc)==string::npos ) continue;
      containsPROCESS = true;
      process = vProc;
      break;
    }
    fullName = ReplaceString("_"+string(process))(fullName);
    if ( process == "WH" || process =="ZH" ) process = "VH";
    //this line imposes all processes to be correlated
    //Its a shortcut for a non needed possilitu of definesystematic
    string processForName = string(process);
    if ( m_debug ) cout << "process : " << process << endl;

    //This imposes all NP from different years to be correlated
    bool containsYear = fullName.find("2015")!=string::npos  || fullName.find("2016" )!=string::npos;
    if ( m_debug ) cout << "contains Year : " << containsYear << endl;
    fullName = ReplaceString("_2015")( ReplaceString("_2016")(fullName));;

    //Do correlation model
    TString NPName = fullName;
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

    //    CleanName( NPName );

    //This int is the functional form of the constraint
    int current_constraint = iter.second;
    RooRealVar *current_syst=0;
    if ( current_constraint == ASYM_CONSTRAINT ) current_syst = GetCurrentSyst( current_constraint, string(NPName), current_err_hi, current_err_lo );
    else current_syst = GetCurrentSyst(  current_constraint, string(NPName), current_value );

    if (NPName.Contains("spurious") || NPName.Contains("BIAS") ) m_mapVar["spurious"]  = current_syst; // the systematics value is the spurious signal
    else if ( NPName.Contains("MSS" ) ) m_mapSet["systematic_mass_"+string(process)]->add(*current_syst);
    else if ( NPName.Contains("MRES" ) ) m_mapSet["systematic_sigma_"+string(process)]->add(*current_syst);
    else if ( NPName.Contains("lhcMass") ) m_mapSet["systematic_mass_"+string(process)]->add(*current_syst);
    else { // means that type == YIELD
      m_mapSet[string("systematic_yield_"+process)]->add(*current_syst);    
    }

  } // end loop on systematics 

  if ( m_debug ) cout << "Category::ReadNuisanceParameter end\n";
}//end ReadNuisanceParameter

//===========================================
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
    if (process=="common") suffix = name;
    else suffix = name+"_"+process;// the name of the systematics needs to take into account its process dependance, so add a suffix to it. E.g. JES_ggH
    cout << "sigmaValues : " << sigma_value_up << " " << sigma_value_down << endl;
    vector<double> vals_up, vals_down; 
    double pdf_mean = 1.; 
    vals_up.push_back(1+(sigma_value_up/100.)); 
    vals_down.push_back(1/(1-(sigma_value_down/100.)));
    // vals_up.push_back(1+sigma_value_up/100.);
    // vals_down.push_back(1-sigma_value_down/100.);
    TString NPName = name;
    NPName.ReplaceAll( "_mean", "" );
    NPName.ReplaceAll( "_sigma", "" );
    RooRealVar* nui = GetNuisanceParameter( NPName, nuisance_parameters, global_parameters, constraints_pdf_list, channel_correlated_np, allConstraints, sigmaRightBifurGauss);
    RooArgList fix;  
    fix.add(*nui);

    HistFactory::FlexibleInterpVar *special_value = new HistFactory::FlexibleInterpVar("asymParam_"+suffix,"asymParam_"+suffix,fix,pdf_mean,vals_down,vals_up); // value that will be asymmetric
    cout << "nominal Value : " << special_value->getVal() << endl;
    nui->setVal(1);
    cout << "upValue : " << special_value->getVal() << endl;
    nui->setVal(-1);
    cout << "downValue : " << special_value->getVal() << endl;
    nui->setVal(0);
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

  RooFormulaVar *x = new RooFormulaVar("x","x", "(@0-100.)/100.", *m_mapVar["invMass"]);  
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
  if ( m_dataset && m_dataset->sumEntries() > 3  ) {
    m_mapPdf["bkg"]->fitTo( *m_dataset, SumW2Error(kFALSE), Verbose(0) );
    nbkg->setRange( nbkg->getVal()/2., nbkg->getVal()*2);
  }
  //DrawPlot( m_mapVar["invMass"], {m_dataset, m_mapPdf["bkg"] }, "/sps/atlas/c/cgoudet/Plots/HgamBkg_"+m_name, {"nComparedEvents=55"} );
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
  m_mapSet["pdfToAdd"]->add( *spuriousPdf );
  m_mapSet["yieldsToAdd"]->add( *m_mapVar["spurious"] );
}

//====================================
void Category::CreateWS() {
  if ( m_debug ) cout << "Category::CreateWS\n";

  GetData();

  //Define datasets of nuisance parameters for ReadNuisanceParameters
  m_mapSet["nuisanceParameters"] = new RooArgSet();
  m_mapSet["globalObservables"] = new RooArgSet();
  m_mapSet["constraintPdf"] = new RooArgSet();
  m_mapSet["allConstraint"] = new RooArgSet();
  vector<string> processes( *m_processes );
  processes.push_back( "common" );
  for ( auto vProc = processes.begin(); vProc != processes.end(); ++vProc ) {
    m_mapSet["systematic_yield_"+*vProc] = new RooArgSet();
    m_mapSet["systematic_mass_"+*vProc] = new RooArgSet();
    m_mapSet["systematic_sigma_"+*vProc] = new RooArgSet();
  }

  string systFileName = m_catProperties.GetAttribute( "systFileName" );
  if ( systFileName.find( ".xml" ) != string::npos ) ReadNuisanceParametersXML();
  else ReadNuisanceParameters();
  
  SignalFromPdf();
  m_mapSet["pdfToAdd"]->add( *m_mapSet["pdfProc"] );
  
  //  m_mapVar["invMass"]->setRange( 105, 160);
  //m_mapVar["invMass"]->setRange( 110, 160);
  
  //If only the pdf all is available (or yield all)
  if ( ( m_mapSet["yieldsToAdd"]->getSize()==1 && m_mapSet["pdfToAdd"]->getSize()!=1  )
       || ( m_mapSet["yieldsToAdd"]->getSize()!=m_mapSet["pdfToAdd"]->getSize() && m_mapSet["pdfToAdd"]->getSize()>1  )
       ) { cout << "ToAdd's sizes do not match" << endl; exit(0); }
  if ( m_mapSet["pdfToAdd"]->getSize()==1 && m_mapSet["yieldsToAdd"]->getSize()>1 ) {
    RooAddition *sumYieldCateg = new RooAddition( "sumYieldCateg", "sumYieldCateg",  *m_mapSet["yieldsToAdd"] );
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

  m_workspace->import( *m_dataset );
  m_dataset = static_cast<RooDataSet*>( m_workspace->data( m_dataset->GetName() ) );

  m_mapPdf["modelSB"] =  new RooAddPdf( "modelSB", "modelSB", *m_mapSet["pdfToAdd"], *m_mapSet["yieldsToAdd"] );
  //  m_workspace->import( *m_mapPdf["modelSB"], RecycleConflictNodes(), RenameAllVariablesExcept( m_name.c_str(), m_correlatedVar.c_str() ), RenameAllNodes( m_name.c_str()), Silence() );
  cout << m_correlatedVar << endl;
  vector<string> parsed;
  ParseVector(m_correlatedVar, parsed, ',' );
  cout << "parsedSize : " << parsed.size() << endl;
  copy( parsed.begin(), parsed.end(), std::ostream_iterator<string>(cout,"\n"));
  cout << endl;
  m_workspace->import( *m_mapPdf["modelSB"], RecycleConflictNodes(), RenameAllVariablesExcept( m_name.c_str(), m_correlatedVar.c_str() ), RenameAllNodes( m_name.c_str()) );
  string wsName = string(m_mapPdf["modelSB"]->GetName()) + "_" + m_name;
  m_mapPdf["modelSB"] = m_workspace->pdf( wsName.c_str() );
  cout << "imported SB" << endl;

  RooArgSet prodPdf;
  prodPdf.add( *m_mapPdf["modelSB"] );
  TIterator* iter = m_mapSet["constraintPdf"]->createIterator();
  RooAbsPdf* parg;
  cout << "importing constraint" << endl;
  while((parg=(RooAbsPdf*)iter->Next()) ) {
    TString name = parg->GetName();
    m_workspace->import( *parg, Silence() );
    prodPdf.add( *m_workspace->pdf( name ) );
  }


  string modelName = "model_" + m_name;
  m_mapPdf["model"] = new RooProdPdf( modelName.c_str(), modelName.c_str(), prodPdf );
  m_workspace->import( *m_mapPdf["model"], RecycleConflictNodes(), Silence() );
  m_mapPdf["model"] = m_workspace->pdf( m_mapPdf["model"]->GetName() );
  cout << "importing model" << endl;
  // cout << m_workspace->var( "mHcomb" ) << endl;
  // m_workspace->var( "mHcomb" )->setConstant(1);
  // m_workspace->var( "mHcomb" )->setVal(125.09);
  // cout << m_workspace->var( "mu" ) << endl;
  // m_workspace->var( "mu" )->setConstant(0);
  // m_workspace->var( "mu" )->setVal(0);

  cout << "seting sets" << endl;
  vector<string> sets = { "nuisanceParameters", "globalObservables", "observables", "parametersOfInterest", "modelParameters" };
  for ( auto set : sets ) DefineSet( set );

  m_workspace->importClassCode();
  if ( m_debug ) cout << "Category::CreateWS done" << endl;

}

//=========================================
void Category::GetData() {
  if ( m_debug ) cout << "Category::GetData\n";

  vector<string> vectNodeNames;
  vector<Arbre> vectNodes;
  Arbre::GetArbresPath( m_catProperties, vectNodes, { "dataFile", "data", "category" } );
  //  if ( m_debug ) vectNodes.front().Dump();

  for ( auto vDataArbre : vectNodes ) {

    string varName = vDataArbre.GetAttribute( "varName" );
    string inFileName = vDataArbre.GetAttribute( "inFileName" );
    string weightName = vDataArbre.IsAttribute( "weightName" ) ? vDataArbre.GetAttribute( "weightName" ) : "";

    RooDataSet *dumDataset = 0;
    if ( inFileName.find(".txt")!=string::npos ) {
      RooDataSet *newData = RooDataSet::read(inFileName.c_str(), *m_mapSet["observables"]);
      dumDataset = newData;
      continue;
    }
    
    //Get the root file
    TFile *inFile = new TFile( inFileName.c_str() );
    if ( !inFile ) throw runtime_error( "Category::GetData : " + inFileName + " does not exist." );

    //Get the TTree or the workspace
    if (  vDataArbre.IsAttribute( "treeName" )  ) {
      TTree *inTree = static_cast<TTree*>(inFile->Get( vDataArbre.GetAttribute( "treeName" ).c_str()) );
      if ( !inTree ) throw runtime_error( "Category::GetData : " + vDataArbre.GetAttribute( "treeName" ) + " not found in " + inFile->GetName() );

      RooRealVar m_yy( varName.c_str(), varName.c_str(), 125 );
      if ( !m_mapSet["observables"]->find( m_yy.GetName() ) ) m_mapSet["observables"]->add( m_yy );
      RooRealVar weight;
      if ( weightName != "" ) {
	weight = RooRealVar( weightName.c_str(), weightName.c_str(), 1 );
	m_mapSet["observables"]->add( weight );
      }
      if ( m_debug ) cout << "WeightName :  " << weight.GetName() << endl;

      string selectionCut = vDataArbre.GetAttribute( "selectionCut" );
      vector<string> selectionVars;
      ParseVector( vDataArbre.GetAttribute( "selectionVars" ), selectionVars, ',');
      for ( auto vVar : selectionVars ) {
	RooRealVar *selectionVar = new RooRealVar( vVar.c_str(), vVar.c_str(), 0 );
	if ( !m_mapSet["observables"]->find( selectionVar->GetName() ) ) m_mapSet["observables"]->add( *selectionVar );
      }
      
      if ( weightName != "" ) dumDataset = new RooDataSet( "dumDataset", "dumDataset", inTree, *m_mapSet["observables"], selectionCut.c_str(), weightName.c_str() );
      else dumDataset = new RooDataSet( "dumDataset", "dumDataset", inTree, *m_mapSet["observables"], selectionCut.c_str() );

      delete inTree; inTree=0;
    }
    else if ( vDataArbre.IsAttribute( "datasetName" ) ) {
      string datasetName = vDataArbre.GetAttribute( "datasetName" );
      RooWorkspace *inWS = static_cast<RooWorkspace*>(inFile->Get( FindDefaultTree( inFile, "RooWorkspace" ).c_str() ));
      if ( !inWS ) throw runtime_error(  "Category::GetData : TTree and Workspace failed." );
      
      dumDataset = static_cast<RooDataSet*>( inWS->data( datasetName.c_str() ));
      if ( !dumDataset ) throw runtime_error( "Category::GetData : dataset failed.");
      dumDataset->get()->first()->SetName(m_mapVar["invMass"]->GetName());

      string catVarName = vDataArbre.IsAttribute( "catName" ) ? vDataArbre.GetAttribute( "catName" ) : "";
      string catIndex = vDataArbre.IsAttribute( "catIndex" ) ? vDataArbre.GetAttribute( "catIndex" ) : "";
      if ( catVarName!="" && catIndex!="") {
	RooDataSet *dumDataSet=dumDataset;
	RooCategory *eventCateg = (RooCategory*) inWS->cat( catVarName.c_str() );
	m_mapSet["observables"]->add( *eventCateg );
	dumDataset = new RooDataSet( datasetName.c_str(), datasetName.c_str(), dumDataSet, *m_mapSet["observables"], catIndex.c_str() );
      }
    }
    else throw runtime_error( "Category::GetData() : No datasetName or treeName given : data could not be read." );

    if ( !m_dataset ) {
      m_dataset = dumDataset;
      m_dataset->SetName( ("obsData_" + m_name).c_str() );
    }
    else m_dataset->append( *dumDataset );

  }//end for dataFile
  m_dataset->SetName( ("obsData_"+m_name).c_str());
  m_dataset->Print();

  if ( m_debug ) cout << "Category::GetData done" << endl;
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
    cout << "Category::GetNuisanceParameter " << name << endl;
    RooRealVar *nui = static_cast<RooRealVar*>(nuisance_parameters->find(name));
    if (!nui) { // if the np does not exist yet in this channel
      RooRealVar* glob_nui = new RooRealVar("glob_nui_"+name, "glob_nui_"+name, 0, -5, 5);
      glob_nui->setConstant();
      nui = new RooRealVar(name, name, 0, -5, 5 );
      nuisance_parameters->add(*nui);
      global_parameters->add(*glob_nui);

      // (Re)create the constraint
      RooRealVar *sigma_gauss_constraint = new RooRealVar("sigma_"+name, "sigma_"+name, 1);
      channel_correlated_np += "," + string(nui->GetName()) + "," + string(glob_nui->GetName());
      //+"," +string(sigma_gauss_constraint->GetName());
      RooAbsPdf *constraint;
      if (sigmaRightBifurGauss>0) {
	RooRealVar *sigma_gauss_constraint_right = new RooRealVar("sigma_right_"+name, "sigma_right_"+name, sigmaRightBifurGauss);
	constraint = new RooBifurGauss("constraint_"+name, "constraint_"+name, *nui, *glob_nui, *sigma_gauss_constraint,  *sigma_gauss_constraint_right);
      }
      else {
	constraint = new RooGaussian("constraint_"+name, "constraint_"+name, *nui, *glob_nui, *sigma_gauss_constraint);
      }
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
    m_mapVar[muName] = new RooRealVar( muName.c_str(), muName.c_str(), 1, -10, 10 );
    m_mapVar[muName]->setConstant(1);
    m_mapSet["parametersOfInterest"]->add( *m_mapVar[muName] );
    m_correlatedVar += "," + muName;
  }
}

//================================================
void Category::SelectInputWorkspace( const string &fileName ) {

  if ( m_readInputFile && fileName == m_readInputFile->GetName() ) return;
  if ( m_readInputFile ) delete m_readInputFile; m_readInputFile=0; 
  if ( m_readInputWorkspace ) delete m_readInputWorkspace; m_readInputFile=0;
  m_readInputFile = new TFile( fileName.c_str() );
  if ( !m_readInputFile ) throw runtime_error( "Category::SelectInputWorkspace : " + fileName + " does not exists.");
  m_readInputWorkspace = static_cast<RooWorkspace*>(m_readInputFile->Get( FindDefaultTree( m_readInputFile, "RooWorkspace" ).c_str() ));
}
//=========================================
void Category::GetPdfFromWS( const Arbre &arbre, map<string, stringstream> &editStr ) {
  if ( m_debug ) cout << "Category::GetPdfFromWS\n";
  if ( arbre.GetNodeName() != "pdf" ) throw runtime_error( "Category:GetPdfFromWS : wrong arbre type." );
  cout << "inFileName : " << arbre.GetAttribute( "inFileName" ) << endl;
  SelectInputWorkspace( arbre.GetAttribute( "inFileName" ) );
  string proc = arbre.GetAttribute( "process" );
  RooRealVar *mass = m_readInputWorkspace->var( arbre.GetAttribute( "invMass" ).c_str() );
  cout << "mass : " << mass << endl;
  mass->Print();

  if ( mass ) mass->SetName( m_mapVar["invMass"]->GetName() );
  RooAbsPdf *pdf = static_cast<RooAbsPdf*>(m_readInputWorkspace->pdf( arbre.GetAttribute( "inVarName" ).c_str() ));
  if ( !pdf ) throw runtime_error( "Category::GetPdfFromWS : pdf not found.");
  pdf->Print();
  m_workspace->import( *pdf, RecycleConflictNodes(), RenameAllVariablesExcept( proc.c_str(), m_correlatedVar.c_str() ), RenameAllNodes( proc.c_str() ), Silence() );
  editStr[proc] << "EDIT::signal_" + proc << "(" << pdf->GetName() << "_" << proc;
  if ( m_debug ) cout << "Category::GetPdfFromWS end\n";
}

//=========================================
void Category::GetYieldFromWS( const Arbre &arbre ) {
  if ( m_debug ) cout << "Category::GetYieldFromWS\n";
  string inFileName = arbre.IsAttribute( "inFileName" ) ? arbre.GetAttribute( "inFileName" ) : "";
  string proc = arbre.GetAttribute( "process" );
  if ( inFileName.find(".root")!=string::npos ) { 
    SelectInputWorkspace( arbre.GetAttribute( "inFileName" ) );
    RooAbsReal *absReal = static_cast<RooAbsReal*>(m_readInputWorkspace->obj( arbre.GetAttribute( "inVarName" ).c_str() ));
    if ( absReal ) m_workspace->import( *absReal, RecycleConflictNodes(), RenameAllVariablesExcept( proc.c_str(), m_correlatedVar.c_str() ), RenameAllNodes( proc.c_str() ), Silence() );
  }
  else if ( inFileName.find(".txt")!=string::npos ) { 
    fstream stream( inFileName.c_str() );
    if ( !stream.is_open() ) throw runtime_error( "Category::FetYieldFromWS : Unknow file " + inFileName );
    double mass{-99}, yield{-99};
    int cat{-1};
    int nCat = distance(m_categoriesNames->begin(), find( m_categoriesNames->begin(), m_categoriesNames->end(), m_name ) );
    while ( stream >> mass >> cat >> yield && cat!=nCat ) {}
    if ( cat!=nCat ) throw runtime_error( "Category::GetYiedFromWS : category not found in file " + inFileName );

    RooRealVar *yieldVar = new RooRealVar( "yieldSignal","yieldSignal", yield );
    m_workspace->import( *yieldVar, RecycleConflictNodes(), RenameAllVariablesExcept( proc.c_str(), m_correlatedVar.c_str() ), Silence() );
  }
  if ( m_debug ) cout << "Category::GetYieldFromWS end\n";
}
//=========================================
void Category::SignalFromPdf() { 
  if ( m_debug ) cout << "Category::SignalFromPdf()" <<endl;

  vector<string> processes  = *m_processes;
  processes.insert( processes.begin(), "all" );

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
  
  for ( auto vArbre : vectNodes ) { //list over all input files
    string inFileName = vArbre.IsAttribute( "inFileName" ) ? vArbre.GetAttribute( "inFileName" ) : "";
    if ( inFileName == "" ) continue;
    if ( vArbre.GetNodeName() == "pdf" ) GetPdfFromWS( vArbre, editStr );
    else GetYieldFromWS( vArbre );
  }

  cout << "passed yield" << endl;
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
      if ( vArbre.GetAttribute( "outName" ) != "" ) var->SetName( tmpName.c_str() );
      if ( vArbre.GetAttribute( "scale" ) != "" ) ((RooRealVar*) var)->setVal( var->getVal() * stod(vArbre.GetAttribute( "scale" )) );
      else if ( vArbre.GetAttribute( "outVal" ) != "" ) ((RooRealVar*) var)->setVal( stod(vArbre.GetAttribute( "outVal" )) );
    }
  }

  //SEcond loop to change name and parametrization of functions
  cout << "mapSet" << endl;
  map<string,RooArgSet> mapSet;
  if ( m_mapSet["systematic_mean_all"] ) mapSet["mean"].add( *m_mapSet["systematic_mean_all"] );
  if ( m_mapSet["systematic_sigma_all"] ) mapSet["sigma"].add( *m_mapSet["systematic_sigma_all"] );
  
  for ( auto vProc  : processes ) {
    cout << vProc << endl;
    if ( vProc != "tWH" && vProc != "bbH" && vProc != "tHjb" ) {
      mapSet["yield"].add( *m_mapVar["mu"] );//globalMu
      mapSet["yield"].add( *m_mapVar["mu_BR_yy"] );//muBR
      if ( m_mapSet["systematic_yield_all"] ) mapSet["yield"].add( *m_mapSet["systematic_yield_all"] );
    }
    if ( vProc != "all"  ) {
      if ( m_mapSet["systematic_mean_" + vProc] ) mapSet["mean_"+vProc].add(*m_mapSet["systematic_mean_" + vProc] );
      if ( m_mapSet["systematic_sigma_" + vProc] ) mapSet["sigma_"+vProc].add( *m_mapSet["systematic_sigma_" + vProc] );
      if ( m_mapSet["systematic_yield_" + vProc] ) mapSet["yield_"+vProc].add( *m_mapSet["systematic_yield_" + vProc] );
      mapSet["yield_"+vProc].add( *m_mapVar["mu_XS_" + vProc ] );//muXS
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
      if ( string(funct->ClassName() ) == "RooRealVar" ) continue;
      string outName = vArbre.GetAttribute( "outName" );
      if ( vProc!="" ) outName += "_" + vProc;

      if ( vArbre.GetAttribute("systNP") == "" && vArbre.GetAttribute("outName") != "" ) funct->SetName( outName.c_str() );
      else {
	if ( vArbre.GetAttribute("systNP") == "mean" ) {
	  string depName = "dependentMean";
	  if ( vProc != "" ) depName += "_" + vProc;
	  RooRealVar *mass  = m_workspace->var( "mHcomb" );
	  if ( !mass ) mass = new RooRealVar( "mHcomb", "mHcomb", 125.09 );
	  RooFormulaVar *depMean = new RooFormulaVar( depName.c_str(), depName.c_str(), "@0+(@1-125)", RooArgSet( *funct, *mass ) );
	  funct = depMean;
	}
	RooArgSet setForProd;
	setForProd.add( *funct );
	setForProd.add( mapSet[vArbre.GetAttribute( "systNP" )] );
	setForProd.add( mapSet[vArbre.GetAttribute( "systNP" )+"_"+vProc] );
	RooProduct *prod = new RooProduct( outName.c_str(), outName.c_str(), setForProd );
	m_workspace->import( *prod, Silence(), RecycleConflictNodes() );

	string tmpEdit = "," + tmpName  + "=" + outName;
	if ( vProc == "" ) for ( map<string,stringstream>::iterator it = editStr.begin(); it != editStr.end(); ++it ) it->second << tmpEdit;
	else if ( editStr.find(vProc) != editStr.end() ) editStr[vProc] << tmpEdit;
      }
    } 
  }
  cout << "editing" << endl;

  //Close the editStr 
  for ( map<string,stringstream>::iterator it = editStr.begin(); it != editStr.end(); ++it ) {
    it->second << ")";
    m_workspace->factory(it->second.str().c_str());      
  }

  //Get from the workspace the pdf and yields and add them 
  for ( auto vProc : processes ) {
    string tmpName = "signal_" + vProc;
    RooAbsPdf *pdf = m_workspace->pdf( tmpName.c_str() );
    if ( pdf ) m_mapSet["pdfProc"]->add( *pdf );

    tmpName = "yieldSignal_"+vProc;
    RooAbsReal *yield = m_workspace->function( tmpName.c_str() );
    if ( yield ) {
      m_mapSet["yieldsToAdd"]->add( *yield );
      m_mapSet["yieldsSpurious"]->add( *yield );
    }
  }
}
 
//=============================================================
void Category::ReadNuisanceParametersXML() {
  if ( m_debug ) cout << "Category::ReadNuisanceParameterXML" << endl;

  //Define variable which will contain the spurious signal
  m_mapVar["spurious"] = 0;

  TDOMParser xmlparser;
  //Check if the xml file is ok                                                                                                                                                                      
  xmlparser.ParseFile( m_catProperties.GetAttribute( "systFileName" ).c_str() );
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
    if ( mapAttr["varName"]=="" ) mapAttr["varName"] = "yield";

    //Loop over all effect in all categories and processes of the current systematic
    while ( systEffectNode != 0 ) {
      map<string,string> mapSystEffect;
      TList *systEffectAttr = systEffectNode->GetAttributes();
      if(systEffectAttr!=0) {
	TIterator *it = systEffectAttr->MakeIterator();
	for ( auto attr = (TXMLAttr*) it->Next(); attr!=0; attr=(TXMLAttr*)it->Next() ) {
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
      if ( mapSystEffect["process"]=="" ) mapSystEffect["process"] = "all";


      int constraint = GAUSS_CONSTRAINT;
      if (  mapSystEffect["constraint"]=="LogNorm" ) constraint = LOGNORM_CONSTRAINT;
      else if ( mapSystEffect["constraint"]=="Asym" ) constraint = ASYM_CONSTRAINT;
      
      //Create the name of the nuisance parameter
      string NPName = mapAttr["Name"];
      if ( (mapAttr["correlation"] == "None" || mapAttr["correlation"] == "Process") && mapSystEffect["process"]!="All" ) NPName += "_" + mapSystEffect["process"];
      if ( mapAttr["correlation"] == "None" || mapAttr["correlation"] == "Category" )  NPName += "_" + m_name;

      RooRealVar *currentSyst = GetCurrentSyst( constraint, NPName, std::stod( mapSystEffect["upVal"] ) , mapSystEffect["downVal"]=="" ? 0 : std::stod( mapSystEffect["downVal"] ) );
      cout << "NPValues : " << mapAttr["correlation"] << " " << NPName << " " << currentSyst->GetName() << endl;      
      if (mapAttr["Name"].find("spurious") != string::npos || mapAttr["Name"].find("BIAS") != string::npos ) m_mapVar["spurious"]  = currentSyst; // the systematics value is the spurious signal

      else {
	string setName = "systematic_" + mapAttr["varName"] + "_" + mapSystEffect["process"];
	if ( !m_mapSet[setName] )  m_mapSet[setName] = new RooArgSet();
	m_mapSet[setName]->add( *currentSyst );
	m_mapSet[setName]->Print();
      }
      
    }//end systEffectNode
  }//end systNode
  if ( m_debug ) cout << "Category::ReadNuisanceParameterXML end" << endl;
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
void Category::ReadConstraintFile() {
  if ( m_debug ) cout << "Category::ReadConstraintFile\n";
  ifstream current_file(m_catProperties.GetAttribute( "systFileName" ).c_str());
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
  if ( m_debug ) cout << "Category::ReadConstraintFile end\n";
}
