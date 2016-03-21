#include "Workspace/Category.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "RooExponential.h"
#include "RooGaussian.h"
#include <iostream>
using std::cout;
using std::endl;
#include "PlotFunctions/SideFunctions.h"
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


using namespace RooStats;
using namespace RooFit;

Category::Category() : m_name( "inclusive" ), m_debug(0)
{
  m_sDef = 0;
  m_dataset = 0;
  m_processes = 0;

  m_mapSet["parametersOfInterest"] = new RooArgSet("parametersOfInterest");
  m_mapVar["lumi"] = new RooRealVar( "lumi", "lumi", 10 );
  m_mapVar["invMass"] = new RooRealVar ("invariant_mass","invariant_mass",126.5, 110.,160.); 
  m_mapVar["mHcomb"] = new RooRealVar("mHcomb","mHcomb",125, 110, 160); // reference is mH = 125 GeV
  m_mapFormula["mHRen"] = new RooFormulaVar("mHRen","mHRen","(@0-100)/100.", RooArgList(*m_mapVar["mHcomb"])); 
  m_mapVar["mu"] = new RooRealVar( "mu", "mu", 1 );
  m_mapVar["mu_BR_yy"] = new RooRealVar( "mu_BR_yy", "mu_BR_yy", 1 );
  m_correlatedVar == "mHcomb,mHRen,mu,mu_BR_yy";
  for ( auto vProc = m_processes->begin(); vProc != m_processes->end(); vProc++ ) {
    string muName = "mu_XS_"+ *vProc;
    m_mapVar[muName] = new RooRealVar( muName.c_str(), muName.c_str(), 1 );
    m_mapSet["parametersOfInterest"]->add( *m_mapVar[muName] );
    m_correlatedVar += "," + muName;
  }

  m_coef = { "a", "b", "c", "d" };
  m_form = { "CB", "GA", "Var" };
  m_param = { "mean", "sigma", "alpha", "yield" };

  m_mapSet["pdfProc"] = new RooArgSet("pdfProc");
  m_mapSet["pdfToAdd"] = new RooArgSet("pdfToAdd");
  m_mapSet["yieldsToAdd"] = new RooArgSet("yieldsToAdd");
  m_mapSet["yieldsSpurious"] = new RooArgSet("yieldsSpurious");
  m_mapSet["observables"] = new RooArgSet("observables");
  m_mapSet["observables"]->add( *m_mapVar["invMass"] );

  m_mapSet["parametersOfInterest"]->add( *m_mapVar["mu"] );
  m_mapSet["parametersOfInterest"]->add( *m_mapVar["mu_BR_yy"] );
  m_mapSet["parametersOfInterest"]->add( *m_mapVar["mHcomb"] );


}

Category::Category( string name ) : Category()
{
  m_name = name;
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

}


void Category::LoadParameters( string configFileName ) {

  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(configFileName, pt);

  for ( auto vCoef = m_coef.begin(); vCoef != m_coef.end(); vCoef++ ) {
    for ( auto vForm = m_form.begin(); vForm != m_form.end(); vForm++ ) {
      for ( auto vParam = m_param.begin(); vParam != m_param.end(); vParam++ ) {
	string name = string( TString::Format("%s_%s%s", vCoef->c_str(), vParam->c_str(), vForm->c_str() ) );
	m_mapVar[name] = new RooRealVar( name.c_str(), name.c_str(), pt.get<double>( m_name + "." + name ) );
      }
    }
  }

  m_dataFileName = pt.get<double>( m_name + ".dataFileName" );
  for ( auto vVar = m_mapVar.begin(); vVar != m_mapVar.end(); vVar++ )
    if ( vVar->second ) vVar->second->setConstant(1);

}


//=========================================
void Category::ReadNuisanceParameters() {

  m_mapSet["systematicValues"] = new RooArgSet();
  for( auto iter = m_sDef->begin(); iter != m_sDef->end(); iter++) {
    string name = iter->first;
    //    cout << iter->first << " " << iter->second << endl;                                         
    m_mapVar[name] = new RooRealVar(name.c_str(),name.c_str(),0); // initialize to 0                  
    m_mapSet["systematicValues"]->add(*m_mapVar[name]);
  }
  
  m_mapSet["systematicValues"]->readFromFile(m_systFileName.c_str(),0,"Common");
  //  m_mapSet["systematicValues"]->readFromFile(m_systFileName,0, TString::Format( "Common_%s", Year.Data() )); // read the common block                     
  m_mapSet["systematicValues"]->readFromFile(m_systFileName.c_str(),0,m_name.c_str()); // read values corresponding to channelname section only.                     

  vector<string> processes( *m_processes );
  processes.push_back( "common" );
  for ( auto vProc = m_processes->begin(); vProc != m_processes->end(); vProc++ ) {
    m_mapSet["systematicYield_"+*vProc] = new RooArgSet();
    m_mapSet["systematicPeak_"+*vProc] = new RooArgSet();
    m_mapSet["systematicResolution_"+*vProc] = new RooArgSet();
  }


  m_mapVar["spurious"] = 0;
  m_mapSet["nuisanceParameters"] = new RooArgSet();
  m_mapSet["globalObservable"] = new RooArgSet();
  m_mapSet["constraintPdf"] = new RooArgSet();
  m_mapSet["allConstraint"] = new RooArgSet();
  
  for( auto iter = m_sDef->begin(); iter != m_sDef->end(); iter++) {
    TString fullName = iter->first;
    if (fullName == "bkg_model") continue;

    TObjArray* Strings = fullName.Tokenize( "_" );
    bool containsTYPE = false;
    TString type =  ((TObjString*) Strings->First())->GetString();
    cout << "type : " << type << endl;
    if (type == "YIELD" || type ==  "PES" || type == "PER" || type ==  "ESS" || type == "MRES" ) containsTYPE = true;
    else type = "YIELD";

    bool containsPROCESS = false;
    TString process =  ((TObjString*) Strings->Last())->GetString();
    if ( SearchVectorBin( string(process), *m_processes ) != m_processes->size())  containsPROCESS = true;
    else process = "common";

    TString NPname = fullName;
    if (containsPROCESS) {
      // remove also the underscore before the process name
      NPname.Replace(NPname.Index("_"+process),NPname.Length(),""); 
      if (((TObjString*) Strings->At(Strings->GetEntries()-2))->GetString()==process)
 	NPname+="_"+process;
    }

    if (containsTYPE) cout  << endl;// NPname.Replace(0, NPname.Index(type)+type.Length()+1,""); // remove also the underscore after the type name

    double current_value =  ((RooRealVar*)m_mapSet["systematicValues"]->find(fullName))->getVal();

    // do not add systematics at 0 (surcharge the workspace without valid reason)
    if (current_value==0) continue; 

    //This int is the functional form of the constraint
    int current_constraint = iter->second;  
    RooRealVar *current_syst;

    switch ( current_constraint ) {
    case GAUSS_CONSTRAINT :
      current_syst  = defineSystematic_Gauss(NPname, current_value,
 					     m_mapSet["nuisanceParameters"],
 					     m_mapSet["globalObservable"],
 					     m_mapSet["constraintPdf"],
 					     m_correlatedVar,
 					     m_mapSet["allConstraint"],
 					     process);
      break;
    case LOGNORM_CONSTRAINT :
      current_syst  = defineSystematic_LogNorm(NPname, current_value,
 					     m_mapSet["nuisanceParameters"],
 					     m_mapSet["globalObservable"],
 					     m_mapSet["constraintPdf"],
 					     m_correlatedVar,
 					     m_mapSet["allConstraint"],
 					       process );
      break;
    case ASYM_CONSTRAINT : {
      /* if (current_value > -99) */
      /* 	cout << "WARNING: the asymmetric systematics will not used the central value for the parameter !" << endl; */
      double current_err_lo =  ((RooRealVar*)m_mapSet["systematicValues"]->find(fullName))->getMin();
      double current_err_hi =  ((RooRealVar*)m_mapSet["systematicValues"]->find(fullName))->getMax();
      if (m_debug) {
 	cout << "asymmetric error with values : +" << current_err_hi << ", and " << current_err_lo << endl;
      }
      current_syst  = defineSystematic_asymmetric(NPname, current_err_hi, current_err_lo,
						  m_mapSet["nuisanceParameters"],
						  m_mapSet["globalObservable"],
						  m_mapSet["constraintPdf"],
						  m_correlatedVar,
						  m_mapSet["allConstraint"],
						  process );
      break;
    }
    default : //if no constraint
      continue;
    }//end switch on constraint

    if (NPname.Contains("spurious") || NPname.Contains("BIAS") ) m_mapVar["spurious"]  = current_syst; // the systematics value is the spurious signal
    else if (type == "ESS" || type == "PES" ) m_mapSet["systematicPeak_"+string(process)]->add(*current_syst);
    else if (type == "MRES" || type == "PER" ) m_mapSet["systematicResolution_"+string(process)]->add(*current_syst);
    else { // means that type == YIELD
      m_mapSet[string("systematicYield_"+process)]->add(*current_syst);    
    }
    cout << "current syst" << endl;
    current_syst->Print();  
  } // end loop on systematics 


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
    RooRealVar *one = new RooRealVar("one_"+suffix, "one_"+suffix, central_val);
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

    //  systematic->Print();
  
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
    RooRealVar *one = new RooRealVar("one_"+suffix, "one_"+suffix, 1);
    RooProduct *systematic = new RooProduct("systematic_"+suffix, "systematic_"+suffix, RooArgSet(*one, *special_value)); // "1+" contained in vals_up and vals_down !
    //  systematic->Print();

    return (RooRealVar*) systematic;

  }

//#################################################
void Category::CreateBackgroundModel() {
  cout << "define the background model : " << endl;
  int bkg_model =  ((RooRealVar*)m_mapSet["systematicValues"]->find("bkg_model"))->getVal();
  
  switch (bkg_model ) {
  case  BKG_MODEL_EXPO : {
    RooRealVar *slope = new RooRealVar("slope", "slope", -2e-2, -5, -5e-6);
    m_mapSet["nuisanceParameters"]->add(*slope);
    RooExponential *exp_pdf = new RooExponential ("bkgExp", "bkgExp", *m_mapVar["invMass"], *slope);
    m_mapPdf["bkg"] = exp_pdf;
    break;
  }
  case  BKG_MODEL_EXPO_POL2 : {
    RooRealVar *p1 = new RooRealVar("p0", "p1", 1,   -100, 100);
    RooRealVar *p2 = new RooRealVar("p1", "p2", 1,   -100, 100);
    RooArgSet *coefficients = new RooArgSet(*p1,*p2);
    m_mapSet["nuisanceParameters"]->add(*coefficients);
    //      RooFormulaVar *x = new RooFormulaVar("x","x", "(@0-100.)/100.", invMass);
    RooGenericPdf *exp_pol2 = new RooGenericPdf("bkgExpPol2", "bkgExpPol2", "exp(@1*@0+@3*@0*@0)", RooArgList(*m_mapVar["mHRen"], *p1, *p2)); // quite long, WHY???
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
    m_mapSet["nuisanceParameters"]->add(*coefficients);
    RooBernstein *bern4 = new RooBernstein("bkgBern4", "bkgBern4", *m_mapVar["invMass"], *coefficients);
    m_mapPdf["bkg"] = bern4;
    break;
  }
  case  BKG_MODEL_BERN0 :   {
    RooRealVar *a0 = new RooRealVar("a0", "a0", 1);
    RooBernstein *bern0 = new RooBernstein("bkgBern0", "bkgBern0", *m_mapVar["invMass"], RooArgList(*a0));
    m_mapPdf["bkg"] = bern0;
    break;
  }
  }
  m_mapPdf["bkg"]->SetNameTitle( "bkg", "bkg" );
  m_mapSet["pdfToAdd"]->add( *m_mapPdf["bkg"] );
}


//##################################
void Category::CreateSignalModel() {

  for ( auto vProcess = m_processes->begin(); vProcess != m_processes->end(); vProcess++ ) {
    //Stores the gaussian and CB pdf of the process in this category   
    RooArgSet signalForms;
    for ( auto vForm = m_form.begin(); vForm != m_form.end(); vForm++ ) {
      map<string, RooProduct *> product;
      for ( auto vPar = m_param.begin(); vPar != m_param.end(); vPar++ ) {
	if ( ( *vPar == "yield"  && *vForm != "Var" )
	     || ( *vPar != "yield"  && *vForm == "Var" )
	     || (  *vPar == "alpha" && *vForm!= "CB" )
	     ) continue;
	//defined formulas for cental values and width of signal parametrization
	RooArgList varSet;
	varSet.add( *m_mapVar["mHcomb"] );
	varSet.add( *m_mapVar["mHRen"] );

	for ( auto vCoef = m_coef.begin(); vCoef != m_coef.end(); vCoef++ ) {
	  string dumString = string(TString::Format( "%s_%s%s", vCoef->c_str(), vPar->c_str(), vForm->c_str() ));
	  if ( m_mapVar[dumString] ) {
	    varSet.add( *m_mapVar[dumString] );
	  }
	}//end vCoef

	if ( varSet.getSize() == 2 ) continue;

	string formulaName = string(TString::Format( "%s%s_%s", vPar->c_str(), vForm->c_str(), vProcess->c_str() ));
	//Define the formula for the parameters
	//Polynome is defined according to varSet size
	string formula = "", base = "";
	for ( int iPlot=2; iPlot < varSet.getSize(); iPlot++ ) {
	  if ( formula != "" ) formula+="+";
	  formula+=TString::Format( "@%d%s%s", iPlot, base=="" ? "" : "*", base.c_str() );
	  base+=( base=="" ? "" : "*" ) + string("@1");
	}
	if ( *vPar == "mean" ) formula += "+@0";
	if ( *vPar == "yield" ) formula = string(TString::Format( "%2.2f*max(0.,%s)", m_mapVar["lumi"]->getVal(), formula.c_str() ) );
	m_mapFormula[formulaName] = new RooFormulaVar(formulaName.c_str(), formulaName.c_str(), formula.c_str(), varSet );

	//Contains mu/sigma and all systematics related
	RooArgSet factors;
	factors.add( *m_mapFormula[formulaName] );

	if ( *vPar == "sigma" ) {
	  factors.add( *m_mapSet["systematicResolution_" + *vProcess] );
	  factors.add( *m_mapSet["systematicResolution_common"] );
	}
	else if ( *vPar == "mean" ) {
	  factors.add( *m_mapSet["systematicPeak_" + *vProcess] );
	  factors.add( *m_mapSet["systematicPeak_common"] );
	}
	formulaName = "prod_" + formulaName;
	product[*vPar] = new RooProduct( formulaName.c_str(), formulaName.c_str(), factors );

      }//end vPar

      if ( ( *vForm == "GA" || *vForm == "CB" ) && ( !product["mean"] || !product["sigma"] ) ) continue;
      
      string pdfName = string( TString::Format( "%s_%s", vForm->c_str(), vProcess->c_str() ) );
      if ( *vForm == "CB" ) {
	string alphaCBName = string(TString::Format("alphaCB_%s", vProcess->c_str() ) );
	string nCBName = string(TString::Format("nCB_%s", vProcess->c_str()) );
	if ( !m_mapFormula[alphaCBName] || !m_mapVar[nCBName] ) continue;
	m_mapPdf[pdfName] =  new RooCBShape( pdfName.c_str(), pdfName.c_str(), *m_mapVar["invMass"],*product["mean"],*product["sigma"],*m_mapFormula[alphaCBName],*m_mapVar[nCBName] );
      }
      else if ( *vForm == "GA" ) m_mapPdf[pdfName] = new RooGaussian( pdfName.c_str(), pdfName.c_str(), *m_mapVar["invMass"], *product["mean"], *product["sigma"] );
      
      if ( *vForm == "GA" || *vForm == "CB" )  signalForms.add( *m_mapPdf[pdfName] );
      else {
	if ( signalForms.getSize() != 2 || ! product["yield"] ) continue;
	//Create a rooArglist with all variables to multiply yields with
	RooArgList yieldsFactors;
	yieldsFactors.add( *product["yield"] );
	yieldsFactors.add( *m_mapVar["mu_XS_" + *vProcess ] );//muXS
	yieldsFactors.add( *m_mapVar["mu"] );//globalMu
	yieldsFactors.add( *m_mapVar["mu_BR_yy"] );//muBR
	yieldsFactors.add( *m_mapSet["systematicYield_" + *vProcess] );
	yieldsFactors.add( *m_mapSet["systematicYield_common"] );
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


  switch ( m_signalModel ) {
  case 1 :
    m_mapPdf["signalSumPdf"] = new RooAddPdf( "signalSumPdf", "signalSumPdf", *m_mapSet["pdfProc"], *m_mapSet["yieldsToAdd"] );
    for ( unsigned int iProc = 0; iProc < m_processes->size(); iProc++ ) 
      m_mapSet["pdfToAdd"]->add( *m_mapPdf["signalSumPdf"] );
    break;
  default :
    m_mapSet["pdfToAdd"]->add( *m_mapSet["pdfProc"] );
    break;

}
}

//===================================
void Category::CreateSpurious() {

  if ( !m_mapVar["spurious"] ) return;

  RooAddPdf *spuriousPdf = new RooAddPdf( "spuriousPdf", "spuriousPdf", *m_mapSet["pdfProc"], *m_mapSet["yieldsSpurious"] );
  m_mapSet["pdfToAdd"]->add( *spuriousPdf );
  m_mapSet["yieldsToAdd"]->add( *m_mapVar["spurious"] );

}

//====================================
void Category::CreateWS() {

  ReadNuisanceParameters();
  CreateSignalModel();
  CreateSpurious();
  CreateBackgroundModel();

  m_mapPdf["modelSB"] =  new RooAddPdf( "modelSB", "modelSB", *m_mapSet["pdfToAdd"], *m_mapSet["yieldsToAdd"] );
  RooArgSet prodPdf;
  prodPdf.add( *m_mapPdf["modelSB"] );
  prodPdf.add( *m_mapSet["constraintPdf"] );
  m_mapPdf["model"] = new RooProdPdf( "model", "model", prodPdf );

  GetData();

  cout << "Correlated NP : " << m_correlatedVar << endl;
  string dumName = "ws_" + m_name;
  RooWorkspace *ws = new RooWorkspace( dumName.c_str(), dumName.c_str() );
  //  ws->import( m_mapPdf["model"], RecycleConflictNodes(), RenameAllVarExcept(m_correlatedVar), ReNameAllVar( m_name ) );
  ws->import( *m_dataset );

  vector<string> sets = { "obserables" };
  for ( auto set = sets.begin(); set != sets.end(); set++ ) DefineSet( *set );
}

//=========================================
 void Category::GetData() {

   TFile *inFile=0;
   if ( TString(m_dataFileName).Contains( ".txt" ) ) m_dataset = RooDataSet::read(m_dataFileName.c_str(), *m_mapSet["observables"]);
  else {

    TTree *inTree=0;
    RooWorkspace *inWS = 0;
    TObjArray *dataNomenclature = TString(m_dataFileName).Tokenize(" ");
    
    //Get the root file
    inFile = new TFile( dynamic_cast<TObjString*>(dataNomenclature->At(0))->GetString() );
    if ( !inFile ) {
      cout << dynamic_cast<TObjString*>(dataNomenclature->At(0))->GetString() << "does not exists" << endl;
      exit(0);
    }


    //Get the TTree or the workspace
    if ( dataNomenclature->GetEntries()>1 ) {
      TString objName = dynamic_cast<TObjString*>(dataNomenclature->At(1))->GetString();
      TString className =  inFile->Get( objName )->ClassName();
      if (  className == "TTree" ) {
	inTree = (TTree*) inFile->Get( objName );
	inTree->Print();
      }
      else {
	inWS = (RooWorkspace*) inFile->Get( objName );
	if ( !inWS ) {
	  cout << objName << " not found in " << inFile->GetName() << endl;
	  exit(0);
	}

	gROOT->cd();
	RooDataSet *dumDataSet=0;

	if ( dataNomenclature->GetEntries() > 2 ) dumDataSet = (RooDataSet*) inWS->data( dynamic_cast<TObjString*>(dataNomenclature->At(2))->GetString() );
	if ( dataNomenclature->GetEntries() > 3 ) m_mapVar["invMass"]->SetName( dynamic_cast<TObjString*>(dataNomenclature->At(3))->GetString() );
	m_mapSet["observables"]->add( *m_mapVar["invMass"] );
	//	observables->add( *m_eventCateg );
	if ( dumDataSet ) m_dataset = new RooDataSet( "obsData", "obsData",  
						    dumDataSet, 
						    *m_mapSet["observables"],
						    m_dataCut.c_str() );
	delete dumDataSet; dumDataSet=0;
	//	obsdata->Print();
      }//end else workspace
    }//end if tokenize

    if ( !m_dataset ) {
      cout << "No data found" << endl;
      exit(0);
    }
    dataNomenclature->Delete();

    delete inWS; inWS=0;
  }//end not txt file

   delete inFile; inFile=0;
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
    RooRealVar *nui = ((RooRealVar*)nuisance_parameters->find("nui_"+name));
    if (!nui) { // if the np does not exist yet in this channel
      RooRealVar* glob_nui = new RooRealVar("glob_nui_"+name, "glob_nui_"+name, 0, -5, 5);
      glob_nui->setConstant();
      nui = new RooRealVar("nui_"+name, "nui_"+name, 0, -5, 5);
      nuisance_parameters->add(*nui);
      global_parameters->add(*glob_nui);
      channel_correlated_np += string( ",nui_"+name+",glob_nui_"+name );
      // (Re)create the constraint
      RooRealVar *sigma_gauss_constraint = new RooRealVar("sigma_"+name, "sigma_"+name, 1);
      RooAbsPdf *constraint;
      if (sigmaRightBifurGauss>0) {
	RooRealVar *sigma_gauss_constraint_right = new RooRealVar("sigma_right_"+name, "sigma_right_"+name, sigmaRightBifurGauss);
	constraint = new RooBifurGauss("constraint_"+name, "constraint_"+name, *nui, *glob_nui, *sigma_gauss_constraint,  *sigma_gauss_constraint_right);
      }
      else {
	constraint = new RooGaussian("constraint_"+name, "constraint_"+name, *nui, *glob_nui, *sigma_gauss_constraint);
      }
      if ( name.Contains( "spurious" )  ) constraint->Print();
      if (! allConstraints->find(constraint->GetName())) { // if constraint not yet considered, should be the case since the NP does not exist yet
	constraints_pdf_list->add(*constraint); // add it to the list of constraints to be applied to this channel
	allConstraints->add(*constraint); // add it to the list of all constraints (for all channels)not to apply it twice

      }
    }

    return nui;
  }