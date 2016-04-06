#include "Workspace/Workspace.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "RooSimultaneous.h"
#include "RooDataSet.h"
using namespace RooFit;
#include <iostream>
using std::cout;
using std::endl;
#include <RooStats/ModelConfig.h>
using namespace RooStats;
using std::stringstream;
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TLatex.h"
#include "PlotFunctions/SideFunctions.h"
using std::ifstream;
#include "PlotFunctions/DrawPlot.h"

Workspace::Workspace() : m_debug(0)
{
  m_workspace=0;
  m_category = new RooCategory( "CombCategory", "combCategory" );
  m_mapSet["observables"] = new RooArgSet( "observables" );
  m_mapSet["globalObservables"] = new RooArgSet( "globalObservables" );
  m_mapSet["parametersOfInterest"] = new RooArgSet( "parametersOfInterest" );
  m_mapSet["nuisanceParameters"] = new RooArgSet( "nuisanceParameters" );
  m_name = "/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/Workspace.root";

}

Workspace::Workspace( string name ) : Workspace()
{
  m_name = name;
}

Workspace::~Workspace() {
  if ( m_workspace ) delete m_workspace;
  m_categories.clear();

}


//=======================================
void Workspace::Configure( string configFileName ) {
  m_configFileName = configFileName;

  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(configFileName, pt);
  string dum = pt.get<string>("General.catNames");
  ParseVector( dum  , m_categoriesNames );

  dum = pt.get<string>("General.process");
  ParseVector( dum, m_processes );

  m_systFileName = pt.get<string>("General.systFileName");

  dum = pt.get<string>("General.outName", "");
  if ( dum != "" ) m_name = dum;
}

//=======================================
void Workspace::CreateWS() {

  map<string,RooDataSet*> datasetMap;
  RooSimultaneous *pdf = new RooSimultaneous( "combinedPdf", "combinedPdf", *m_category );
  vector<string> sets = { "nuisanceParameters", "globalObservables", "observables", "parametersOfInterest" };

  readConstraintFile();

  for ( auto vName = m_categoriesNames.begin(); vName != m_categoriesNames.end(); vName++ ) {
    m_category->defineType( vName->c_str() );
    m_category->setLabel( vName->c_str() );
    m_categories.push_back( 0);
    m_categories.back() = new Category( *vName );
    m_categories.back()->SetSDef( &m_sDef );
    m_categories.back()->SetSystFileName( m_systFileName );
    m_categories.back()->SetProcesses( &m_processes );
    m_categories.back()->LoadParameters( m_configFileName );
    m_categories.back()->CreateWS();

    RooWorkspace *workspace = m_categories.back()->GetWorkspace();
    string pdfName = "model_" + *vName;
    pdf->addPdf( *workspace->pdf( pdfName.c_str() ), m_category->getLabel() );


    for ( auto vSet = sets.begin(); vSet != sets.end(); vSet++ )  {
      cout << *vSet << " " << m_mapSet[*vSet] << endl;
      m_mapSet[*vSet]->add( *workspace->set(vSet->c_str()) );
    }
    string dataName = "obsData_" + *vName;
    datasetMap[*vName] = (RooDataSet*) workspace->data( dataName.c_str() );

  }

  m_workspace = new RooWorkspace( "combination", "combination" );
  m_workspace->import( *pdf, RecycleConflictNodes() );
  cout << "importedPdf" << endl;
  for ( auto vSet = sets.begin(); vSet != sets.end(); vSet++ ) m_workspace->defineSet( vSet->c_str(), *m_mapSet[*vSet], kTRUE );


  RooDataSet* obsData = new RooDataSet("obsData","combined data ",*m_mapSet["observables"], Index(*m_category), Import(datasetMap)); // ,WeightVar(wt));
  m_workspace->import(*obsData);
  m_workspace->var( "mHcomb" )->setConstant(1);
  m_workspace->var( "mu" )->setConstant(0);



  cout << "Creating dataset with ghosts..." << endl;
  RooDataSet* newData = addGhosts(obsData,m_workspace->set("observables"));
  newData->SetNameTitle("obsData_G","obsData_G");
  m_workspace->import(*newData);
  cout << "... sucessfully imported." << endl;
  //  m_workspace->pdf("combinedPdf")->fitTo( *newData );

  cout << "Defining mconfig..." << endl; 
  ModelConfig *mconfig = new ModelConfig("mconfig", m_workspace);
  mconfig->SetPdf(*m_workspace->pdf("combinedPdf"));
  
  mconfig->SetObservables( *m_workspace->set("observables"));
  mconfig->SetParametersOfInterest( *m_workspace->set("parametersOfInterest") );
  mconfig->SetNuisanceParameters( *m_workspace->set("nuisanceParameters") );
  mconfig->SetGlobalObservables( *m_workspace->set("globalObservables") );


  m_workspace->import(*mconfig);

  cout << "Saving workspace to file... '" << m_name << "'" << endl;
  m_workspace->writeToFile(m_name.c_str(), 1);

}


//==================================================
RooDataSet* Workspace::addGhosts(RooDataSet* orig,  const RooArgSet *observables ) {
  if ( m_debug ) cout << "addGhosts" << endl;
  map<string, RooDataSet*> data_map;
  TList* datalist = orig->split(*m_category, true);
  TIterator* dataItr = datalist->MakeIterator();
  RooAbsData* ds;
  RooRealVar* weightVar = new RooRealVar("wt","wt",1);
  RooRealVar* firstObs = (RooRealVar*) observables->first();  
  RooArgSet obsWeight( *firstObs, *m_category, *weightVar );
  //  obsWeight.add( *weightVar );
  while ((ds = (RooAbsData*)dataItr->Next())) { // loop over all channels
    m_category->setLabel( ds->GetName() );
    int nrEntries = ds->numEntries();
    string typeName(ds->GetName());
    stringstream datasetName;
    datasetName << "newData_" << typeName << endl;

    RooDataSet* thisData = new RooDataSet(datasetName.str().c_str(),datasetName.str().c_str(), obsWeight, WeightVar(*weightVar) );
    //    cout << "entries : " << nrEntries << endl;

    for (int ib=0;ib<nrEntries;ib++) {
      firstObs->setVal( ((RooRealVar*) ds->get(ib)->first())->getVal() );
      thisData->add(obsWeight, 1);  
    }

    TString string_tn(typeName);
    if ( nrEntries <20) {
      for (double mgg =  firstObs->getMin(); mgg < firstObs->getMax(); mgg += 0.1) {
	firstObs->setVal( mgg );
	thisData->add(obsWeight, 1e-9);
	//	  Cout << "Adding evt to ds " << thisData->GetName() << " at mgg=" << mgg << endl;
      }
      //      int nrEntries2 = thisData->numEntries();
      //      cout << " nb entries 2 " << nrEntries2 << endl; 
    }

      if ( m_debug ) {
	TCanvas *c = new TCanvas();
	RooPlot* frame=firstObs->frame(40);
       frame->SetTitle(""); //empty title to prevent printing "A RooPlot of ..."
       frame->SetXTitle("m_{H}");
       char buffer_dummy[50];
       sprintf(buffer_dummy,"Events / %4.1f GeV",(frame->GetXaxis()->GetXmax()-frame->GetXaxis()->GetXmin())/frame->GetXaxis()->GetNbins());
       frame->SetXTitle("m_{#gamma#gamma} [GeV]");
       frame->SetYTitle(buffer_dummy);
       thisData->plotOn(frame);    
       frame->Draw();
       TLatex latex;
       latex.DrawLatexNDC( 0.15, 0.9, TString::Format( "sumentries : %2.2f ", thisData->sumEntries() ).Data() );
       latex.DrawLatexNDC( 0.15, 0.85, TString::Format( "entries : %d ", thisData->numEntries() ).Data() );
       c->SaveAs("Figures/ghost_"+TString(m_category->getLabel())+".pdf");
       delete frame;
       delete c;
     }

      //    cout << "NAME = " <<  ds->GetName() << endl;
    data_map[string(ds->GetName())] = (RooDataSet*)thisData;
  }

  RooDataSet* newData = new RooDataSet("newData","newData", obsWeight,
				       Index(*m_category), Import(data_map), WeightVar( *weightVar) );

  newData->Print();
  return newData;
  
}


//=================================
void Workspace::readConstraintFile()
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
	switch ( tmpAr.GetEntries() ) {
	case 1 : 
	  continue;
	case 2 : {
	  TString tmpStrDefConst= ((TObjString*) tmpAr.At(1))->GetString();
	  if(tmpStrDefConst == "NO_CONSTRAINT") defConstraint = NO_CONSTRAINT;
	  else if(tmpStrDefConst == "GAUSS_CONSTRAINT") defConstraint = GAUSS_CONSTRAINT;
	  else if(tmpStrDefConst == "LOGNORM_CONSTRAINT") defConstraint = LOGNORM_CONSTRAINT;
	  else if(tmpStrDefConst == "ASYM_CONSTRAINT") defConstraint = ASYM_CONSTRAINT;
	  else {
	    cout << "Unknown constraint for syst. " << ((TObjString*) tmpAr.First())->GetString() << endl;
	    continue;
	  }
	}// end case 2
	  break;
	case 8 :
	  defConstraint = ASYM_CONSTRAINT;
	  break;
	case 3 :
	  defConstraint = GAUSS_CONSTRAINT;
	  break;
	default :
	  defConstraint = NO_CONSTRAINT;
	}//end switch
	m_sDef[string(((TObjString*) tmpAr.First())->GetString())] = defConstraint;
      } while(!current_file.eof());


  }

