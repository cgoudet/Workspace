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


Workspace::Workspace() : m_debug(0)
{
  m_workspace=0;
  m_category = new RooCategory( "category", "category" );
  m_mapSet["observables"] = new RooArgSet( "observables" );
  m_mapSet["globalObservables"] = new RooArgSet( "globalObservables" );
  m_mapSet["parametersOfInterest"] = new RooArgSet( "parametersOfInterest" );
  m_mapSet["nuisanceParameters"] = new RooArgSet( "nuisanceParameters" );
  m_name = "Workspace";

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
  std::cout << pt.get<std::string>("Section1.Value1") << std::endl;
  std::cout << pt.get<std::string>("Section1.Value2") << std::endl;
}

//=======================================
void Workspace::CreateWS() {

  map<string,RooDataSet*> datasetMap;
  RooSimultaneous *pdf = new RooSimultaneous( "combinedPdf", "combinedPdf", *m_category );
  vector<string> sets = { "nuisanceParameters", "globalObservables", "observables", "parameterOfInterest" };

  for ( auto vName = m_categoriesNames.begin(); vName != m_categoriesNames.end(); vName++ ) {
    m_category->defineType( vName->c_str() );
    m_categories.push_back( Category( *vName ) );
    m_categories.back().SetSDef( &m_sDef );
    m_categories.back().SetSystFileName( m_systFileName );
    m_categories.back().SetProcesses( &m_processes );
    m_categories.back().LoadParameters( m_configFileName );
    m_categories.back().CreateWS();

    RooWorkspace *workspace = m_categories.back().GetWorkspace();
    pdf->addPdf( *workspace->pdf( string("model_" + *vName).c_str() ), m_category->getLabel() );


    for ( auto vSet = sets.begin(); vSet != sets.end(); vSet++ ) 
      m_mapSet[*vSet]->add( *workspace->set(vSet->c_str()) );

    datasetMap[*vName] = (RooDataSet*) workspace->data("obsData");
  }

  m_workspace = new RooWorkspace( "combination", "combination" );
  for ( auto vSet = sets.begin(); vSet != sets.end(); vSet++ ) 
    m_workspace->defineSet( vSet->c_str(), *m_mapSet[*vSet], kTRUE );

  RooDataSet* obsData = new RooDataSet("obsData","combined data ",*m_mapSet["observables"], Index(*m_category), Import(datasetMap)); // ,WeightVar(wt));
  m_workspace->import(*obsData);
  
  cout << "Creating dataset with ghosts..." << endl;
  RooDataSet* newData = addGhosts(obsData,m_workspace->set("observables"));
  newData->SetNameTitle("obsData_G","obsData_G");
  m_workspace->import(*newData);
  cout << "... sucessfully imported." << endl;

  cout << "Defining mconfig..." << endl; 
  ModelConfig *mconfig = new ModelConfig("mconfig", m_workspace);
  // mconfig->SetPdf(*m_workspace->pdf("combinedPdf"));

  // mconfig->SetObservables( *m_workspace->set("Observables"));
  // mconfig->SetParametersOfInterest( *m_workspace->set("ParametersOfInterest") );
  // mconfig->SetNuisanceParameters( *m_workspace->set("nuisanceParameters") );
  // mconfig->SetGlobalObservables( *m_workspace->set("globalObservables") );


  m_workspace->import(*mconfig);

  TString myfile = "/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/"+m_name+".root";
  cout << "Saving workspace to file... '" << myfile << "'" << endl;
  m_workspace->writeToFile(myfile, 1);

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
