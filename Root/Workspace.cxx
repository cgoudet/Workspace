#include "Workspace/Workspace.h"
#include "PlotFunctions/SideFunctions.h"
#include "PlotFunctions/SideFunctionsTpp.h"
#include "PlotFunctions/DrawPlot.h"
#include "PlotFunctions/Arbre.h"
using namespace ChrisLib;

#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include <RooStats/ModelConfig.h>
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooFitResult.h"
#include "TDOMParser.h"
#include "TXMLDocument.h"
#include "TXMLNode.h"
#include "TXMLAttr.h"
#include "TIterator.h"
using namespace RooFit;
using namespace RooStats;

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <exception>
using std::runtime_error;
#include <iostream>
using std::cout;
using std::endl;
using std::stringstream;
using std::ifstream;
using std::list;

using std::string;
using std::vector;
using std::map;

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
  if ( m_debug ) cout << "Workspace::Configure" << endl;
  m_configFileName = configFileName;

  Arbre arbre( Arbre::ParseXML(configFileName) );
  m_name = arbre.GetAttribute( "Name" );

  list<Arbre> children = arbre.GetChildren();
  for ( list<Arbre>::const_iterator it=children.begin(); it!=children.end(); ++it ) {
    string nodeName = it->GetNodeName();
    if ( nodeName == "processes" ) ParseVector( it->GetAttribute( "text" ), m_processes );
    else if ( nodeName == "category" ) m_categoriesNames.push_back( it->GetAttribute("Name") );
  }
  if ( m_debug ) cout << "Workspace::Configure Done" << endl;
}

//=======================================
void Workspace::CreateWS() {

  map<string,RooDataSet*> datasetMap;
  RooSimultaneous *pdf = new RooSimultaneous( "combinedPdf", "combinedPdf", *m_category );
  vector<string> sets = { "nuisanceParameters", "globalObservables", "observables", "parametersOfInterest" };

  for ( auto vName : m_categoriesNames ) {
    m_category->defineType( vName.c_str() );
    m_category->setLabel( vName.c_str() );
    m_categories.push_back( 0);
    m_categories.back() = new Category( vName );
    if ( m_debug ) m_categories.back()->SetDebug( m_debug );
    m_categories.back()->SetProcesses( &m_processes );
    m_categories.back()->SetCategoriesNames( &m_categoriesNames );
    m_categories.back()->LoadParameters( m_configFileName );
    m_categories.back()->CreateWS();

    RooWorkspace *workspace = m_categories.back()->GetWorkspace();
    string pdfName = "model_" + vName;
    workspace->pdf( pdfName.c_str() )->Print();
    pdf->addPdf( *workspace->pdf( pdfName.c_str() ), m_category->getLabel() );


    for ( auto vSet : sets ) m_mapSet[vSet]->add( *workspace->set(vSet.c_str()), 1 );
    string dataName = "obsData_" + vName;
    datasetMap[vName] = (RooDataSet*) workspace->data( dataName.c_str() );
    if ( !datasetMap[vName] ) { cout << "Null dataset : " << vName << endl; exit(0);}
  }

  m_workspace = new RooWorkspace( "combination", "combination" );

  m_workspace->import( *pdf, RecycleConflictNodes(), Silence() );
  cout << "importedPdf" << endl;
  for ( auto vSet  : sets ) m_workspace->defineSet( vSet.c_str(), *m_mapSet[vSet], kTRUE );
  RooRealVar wt( "weight", "weight", 1 );
  m_mapSet["observables"]->add( wt, 1 );

  RooDataSet* obsData = new RooDataSet("obsData","combined data ",*m_mapSet["observables"], Index(*m_category), Import(datasetMap) ,WeightVar(wt ));
  cout << "obsData : " << obsData << endl;
  cout << obsData->GetName() << endl;
  m_workspace->import(*obsData, Silence());

  cout << "Creating dataset with ghosts..." << endl;
  // RooDataSet* newData = addGhosts(obsData,m_workspace->set("observables"));
  // newData->SetNameTitle("obsData_G","obsData_G");
  // m_workspace->import(*newData, Silence());
  // cout << "... sucessfully imported." << endl;
  
  cout << "fitted constrained" << endl;


  cout << "Defining mconfig..." << endl; 
  ModelConfig *mconfig = new ModelConfig("mconfig", m_workspace);
  mconfig->SetPdf(*m_workspace->pdf("combinedPdf"));
  
  mconfig->SetObservables( *m_workspace->set("observables"));
  mconfig->SetParametersOfInterest( *m_workspace->set("parametersOfInterest") );
  mconfig->SetNuisanceParameters( *m_workspace->set("nuisanceParameters") );
  mconfig->SetGlobalObservables( *m_workspace->set("globalObservables") );



  m_workspace->import(*mconfig);
  RooDataSet* asimovData = MakeAsimovData();
  cout << "Saving workspace to file... '" << m_name << "'" << endl;
  m_workspace->importClassCode();
  m_workspace->writeToFile(m_name.c_str(), 1);

  m_workspace->Print();
  m_workspace->obj("yield_common_Inclusive")->Print();
  //  mconfig->GetPdf()->fitTo( *m_workspace->data( obsData->GetName() ), SumW2Error(kFALSE) );
  // RooArgSet dumSet;
  // dumSet.add( *mconfig->GetNuisanceParameters() );
  // dumSet.add( *mconfig->GetParametersOfInterest() );
  // m_workspace->saveSnapshot("miminum",dumSet);
}


//==================================================
RooDataSet* Workspace::addGhosts(RooDataSet* orig,  const RooArgSet *observables ) {
  if ( m_debug ) cout << "addGhosts" << endl;

  map<string, RooDataSet*> data_map;
  TList* datalist = orig->split(*m_category, true);
  TIterator* dataItr = datalist->MakeIterator();
  RooAbsData* ds;
  RooRealVar* weightVar = new RooRealVar("weight","weight",1);
  RooArgSet obsWeight;

  while ((ds = (RooAbsData*)dataItr->Next())) { // loop over all channels
    m_category->setLabel( ds->GetName() );
    int nrEntries = ds->numEntries();
    string typeName(ds->GetName());
    stringstream datasetName;
    datasetName << "newData_" << typeName << endl;

    RooRealVar* firstObs = 0;
    TIterator* iter = observables->createIterator();
    RooAbsPdf* parg;
    while((parg=(RooAbsPdf*)iter->Next()) ) {
      if ( TString( parg->GetName() ).Contains( ds->GetName() ) ) {
	firstObs = (RooRealVar*) parg;
	break;
      }
    }
    RooDataSet* thisData = new RooDataSet( *(RooDataSet*) ds, datasetName.str().c_str() );
    obsWeight = RooArgSet( *firstObs, *m_category, *weightVar );
	
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

  obsWeight = RooRealVar();
  obsWeight.add( *weightVar );
  obsWeight.add( *orig->get(1) );
  RooDataSet* newData = new RooDataSet("newData","newData", obsWeight,
				       Index(*m_category), Import(data_map), WeightVar( *weightVar) );

  newData->Print();
  return newData;
  
}


//=================================
RooDataSet* Workspace::MakeAsimovData() {
  
  cout << " Make Asimov Data Beginning " << endl; 
  unsigned int bkgOnly=0;
  unsigned int doConditional = 0;
  
  // ////////////////////
  // //make asimov data//
  // ////////////////////
  
  stringstream muStr;
  muStr << "_" << !bkgOnly;

  //Fix mus to some value
  for ( auto vProc : m_processes ) {
    string muName = "mu_XS_"+vProc;
    RooRealVar * mu = m_workspace->var( muName.c_str() );
    if ( !mu ) continue;
    mu->setConstant(0);
    mu->setVal( bkgOnly ? 0 : 1 );
  }

  RooRealVar* mu = (RooRealVar*) m_workspace->var("mu");  
  RooArgSet mc_obs = *m_workspace->set("observables");
  RooArgSet mc_globs = *m_workspace->set("globalObservables");
  RooArgSet mc_nuis = *m_workspace->set("nuisanceParameters");
  
  // set all global variables to constant
  //const RooArgSet* globs_ptr = mcInWs->GetGlobalObservables();
  TIterator* glob_itr = mc_globs.createIterator();
  RooRealVar* glob_var;
  while ((glob_var = (RooRealVar*)glob_itr->Next())) glob_var->setConstant(1);

  RooSimultaneous *combPdf = (RooSimultaneous*) m_workspace->pdf( "combinedPdf" );
  //pair the nuisance parameter to the global observable
  RooArgSet mc_nuis_tmp = mc_nuis;
  RooArgList nui_list("ordered_nuis");
  RooArgList glob_list("ordered_globs");
  RooArgSet constraint_set_tmp(*combPdf->getAllConstraints(mc_obs, mc_nuis_tmp, false));
  RooArgSet constraint_set;
  int counter_tmp = 0;
  UnfoldConstraints(constraint_set_tmp, constraint_set, mc_obs, mc_nuis_tmp, counter_tmp);
  cout << " After Unfold Constraints " << endl; 

  TIterator* cIter = constraint_set.createIterator();
  RooAbsArg* arg;
  //For each constraint
  while ((arg = (RooAbsArg*)cIter->Next())) {
    RooAbsPdf* pdf = (RooAbsPdf*)arg;
    if (!pdf) continue;
    //finds a variable from which the pdf depends on
    TIterator* nIter = mc_nuis.createIterator();
    RooRealVar* thisNui = NULL;
    RooAbsArg* nui_arg;
    while ((nui_arg = (RooAbsArg*)nIter->Next())) {
      if (pdf->dependsOn(*nui_arg)) {
  	thisNui = (RooRealVar*)nui_arg;
  	break;
      }
    }
    delete nIter;

    //Delete components that depends on other variables in the set
    RooArgSet* components = pdf->getComponents();
    components->remove(*pdf);
    if (components->getSize()) {
      TIterator* itr1 = components->createIterator();
      RooAbsArg* arg1;
      while ((arg1 = (RooAbsArg*)itr1->Next())) {
  	TIterator* itr2 = components->createIterator();
  	RooAbsArg* arg2;
  	while ((arg2 = (RooAbsArg*)itr2->Next())) {
  	  if (arg1 == arg2) 
  	    continue;
  	  if (arg2->dependsOn(*arg1)) {
	    components->remove(*arg1);
  	  }
  	}
  	delete itr2;
      }
      delete itr1;
    }//end components getsize

    if (components->getSize() > 1) {
      cout << "ERROR::Couldn't isolate proper nuisance parameter" << endl;
      return 0;
    }
    else if (components->getSize() == 1) {
      thisNui = (RooRealVar*)components->first();
    }
    
    //Get the global observable corresponding on which the constraint depend.    
    TIterator* gIter = mc_globs.createIterator();
    RooRealVar* thisGlob = NULL;
    RooAbsArg* glob_arg;
    while ((glob_arg = (RooAbsArg*)gIter->Next())) {
      if (pdf->dependsOn(*glob_arg)) {
  	thisGlob = (RooRealVar*)glob_arg;
  	break;
      }
    }
    delete gIter;

    if (!thisNui || !thisGlob) {
      cout << "WARNING::Couldn't find nui or glob for constraint: " << pdf->GetName() << endl;
      continue;
    }
    
    cout << "Pairing nui: " << thisNui->GetName() << ", with glob: " << thisGlob->GetName() << ", from constraint: " << pdf->GetName() << endl;
    
    nui_list.add(*thisNui);
    glob_list.add(*thisGlob);  
  } // end loop on cIter (constraints)
  delete cIter;
  
  
  m_workspace->saveSnapshot("nominalGlobs",glob_list);
  m_workspace->saveSnapshot("nominalNuis", nui_list);
  

  //Set the spurious signal to 0
  TIterator *iter = nui_list.createIterator();
  RooRealVar *dumVar =0;
  while ((dumVar = (RooRealVar*)iter->Next())) {
    TString nui_name = dumVar->GetName();
    if ( !nui_name.Contains("spurious_") && !nui_name.Contains("BIAS")) continue;
    cout << " Fixing to constant nuisance parameters" << endl; 
    dumVar->setVal(0);
    dumVar->setConstant(1);
  }
  delete iter; iter=0;
  
  //fit data with backgroud only model.
  // RooArgSet nuiSet_tmp(nui_list);
  // RooDataSet *combData = (RooDataSet*) m_workspace->data( "obsData_G" );
  RooFitResult *result2 = 0;//combPdf->fitTo(*combData,Hesse(false),Minos(false),PrintLevel(0),Extended(), Constrain(nuiSet_tmp), Save(),SumW2Error(kFALSE));
  //  combPdf->fitTo(*combData, Hesse(false),Minos(false),PrintLevel(0),Extended(), Constrain(nuiSet_tmp));
  //  cout << "Done" << endl;
  //  mu->setConstant(0);

  int nrNuis = nui_list.getSize();   
  if (nrNuis != glob_list.getSize())  {
    cout << "ERROR::nui_list.getSize() != glob_list.getSize()!" << endl;
    return 0;
  }


  //If randomization of the globals do it here
  for (int i=0;i<nrNuis;i++) {
    RooRealVar* nui = (RooRealVar*)nui_list.at(i);
    RooRealVar* glob = (RooRealVar*)glob_list.at(i);
    cout << "nui: " << nui << ", glob: " << glob << endl;
    cout << "Setting glob: " << glob->GetName() << ", which had previous val: " << glob->getVal() << ", to conditional val: " << nui->getVal() << endl;
    glob->setVal(nui->getVal());  
  }//end for i


  //save the snapshots of conditional parameters
  cout << "Saving conditional snapshots" << endl;
  m_workspace->saveSnapshot(("conditionalGlobs"+muStr.str()).c_str(),glob_list);
  m_workspace->saveSnapshot(("conditionalNuis" +muStr.str()).c_str(), nui_list);
  
  if (!doConditional) {
    m_workspace->loadSnapshot("nominalGlobs");
    m_workspace->loadSnapshot("nominalNuis");
  }
  
  
  iter = nui_list.createIterator();
  dumVar =0;
  while ((dumVar = (RooRealVar*)iter->Next())) {
    TString nui_name = dumVar->GetName();
    if ( !nui_name.Contains("spurious_") && !nui_name.Contains("BIAS")) continue;
    dumVar->setConstant(0);
  }
  delete iter; iter=0;

  
  cout << "Making asimov" << endl;
  ModelConfig* mc = static_cast<ModelConfig*>( m_workspace->obj( "mconfig" ));
  if ( !mc ) throw runtime_error( "Workspace::MakeAsimovData : ModelConfig not found" );


  int iFrame=0;
  const char* weightName="weightVar";
  RooArgSet obsAndWeight;
  obsAndWeight.add(*mc->GetObservables());
  
  RooRealVar* weightVar = NULL;
  if (!(weightVar = m_workspace->var(weightName))) {
    m_workspace->import(*(new RooRealVar(weightName, weightName, 1,0,10000000)));
    //   m_workspace->import(*(new RooRealVar(weightName, weightName, 1,0,1000)));
    weightVar = m_workspace->var(weightName);    
  }
  //  cout << "weightVar: " << weightVar << endl;
  obsAndWeight.add(*m_workspace->var(weightName));
  
  //  cout << "defining set" << endl;
  m_workspace->defineSet("obsAndWeight",obsAndWeight);

  
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  // MAKE ASIMOV DATA FOR OBSERVABLES

  // dummy var can just have one bin since it's a dummy
  // if(m_workspace->var("dummyX"))  
  //   m_workspace->var("dummyX")->setBins(1);

  cout << endl << "Check expectedData by category" << endl;
  //RooDataSet* simData=NULL;
  RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
  map<string, RooDataSet*> asimovDataMap;
  //try fix for sim pdf
  RooCategory* channelCat = (RooCategory*) (&simPdf->indexCat());

  // TIterator* iter = simPdf->indexCat().typeIterator() ;
  // TIterator* iter = channelCat->typeIterator() ;
  // RooCatType* tt = NULL;
  int nrIndices = channelCat->numTypes();

  for (int i=0;i<nrIndices;i++) {
    channelCat->setIndex(i);
    iFrame++;
    // Get pdf associated with state from simpdf
    RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;
	
    // Generate observables defined by the pdf associated with this state
    RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;
    obstmp->Print();
    
    //channelCat->setIndex(string(tt->GetName()).c_str());
    cout << "on type " << channelCat->getLabel() << " " << iFrame << endl;
    
    RooDataSet* obsDataUnbinned = new RooDataSet(Form("combAsimovData%d",iFrame),Form("combAsimovData%d",iFrame),RooArgSet(obsAndWeight,*channelCat),WeightVar(*weightVar));
    RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
    thisObs->setBins( (thisObs->getMax() - thisObs->getMin()) /0.1);
    cout << "thisObsBins : " << thisObs->getBins() << " " << thisObs->getMax() << " " << thisObs->getMin() << endl;
    thisObs->Print();
    double expectedEvents = pdftmp->expectedEvents(*obstmp);
    cout << "expectedEvents : " << expectedEvents << endl;

    double thisNorm = 0;
    cout << "nBins : " << thisObs->numBins() << endl;
    for(int jj=0; jj<thisObs->numBins(); ++jj) {
      thisObs->setBin(jj);
      //      thisObs->Print();
      //	cout << "pdf = "<<pdftmp->getVal(obstmp) <<endl;
      thisNorm = pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
      //      cout << "thisNorm : " << thisNorm*expectedEvents << endl;
      //      if (thisNorm*expectedEvents > pow(10., -2) && thisNorm*expectedEvents < pow(10., 9)) 
      obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
    }
    
    obsDataUnbinned->Print();
    cout << "sum entries " << obsDataUnbinned->sumEntries() << endl;

    if(obsDataUnbinned->sumEntries() != obsDataUnbinned->sumEntries()){
      cout << "sum entries is nan"<<endl;
      exit(1);
    }
    
    ((RooRealVar*)obstmp->first())->Print();
    cout << "expected events " << pdftmp->expectedEvents(*obstmp) << endl;
    
    asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;//tempData;
    
    cout << "channel: " << channelCat->getLabel() << ", data: ";
    obsDataUnbinned->Print();
    cout << endl;
  }
  
  
  RooDataSet* asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));

  if (m_workspace->data(asimovData->GetName())) {
    m_workspace->import(*asimovData, true ); // Bool_t import(TObject& object, const char* aliasName, Bool_t replaceExisting = kFALSE)
  }
  else {
    m_workspace->import(*asimovData );
  }
  
  
  //  RooDataSet* asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));
  //m_workspace->import(*asimovData);
  
  
  //bring us back to nominal for exporting
  m_workspace->loadSnapshot("nominalNuis");
  m_workspace->loadSnapshot("nominalGlobs");
  
  asimovData->Print();
  if ( result2 )  delete result2;  
  return asimovData;
}

//=============================================================
void Workspace::UnfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter)
{
  if (counter > 50) {
    cout << "ERROR::Couldn't unfold constraints!" << endl;
    cout << "Initial: " << endl;
    initial.Print("v");
    cout << endl;
    cout << "Final: " << endl;
    final.Print("v");
    exit(1);
  }

  TIterator* itr = initial.createIterator();
  RooAbsPdf* pdf;
  while ((pdf = (RooAbsPdf*)itr->Next()))
    {
      RooArgSet nuis_tmp = nuis;
      RooArgSet constraint_set(*pdf->getAllConstraints(obs, nuis_tmp, false));
      string className(pdf->ClassName());
      if (className != "RooGaussian" && className != "RooLognormal" && className != "RooGamma" && className != "RooPoisson" && className != "RooBifurGauss")
	{
	  counter++;
	  UnfoldConstraints(constraint_set, final, obs, nuis, counter);
	}
      else
	{
	  final.add(*pdf);
	}
    }
  delete itr;
}
