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

    cout << "fit 2 : " << endl;
    //    workspace->pdf( pdfName.c_str() )->fitTo( *workspace->data( dataName.c_str() ) );
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

  //  mconfig->GetNuisanceParameters()->Print("v");
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

// void Workspace::makeAsimovData( RooRealVar* mH, ModelConfig* mcInWs, bool doConditional, RooWorkspace* combWs, RooAbsPdf* combPdf, RooDataSet* combData, bool b_only)
// {
  // cout << " Make Asimov Data Beginning " << endl; 
  
  // ////////////////////
  // //make asimov data//
  // ////////////////////
  
  // stringstream muStr;
  // muStr << "_" << !b_only;
  // RooRealVar* mu = (RooRealVar*)mcInWs->GetParametersOfInterest()->first();
  // mu->setVal(!b_only);  

  // cout <<  mu->getVal() << endl; 


  // // RooRealVar* mH = NULL; 
  // // mH = combWs->var("mH");
  // mH->setConstant();
  // mH->setVal(125.6);

  // RooRealVar* mu_ggF = NULL; 
  // mu_ggF = combWs->var("mu_ggH");
  // mu_ggF->setConstant();
  // mu_ggF->setVal(1);

  // RooRealVar* mu_VBF = NULL; 
  // mu_VBF= combWs->var("mu_VBF");
  // mu_VBF->setConstant();
  // mu_VBF->setVal(1);

  // cout << "mu ggF Asimov == " << mu_ggF->getVal() << endl; 
  // cout << "mu VBF Asimov == " << mu_VBF->getVal() << endl; 

  // RooRealVar* mu_VH = NULL; 
  // if (combWs->var("mu_VH")) {
  //   mu_VH = combWs->var("mu_VH");
  //   mu_VH->setConstant();
  //   mu_VH->setVal(1);
  //   cout << "mu VH Asimov == " << mu_VH->getVal() << endl; 
  // }

  // RooRealVar* mu_WH = NULL; 
  // if (combWs->var("mu_WH")) {
  //   mu_WH = combWs->var("mu_WH");
  //   mu_WH->setConstant();
  //   mu_WH->setVal(1);
  //   cout << "mu WH Asimov == " << mu_WH->getVal() << endl; 
  // }

  // RooRealVar* mu_ZH = NULL; 
  // if (combWs->var("mu_ZH")) {
  //   mu_ZH = combWs->var("mu_ZH");
  //   mu_ZH->setConstant();
  //   mu_ZH->setVal(1);
  //   cout << "mu ZH Asimov == " << mu_ZH->getVal() << endl; 
  // }

  // RooRealVar* mu_ttH = NULL; 
  // if (combWs->var("mu_ttH")) {
  //   mu_ttH = combWs->var("mu_ttH");
  //   mu_ttH->setConstant();
  //   mu_ttH->setVal(1);
  //   cout << "mu ttH Asimov == " << mu_ttH->getVal() << endl; 
  // }

  // RooRealVar* mu_bbH = NULL; 
  // if (combWs->var("mu_bbH")) {
  //   mu_bbH = combWs->var("mu_bbH");
  //   mu_bbH->setConstant();
  //   mu_bbH->setVal(1);
  //   cout << "mu bbH Asimov == " << mu_bbH->getVal() << endl; 
  // }
  
  // RooArgSet mc_obs = *mcInWs->GetObservables();
  // RooArgSet mc_globs = *mcInWs->GetGlobalObservables();
  // RooArgSet mc_nuis = *mcInWs->GetNuisanceParameters();
  
  // // set all global variables to constant
  // const RooArgSet* globs_ptr = mcInWs->GetGlobalObservables();
  // TIterator* glob_itr = globs_ptr->createIterator();
  // RooRealVar* glob_var;
  // while ((glob_var = (RooRealVar*)glob_itr->Next())) {
  //   glob_var->setConstant(1);
  // }
  
  // //pair the nuisance parameter to the global observable
  // RooArgSet mc_nuis_tmp = mc_nuis;
  // RooArgList nui_list("ordered_nuis");
  // RooArgList glob_list("ordered_globs");
  // RooArgSet constraint_set_tmp(*combPdf->getAllConstraints(mc_obs, mc_nuis_tmp, false));
  // RooArgSet constraint_set;
  // int counter_tmp = 0;
  // unfoldConstraints(constraint_set_tmp, constraint_set, mc_obs, mc_nuis_tmp, counter_tmp);
  // cout << " After Unfold Constraints " << endl; 

  // TIterator* cIter = constraint_set.createIterator();
  // RooAbsArg* arg;
  // //For each constraint
  // while ((arg = (RooAbsArg*)cIter->Next())) {
  //   RooAbsPdf* pdf = (RooAbsPdf*)arg;
  //   if (!pdf) continue;

  //   //finds a variable from which the pdf depends on
  //   TIterator* nIter = mc_nuis.createIterator();
  //   RooRealVar* thisNui = NULL;
  //   RooAbsArg* nui_arg;
  //   while ((nui_arg = (RooAbsArg*)nIter->Next())) {
  //     if (pdf->dependsOn(*nui_arg)) {
  // 	thisNui = (RooRealVar*)nui_arg;
  // 	break;
  //     }
  //   }
  //   delete nIter;

  //   //Delete components that depends on other variables in the set
  //   RooArgSet* components = pdf->getComponents();
  //   components->remove(*pdf);
  //   if (components->getSize()) {
  //     TIterator* itr1 = components->createIterator();
  //     RooAbsArg* arg1;
  //     while ((arg1 = (RooAbsArg*)itr1->Next())) {
  // 	TIterator* itr2 = components->createIterator();
  // 	RooAbsArg* arg2;
  // 	while ((arg2 = (RooAbsArg*)itr2->Next())) {
  // 	  if (arg1 == arg2) 
  // 	    continue;
  // 	  if (arg2->dependsOn(*arg1)) {
  // 	    components->remove(*arg1);
  // 	  }
  // 	}
  // 	delete itr2;
  //     }
  //     delete itr1;
  //   }//end components getsize
  //   if (components->getSize() > 1) {
  //     cout << "ERROR::Couldn't isolate proper nuisance parameter" << endl;
  //     return;
  //   }
  //   else if (components->getSize() == 1) {
  //     thisNui = (RooRealVar*)components->first();
  //   }
    
  //   //Get the global observable corresponding on which the constraint depend.    
  //   TIterator* gIter = mc_globs.createIterator();
  //   RooRealVar* thisGlob = NULL;
  //   RooAbsArg* glob_arg;
  //   while ((glob_arg = (RooAbsArg*)gIter->Next())) {
  //     if (pdf->dependsOn(*glob_arg)) {
  // 	thisGlob = (RooRealVar*)glob_arg;
  // 	break;
  //     }
  //   }
  //   delete gIter;
    
  //   if (!thisNui || !thisGlob) {
  //     cout << "WARNING::Couldn't find nui or glob for constraint: " << pdf->GetName() << endl;
  //     continue;
  //   }
    
  //   cout << "Pairing nui: " << thisNui->GetName() << ", with glob: " << thisGlob->GetName() << ", from constraint: " << pdf->GetName() << endl;
    
  //   nui_list.add(*thisNui);
  //   glob_list.add(*thisGlob);  
  // } // end loop on cIter (constraints)
  // delete cIter;
  
  
  // combWs->saveSnapshot("nominalGlobs",glob_list);
  // combWs->saveSnapshot("nominalNuis", nui_list);

  
  // ///  RooArgSet nuiSet_tmp_before(nui_list);
  // //mu->setRange(-100,100.);
  // mu->setVal(0.0);
  // // combPdf->fitTo(*combData, Hesse(false), Minos(false), PrintLevel(0), Constrain(nuiSet_tmp_before));
  // // double mu_hat = mu->getVal(); 
  // // mu->setVal(mu_hat);
  // mu->setConstant(1);
  
  // //Set the spurious signal to 0
  // int nrNuis = nui_list.getSize();
  // cout << " Number of NP = "<< nrNuis << endl; 
  // for (int i=0;i<nrNuis;i++) {
  //   RooRealVar* nui = (RooRealVar*)nui_list.at(i);
  //   TString nui_name = nui->GetName();
  //   if ( !nui_name.Contains("spurious_") && !nui_name.Contains("BIAS")) continue;
  //   cout << " Fixing to constant nuisance parameters" << endl; 
  //   nui->setVal(0);
  //   nui->setConstant(1);
  // }
  
  // //fit data with backgroud only model.
  // RooArgSet nuiSet_tmp(nui_list);
  // RooFitResult *result2 = combPdf->fitTo(*combData,Hesse(false),Minos(false),PrintLevel(0),Extended(), Constrain(nuiSet_tmp), Save(),SumW2Error(kFALSE));
  // //  combPdf->fitTo(*combData, Hesse(false),Minos(false),PrintLevel(0),Extended(), Constrain(nuiSet_tmp));
  // cout << "Done" << endl;
  
  // mu->setConstant(0);
  
  // if (nrNuis != glob_list.getSize())  {
  //   cout << "ERROR::nui_list.getSize() != glob_list.getSize()!" << endl;
  //   return;
  // }
  
  // for (int i=0;i<nrNuis;i++) {
  //   RooRealVar* nui = (RooRealVar*)nui_list.at(i);
  //   RooRealVar* glob = (RooRealVar*)glob_list.at(i);
    
  //   cout << "nui: " << nui << ", glob: " << glob << endl;
  //   cout << "Setting glob: " << glob->GetName() << ", which had previous val: " << glob->getVal() << ", to conditional val: " << nui->getVal() << endl;
  //   glob->setVal(nui->getVal());  
  // }//end for i
  
  // //save the snapshots of conditional parameters
  // cout << "Saving conditional snapshots" << endl;
  // combWs->saveSnapshot(("conditionalGlobs"+muStr.str()).c_str(),glob_list);
  // combWs->saveSnapshot(("conditionalNuis" +muStr.str()).c_str(), nui_list);
  
  // if (!doConditional) {
  //   combWs->loadSnapshot("nominalGlobs");
  //   combWs->loadSnapshot("nominalNuis");
  // }
  
  // for (int i=0;i<nrNuis;i++) {
  //   RooRealVar* nui = (RooRealVar*)nui_list.at(i);
  //   TString nui_name2 = nui->GetName();
  //   if (!nui_name2.Contains("spurious_") && !nui_name2.Contains("BIAS")) continue;
  //   nui->setConstant(0);
  // }
  
  
  // cout << "Making asimov" << endl;
  // //  mu->setVal(!b_only);
  // // mH->setConstant();
  // // mH->setVal(125.5);

  // if (!b_only) {
  //   mu->setVal(1.);
  //   //    mu->setVal(1.);
  //   // mu_ggF->setVal(1.62);
  //   // mu_VBF->setVal(1.55);
  //   // mu_WH->setVal(1.47);
  //   // mu_ZH->setVal(2.92);
  // }
  // else 
  //   mu->setVal(0);
  // ModelConfig* mc = mcInWs;
  
  // int iFrame=0;
  
  // const char* weightName="weightVar";
  // RooArgSet obsAndWeight;
  // cout << "adding obs" << endl;
  // obsAndWeight.add(*mc->GetObservables());
  // cout << "adding weight" << endl;
  
  // RooRealVar* weightVar = NULL;
  // if (!(weightVar = combWs->var(weightName))) {
  //   combWs->import(*(new RooRealVar(weightName, weightName, 1,0,10000000)));
  //   //   combWs->import(*(new RooRealVar(weightName, weightName, 1,0,1000)));
  //   weightVar = combWs->var(weightName);    
  // }
  // //  cout << "weightVar: " << weightVar << endl;
  // obsAndWeight.add(*combWs->var(weightName));
  
  // //  cout << "defining set" << endl;
  // combWs->defineSet("obsAndWeight",obsAndWeight);

  
  // //////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////
  // // MAKE ASIMOV DATA FOR OBSERVABLES

  // // dummy var can just have one bin since it's a dummy
  // if(combWs->var("dummyX"))  
  //   combWs->var("dummyX")->setBins(1);

  // cout << endl << "Check expectedData by category" << endl;
  // //RooDataSet* simData=NULL;
  // RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
  // map<string, RooDataSet*> asimovDataMap;
  // //try fix for sim pdf
  // // RooCategory* channelCat = (RooCategory*)combWs->cat("master_channel");//(RooCategory*) (&simPdf->indexCat());
  // RooCategory* channelCat = (RooCategory*) (&simPdf->indexCat());

  // // TIterator* iter = simPdf->indexCat().typeIterator() ;
  // // TIterator* iter = channelCat->typeIterator() ;
  // // RooCatType* tt = NULL;
  // int nrIndices = channelCat->numTypes();

  // for (int i=0;i<nrIndices;i++) {
  //   channelCat->setIndex(i);
  //   iFrame++;
  //   // Get pdf associated with state from simpdf
  //   RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;
	
  //   // Generate observables defined by the pdf associated with this state
  //   RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;
  //   obstmp->Print();
    
  //   //channelCat->setIndex(string(tt->GetName()).c_str());
  //   cout << "on type " << channelCat->getLabel() << " " << iFrame << endl;
    
  //   RooDataSet* obsDataUnbinned = new RooDataSet(Form("combAsimovData%d",iFrame),Form("combAsimovData%d",iFrame),RooArgSet(obsAndWeight,*channelCat),WeightVar(*weightVar));
  //   RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
  //   double expectedEvents = pdftmp->expectedEvents(*obstmp);
  //   double thisNorm = 0;
  //   for(int jj=0; jj<thisObs->numBins(); ++jj) {
  //     thisObs->setBin(jj);
  //     //	cout << "pdf = "<<pdftmp->getVal(obstmp) <<endl;
  //     thisNorm = pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
  //     if (thisNorm*expectedEvents > pow(10., -2) && thisNorm*expectedEvents < pow(10., 9)) 
  // 	obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
  //   }
    
  //   obsDataUnbinned->Print();
  //   cout << "sum entries " << obsDataUnbinned->sumEntries() << endl;
  //   if(obsDataUnbinned->sumEntries() != obsDataUnbinned->sumEntries()){
  //     cout << "sum entries is nan"<<endl;
  //     exit(1);
  //   }
    
  //   ((RooRealVar*)obstmp->first())->Print();
  //   cout << "expected events " << pdftmp->expectedEvents(*obstmp) << endl;
    
  //   asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;//tempData;
    
  //   cout << "channel: " << channelCat->getLabel() << ", data: ";
  //   obsDataUnbinned->Print();
  //   cout << endl;
  // }
  
  
  // RooDataSet* asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));
  // if (combWs->data(asimovData->GetName())) {
  //   combWs->import(*asimovData, true); // Bool_t import(TObject& object, const char* aliasName, Bool_t replaceExisting = kFALSE)
  // }
  // else {
  //   combWs->import(*asimovData);
  // }
  
  
  // //  RooDataSet* asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));
  // //combWs->import(*asimovData);
  
  
  // //bring us back to nominal for exporting
  // combWs->loadSnapshot("nominalNuis");
  // combWs->loadSnapshot("nominalGlobs");
  
  // asimovData->Print();

  // delete result2;  


  // // TCanvas *c = new TCanvas();
  // // RooPlot* frame=combWs->var("invMass")->frame(100);
  // // frame->SetTitle(""); //empty title to prevent printing "A RooPlot of ..."
  // // frame->SetXTitle("m_{H}");
  // // char buffer_dummy[200]; 
  // // sprintf(buffer_dummy,"Events / %4.1f GeV",(frame->GetXaxis()->GetXmax()-frame->GetXaxis()->GetXmin())/frame->GetXaxis()->GetNbins());
  // // frame->SetYTitle(buffer_dummy);
  // // frame->SetXTitle("m_{#gamma#gamma} [GeV]");
  
  // // asimovData->plotOn(frame, MarkerColor(kRed));
  // // simPdf->plotOn(frame, LineColor(kBlue));
  
  // // frame->Draw();      
  // // c->SaveAs("test_asimov.png");

//}
