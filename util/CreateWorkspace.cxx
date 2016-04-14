#include "Workspace/Workspace.h"
#include <iostream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using std::cout;
using std::endl;

int main( int argc, char* argv[] ) {

    //Check validity of Arguments                                                      
    po::options_description desc("LikelihoodProfiel Usage");
    string configFile;
    //define all options in the program                                                
    desc.add_options()
      ("help", "Display this help message")
      ( "configFile", po::value<string>( &configFile ), "" )
      ;

    //Define options gathered by position                                              
    po::positional_options_description p;
    p.add("configFile", -1);

    // create a map vm that contains options and all arguments of options              
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
    po::notify(vm);
    if(vm.count("help")) {cout << desc; return 0;}


    Workspace ws;
    ws.Configure( configFile.c_str() );
    //    ws.Configure( "/afs/in2p3.fr/home/c/cgoudet/private/Couplings/Workspace/python/StatChallenge011.boost" );
    //ws.Configure( "/afs/in2p3.fr/home/c/cgoudet/private/Couplings/Workspace/python/StatChallenge011_pdf.boost" );
    ws.CreateWS();
    return 0;
}
