import argparse
import subprocess as sub
import os
import sys
from xml.dom import minidom
import xml.etree.cElementTree as ET
sys.path.append(os.path.abspath("/afs/in2p3.fr/home/c/cgoudet/private/Calibration/PlotFunctions/python"))
from SideFunction import *

categoriesNames = ['ggH_CenLow', 'ggH_CenHigh', "ggH_FwdLow", 'ggH_FwdHigh', "VBFloose", "VBFtight", "VHMET", "VHlep", "VHdilep", "VHhad_loose", "VHhad_tight",  "ttHhad", "ttHlep"]

coefNames=['a', 'b', 'c', 'd' ]
paramNames=['mean', 'sigma' ]
formNames=[ 'CB' ]
processes=['ggH', 'VBF', 'VH', 'ttH' ]
subProcesses = ['bbH', 'tHjb', 'tWH' ]
#processes=['ggH', 'VBF', 'WH', 'ZH', 'ttH' ]

inputsFile='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/h012/StatisticsChallenge/h013/inputs/'
#inputsFile='/sps/atlas/e/escalier/HGamma/ProductionModes/FromMxAOD/h013/'
newProcesses = [ 'bbH', 'tHjb', 'tWH' ]

#====================================
def CategoryNode( catName, mode = 0 ) : 
    catIndex = categoriesNames.index( catName ) 
    xmlObj = CreateNode( 'category', { 'Name':catName, 'systFileName':inputsFile+'datacard_ICHEP.txt' } )

    for vProc in processes+subProcesses : xmlObj.append( CreateNode( 'yield', { 'process':vProc, 'inFileName':inputsFile+'workspace_signal_yields_categories.root', 'inVarName':'Yield_Signal_'+vProc+'_SM_'+catName } ) )
    xmlObj.append( CreateNode( 'pdf', {'process':'All', 'inFileName':inputsFile+'ModelSignal/RAW/SigSimple_all_shape_categories_DBCB/Individual/SM/res_SM_DoubleCB_workspace.root', 'inVarName':'sigPdf_SM_m125000_c'+str( catIndex ) } ) )

#Variables which have to be renamed
    varChanges = {}
    if mode == 0 :  
        varChanges = {
            'meanCB_all': 'muCB_SM_m125000_c'+str(catIndex),
            'sigmaCB_all':'sigmaCB_SM_m125000_c'+str(catIndex),
            'invMass':'m_yy_m125000_c'+str(catIndex),
            'mHcomb': 'mResonnance'
            }
        for vVarName in varChanges :  xmlObj.append( CreateNode( 'changeVar', { 'outName':vVarName, 'inName':varChanges[vVarName] } ) )

#Marc's asimove have been generated simulating 10fb of data but with luminosity of 13fb. Need to rescale luminosity to 10
        for year in [ '20015', '2016' ] : xmlObj.append( CreateNode( 'changeVar', { 'inName':'lumi_'+year, 'scale':str(10/13.27676) } ) )

        dataNode = CreateNode( 'data' )
        dataNode.append( CreateNode('dataFile', { 'inFileName':inputsFile+'PseudoData/ws_challenge_pseudo_data_'+catName+'.root', 'varName':'m_yy_'+catName, 'weightName':'weight', 'datasetName':'absdata_data_'+catName} ) )
        xmlObj.append( dataNode )
    return xmlObj
#================================================================
def ConfigFile( inFileName ) :
    if  inFileName == '' : print( 'No input name for config file.' ); exit(1);

    coreName = StripString(inFileName,0, 1)
    xmlObj = CreateNode( 'CreateWorkspace', { 'Name':'/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/' + coreName + '.root' } )    

    processNode = CreateNode( 'processes' )
    processNode.text = ' '.join( processes+subProcesses )
    xmlObj.append( processNode )

    for vCatName in categoriesNames : xmlObj.append( CategoryNode(vCatName ) )
        

    configFile = open( coreName + '.xml', 'w+' )
    docTypeLine =  '<!DOCTYPE CreateWorkspace  SYSTEM "/afs/in2p3.fr/home/c/cgoudet/private/Couplings/Workspace/python/CreateWorkspace.dtd">'
    stringToWrite = prettify( xmlObj ).split('\n')
    stringToWrite.insert( 1, docTypeLine )
    configFile.write( '\n'.join( stringToWrite ) )

    return inFileName

#================================================================
def ConfigFileBoost( inFileName ) :
    if inFileName == '' : 
        print( 'ConfigFile is missing name' )
        exit(0)

    with open( inFileName, 'w+' ) as configFile :
        configFile.write( "[General]\n" )
        configFile.write( 'catNames=' + ' '.join( categoriesNames ) + '\n' )
        configFile.write( 'process=' + ' '.join( processes + newProcesses) + '\n' )

 #
#        configFile.write( 'systFileName=' + inFileName.replace('.boost', '.xml' ) + '\n' )
        configFile.write( 'outName=/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/' + inFileName.replace('.boost', '.root' )+'\n' )    
#        configFile.write( 'yieldScale=' + str( 10/13.27676 ) )
        
        for iCat in range( 0, len( categoriesNames ) ) : 
            configFile.write( '\n' )
            configFile.write( '[' + categoriesNames[iCat] + ']\n' )
            configFile.write( 'systFileName=/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/ICHEP_2016/config/category_run2_'+ categoriesNames[iCat] + '.xml\n' )            
#           configFile.write( 'systFileName=' + inputsFile + 'datacard_ICHEP.txt\n' )
            #            configFile.write( 'signalModel=0\n' )
            
            configFile.write( '\n'.join( [ "signal_" + proc 
                                           +'=' + inputsFile + 'ModelSignal/RAW/SigSimple_all_shape_categories_DBCB/Individual/SM/res_SM_DoubleCB_workspace.root sigPdf_SM_m125000_c' + str( iCat)
                                           for proc in ['all']  ] ) )
            configFile.write( '\n' )
            configFile.write( '\n'.join( [ "yield_" + proc 
                                           +'=' + inputsFile + 'workspace_signal_yields_categories.root Yield_Signal_'+proc+'_SM_'+categoriesNames[iCat]
                                           for proc in processes+newProcesses  ] ) )
        
            configFile.write( '\n' )
            configFile.write( '\n'.join( [ param + form + '_' + proc 
                                           +'=' + (param if param != "mean" else 'mu')  + form + '_SM_m125000_c' + str( iCat )
                                           for param in paramNames for form in formNames for proc in ['all' ] ] ) )
            
            configFile.write( '\n' )
            configFile.write( 'invMass=m_yy_m125000_c' + str(iCat) )
            configFile.write( '\n' )
            configFile.write( 'mHcomb=mResonance' )
            configFile.write( '\n' )
        
            configFile.write( '\n' )
 #           configFile.write( 'dataFileName ='+ inputsFile + 'PseudoData/ws_challenge_pseudo_data_' + categoriesNames[iCat] + '.root ws_challenge_pseudo_data_' + categoriesNames[iCat] + ' absdata_data_' + categoriesNames[iCat] + ' m_yy_' + categoriesNames[iCat] +'\n' )
            configFile.write( 'dataFileName =/sps/atlas/e/escalier/HGamma/ProductionModes/FromMxAOD/h013/data15/data15/hist-data.root tree_selected m_yy\n' )
            configFile.write( 'dataFileName =/sps/atlas/e/escalier/HGamma/ProductionModes/FromMxAOD/h013/data15/data16/hist-data.root tree_selected m_yy\n' )
            configFile.write( 'dataWeight=weight\n' )

            configFile.write( 'changeVarName=' + ' '.join( ['lumi_2015', 'lumi_2016' ] ) +'\n')
            configFile.write( 'changeVarVal=' + ' '.join( ['2.42', '7.58'] ) +'\n')

    return inFileName


#=================================================
def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """

    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")
#================================================
def CreateNode( nodeName, options={} ) :
    xmlObj = ET.Element( nodeName ) 
    for opt in options : xmlObj.set( opt, options[opt] )
    return xmlObj

#=================================================
def SystEffectNode( varName='yield', upVal=0, downVal=0, constraint='Gaus') :
    if downVal==0: downVal = upVal
    if constraint == 'Gaus' and upVal!=downVal : print( 'Gaussian constraint but asym values' ); exit(0)

    xmlObj = ET.Element('systEffect')
    xmlObj.set( 'varName', varName )
    xmlObj.set( 'downVal', str(downVal) )
    xmlObj.set( 'upVal', str(upVal) )
    xmlObj.set( 'category', 'Common' )
    xmlObj.set( 'process', 'All' )


    return xmlObj

#=================================================
def AppendProc( node ) :
    for procName in ['All'] + processes + subProcesses : 
        node.append( CreateNode( 'process', { 'Name' : procName } ) )
                     
    return node
#=================================================
def CreateStruct(node) :
    for catName in ['Common', 'Inclusive' ]+ categoriesNames : 
        catNode = ET.Element('category')
        catNode.set( 'name', catName )
        AppendProc( catNode )
        node.append( catNode )

#=================================================
def CreateXMLSystFromDataCard( inFileName, outFileName ) :
    print( 'CreateXMLSyst' )

    xmlObj = ET.Element("NPCorrelation")

    lumi = CreateNode( 'systematic', { 'correlation':'All', 'centralValue':'1', 'constraint':'Gaus', 'Name':'ATLAS_LUMI' } )
    lumi.append( CreateNode( 'systEffect', {'varName':'yield', 'upVal':'0.05', 'category':'Common', 'process':'All', 'constraint':'Gaus'} ) )
    ET.dump( lumi )

    xmlObj.append( lumi )

    inFile = open( inFileName )
    currentCat = ''
    for line in inFile :
        if line[0]=='#' : continue 
        if '[' in line :
            line = line[1:line.rfind(']')]
            currentCat = line
            continue
        if '=' not in line : continue

        name = line.split("=")[0]


    outFile=open( outFileName, 'w' )
    outFile.write( '<!DOCTYPE NPCorrelation  SYSTEM "/afs/in2p3.fr/home/c/cgoudet/private/Couplings/Workspace/python/xmlCard.dtd">\n' )
    outFile.write( prettify( xmlObj ) )
    outFile.close()


#==========================================
def parseArgs():
    """
    ArgumentParser.
    Return
    ------
        args: the parsed arguments.
    """
    # First create a parser with a short description of the program.
    # The parser will automatically handle the usual stuff like the --help messages.
    parser = argparse.ArgumentParser(
        description="This text will be displayed when the -h or --help option are given"
                    "It should contain a short description of what this program is doing.")

    # Retrieve the optionnal arguments. Optional argument should start with '-' or '--'
    # Simply give to the parser.add_argument() function the expected names 
    # (one short '-' and/or one long '--')

    # The parser can cast the given arguments to some usual types and will return 
    # standard error message if the given argument is not of the good type.

    # Integers
    # Here I give the short and the long argument name
    parser.add_argument(
        '--doMode', help=( 'Tag for workspace mode\n' +
                           ' 0 : Marc asimov \n '+
                           ' 1 : ICHEP2016 data \n '
                           ),
        default=0, type=int )
    # parser.add_argument(
    #     '--doLatex', help='Tag for recreating plots',
    #     default=1, type=int )
    # parser.add_argument(
    #     '--doSyst', help='Tag for recreating systematics histos and plots',
    #     default=0, type=int )
    # parser.add_argument(
    #     '--doCorrection', help='Tag for recreating systematics histos and plots',
    #     default=0, type=int )

    # Floats
    # Here I give only the long argument name

    # Choices
    # Here is illustrated the possibility of having to choose between several values
    # Again help and error messages are automatically built.
    # Moreover when you do not precise the dest= argument it choose in order of priority
    # the long name and the short name (if the long one does not exist)
    # parser.add_argument("--verbosity", "-v", type=int, choices=[0, 1, 2],
    #                 help="increase output verbosity")

    # Add an argument to handle the output file.  If we use argparse.FileType it will
    # handle opening a writeable file (and ensuring we can write to it).
    # (Note: the str argument follow teh standard names 
    #   'w'=write, 'r'=read, 'r+'= read and write, 'a'=append,
    # For Windows there is a distinction between text and binary files: use the 'b' option
    #   'rb'=read binary, 'wb'=write binary, 'r+b'= read and write, etc)
    # Note that the argument name 'file' does not beging with a '-' or '--'; this indicates
    # to argparse that it is a *positional* argument
    # Whoa 2 new things in a single example ! This is crazy !
    # parser.add_argument('file', type=argparse.FileType('w'),
    #                     help='the name of the output file')

    parser.add_argument('--outFileName', type=str, default='testConfig', help="Directory where all inputs are stored" )
    # Now do your job and parse my command line !
    args = parser.parse_args()

    return args

#========================================
def main():
    """
    The main script
    """
    # Parsing the command line arguments
    args = parseArgs()

    inFileName = 'StatChallenge013'
#    CreateXMLSystFromDataCard( '/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/h012/StatisticsChallenge/h013/inputs/datacard_ICHEP.txt', args.outFileName + '.xml' )    
    print( ConfigFile( args.outFileName + '.boost' ) )


# The program entrance
if __name__ == '__main__':
    main()
