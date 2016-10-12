import argparse
import subprocess as sub
import os
import sys
from xml.dom import minidom
import xml.etree.cElementTree as ET
sys.path.append(os.path.abspath("/afs/in2p3.fr/home/c/cgoudet/private/Couplings/PlotFunctions/python"))
from SideFunction import *

categoriesNames = ['ggH_CenLow', 'ggH_CenHigh', "ggH_FwdLow", 'ggH_FwdHigh', "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep",   "ttHhad", "ttHlep"]

coefNames=['a', 'b', 'c', 'd' ]
paramNames=['mean', 'sigma' ]
formNames=[ 'CB' ]
processes=['ggH', 'VBF', 'VH', 'ttH' ]
subProcesses = ['bbH', 'tHjb', 'tWH' ]
#processes=['ggH', 'VBF', 'WH', 'ZH', 'ttH' ]

inputsFile='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/h012/StatisticsChallenge/h013/inputs/'
#inputsFile='/sps/atlas/e/escalier/HGamma/ProductionModes/FromMxAOD/h013/'

#====================================
def CategoryNode( catName, mode = 0 ) : 
    catIndex = categoriesNames.index( catName ) 
    xmlObj = CreateNode( 'category', { 'Name':catName, 'systFileName':inputsFile+'datacard_ICHEP.txt' } )

    for vProc in processes+subProcesses : xmlObj.append( CreateNode( 'yield', { 'process':vProc, 'inFileName':inputsFile+'workspace_signal_yields_categories.root', 'inVarName':'Yield_Signal_'+vProc+'_SM_'+catName } ) )
    xmlObj.append( CreateNode( 'pdf', {'process':'all', 'inFileName':inputsFile+'ModelSignal/RAW/SigSimple_all_shape_categories_DBCB/Individual/SM/res_SM_DoubleCB_workspace.root', 'inVarName':'sigPdf_SM_m125000_c'+str( catIndex ), 'invMass' : 'm_yy_m125000_c'+str(catIndex) } ) )

#Variables which have to be renamed
    varChanges = {}
    if mode == 0 :  
        varChanges = [
            { 'outName':'meanCB', 'inName': 'muCB_SM_m125000_c'+str(catIndex), 'systNP':'mean' },
            { 'outName':'sigmaCB', 'inName':'sigmaCB_SM_m125000_c'+str(catIndex), 'systNP':'sigma' },
#            { 'outName':'invMass', 'inName':'m_yy_m125000_c'+str(catIndex)} ,
            { 'outName':'mHcomb', 'inName':'mResonance', 'outVal':'125.09' }
            ]
        varChanges +=  [ { 'inName' : 'Yield_Signal_'+vProc+'_SM_'+catName, 'outName': 'yieldSignal', 'systNP':'yield' } for vProc in processes+subProcesses ]

        for vVarName in varChanges :  xmlObj.append( CreateNode( 'changeVar', vVarName ) )

#Marc's asimove have been generated simulating 10fb of data but with luminosity of 13fb. Need to rescale luminosity to 10
        for year in [ '2015', '2016' ] : xmlObj.append( CreateNode( 'changeVar', { 'inName':'lumi_'+year, 'scale':str(10/13.27676) } ) )

        dataNode = CreateNode( 'data' )

        # for year in ['15','16'] : dataNode.append( CreateNode('dataFile', 
        #                               { 'inFileName':'/sps/atlas/e/escalier/HGamma/ProductionModes/FromMxAOD/h013/data'+year+'/hist-data.root', 
        #                                 'varName':'m_yy', 
        #                                 'treeName':'tree_selected',
        #                                 'selectionCut':'catCoup_dev=='+str(catIndex+1),
        #                                 'selectionVars':'catCoup_dev'
        #                                 } ) )
        dataNode.append( CreateNode('dataFile', { 'inFileName':inputsFile+'PseudoData/ws_challenge_pseudo_data_'+catName+'.root', 'varName':'m_yy_'+catName, 'weightName':'weight', 'datasetName':'absdata_data_'+catName} ) )
        xmlObj.append( dataNode )

        correlatedVarNode = CreateNode( 'correlatedVar' )
        correlatedVarNode.text = 'mHcomb,lumi_2015,lumi_2016,' + ','.join( [ 'XS13_' + vProc + '_yy' for vProc in processes+subProcesses ] )
        xmlObj.append( correlatedVarNode )

        # createMass = CreateNode( 'createVar', {"Name" : "dependentMean", "formula": 'muCB_SM_m125000_c'+str(catIndex) + '+(mHcomb-125)'} )
        # createMass.text = 'mHcomb muCB_SM_m125000_c'+str(catIndex)

        bkg = CreateNode( "bkg", {'form':'expPol2' if 'ggH' in catName else "exp" } );
        xmlObj.append( bkg )

    return xmlObj
#================================================================
def ConfigFile( inFileName ) :
    if  inFileName == '' : print( 'No input name for config file.' ); exit(1);

    coreName = StripString(inFileName)
    xmlObj = CreateNode( 'CreateWorkspace', { 'Name':'/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/' + coreName + '.root' } )    

    processNode = CreateNode( 'processes' )
    processNode.text = ' '.join( processes+subProcesses )
    xmlObj.append( processNode )

    for vCatName in categoriesNames : xmlObj.append( CategoryNode(vCatName ) )
        

    configFile = open( StripString(inFileName, 0, 1 ) + '.xml', 'w+' )
    docTypeLine =  '<!DOCTYPE CreateWorkspace  SYSTEM "/afs/in2p3.fr/home/c/cgoudet/private/Couplings/Workspace/config/CreateWorkspace.dtd">'
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
        configFile.write( 'process=' + ' '.join( processes + subProcesses) + '\n' )

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
                                           for proc in processes+subProcesses  ] ) )
        
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
def SystEffectNode( varName='yield', upVal=0, downVal=0, constraint='Gaus') :
    if downVal==0: downVal = upVal
    if constraint == 'Gaus' and upVal!=downVal : print( 'Gaussian constraint but asym values' ); exit(0)

    xmlObj = ET.Element('systEffect')
    xmlObj.set( 'varName', varName )
    xmlObj.set( 'downVal', str(downVal) )
    xmlObj.set( 'upVal', str(upVal) )
    xmlObj.set( 'category', 'Common' )
    xmlObj.set( 'process', 'all' )


    return xmlObj

#=================================================
def AppendProc( node ) :
    for procName in ['all'] + processes + subProcesses : 
        node.append( CreateNode( 'process', { 'Name' : procName } ) )
                     
    return node
#==============================================
def FindProcInName( name ) : 
    for vProc in processes+subProcesses : 
        if vProc in name : return vProc
    return 'all'
#==============================================
def ParseLine( line, mapSyst, currentCat ) :
    line = line.replace( ' ', '' ).replace( '\n', '' )
 #   print('line : ' + line)
    name = line.split("=")[0]
 #   print( 'name : ' + name )

    systName = name.replace( currentCat, '' )
    if systName not in name : print( 'category name has troubles : ' + name ); return

#Check if the systematic 
    process = 'all'
    if 'ggH_ptH' in systName or 'ggH2in' in systName or 'ggH_m23' in systName : process = 'ggH'
    else :
        if 'ggH2j' in systName : process = 'VBF'
        else : process = FindProcInName( systName )
        systName = systName.replace( process, '' )

    systName = systName.replace( '__', '_' )
    while systName[-1]=='_' : systName = systName[:-1]
    if systName not in name : print( 'category name has troubles for process : ' + systName + ' ' + name  ); return

#Clean the name
    
#    print( 'systName : ' + systName )    
    dictOptionsSystEffect = { 'category':currentCat, 'process':process }    
    systValue = line.split('=')[1]
    dictOptionsSystEffect['constraint'] = 'Asym' if '-100L' in systValue else 'Gauss'

    if dictOptionsSystEffect['constraint'] != "Asym" : dictOptionsSystEffect['upVal'] = systValue
    else :
        systValue = systValue.replace( '-100L(', '' ).replace( ')', '' )
        minusSplitted = systValue.split('-')
        if systValue[0]=='-' : dictOptionsSystEffect['downVal'] = '-'+minusSplitted[1]
        else : dictOptionsSystEffect['downVal'] = systValue.minusSplitted[0]

        dictOptionsSystEffect['upVal'] = systValue.replace( dictOptionsSystEffect['downVal'], '' )[1:]

    systEffectNode = CreateNode( 'systEffect', dictOptionsSystEffect )

    if systName not in mapSyst : mapSyst[systName] = CreateNode( 'systematic', { 'Name':systName, 'centralValue': '0' if 'BIAS' in systName else '1' } )
    mapSyst[systName].append( systEffectNode )

    return 

#=================================================
def CreateXMLSystFromDataCard( inFileName, outFileName ) :
    print( 'CreateXMLSyst' )

    xmlObj = ET.Element("NPCorrelation")

    lumi = CreateNode( 'systematic', { 'correlation':'All', 'centralValue':'1', 'constraint':'Gaus', 'Name':'ATLAS_LUMI' } )
    lumi.append( CreateNode( 'systEffect', {'varName':'yield', 'upVal':'0.05', 'category':'Common', 'process':'all', 'constraint':'Gaus'} ) )
    ET.dump( lumi )
    xmlObj.append( lumi )

    dictSystNodes = {}

    inFile = open( inFileName )
    currentCat = ''
    for line in inFile : 
        if line[0]=='#' : continue 
        if '[' in line :
            currentCat = line[1:line.rfind(']')]
            continue
        if '=' not in line or 'bkg_model' in line : continue
        ParseLine( line, dictSystNodes, currentCat )
        

    for syst in dictSystNodes : xmlObj.append( dictSystNodes[syst] )

    outFile=open( outFileName, 'w' )
    outFile.write( '<!DOCTYPE NPCorrelation  SYSTEM "/afs/in2p3.fr/home/c/cgoudet/private/Couplings/Workspace/python/xmlCard.dtd">\n' )
    outFile.write( prettify( xmlObj ) )
    outFile.close()


#==========================================
def parseArgs():
    parser = argparse.ArgumentParser(
        description="This text will be displayed when the -h or --help option are given"
                    "It should contain a short description of what this program is doing.")

    parser.add_argument(
        '--doMode', help=( 'Tag for workspace mode\n' +
                           ' 0 : Marc asimov \n '+
                           ' 1 : ICHEP2016 data \n '
                           ),
        default=0, type=int )
    parser.add_argument('--outFileName', type=str, default='testConfig', help="Directory where all inputs are stored" )
    args = parser.parse_args()

    return args

#========================================
def main():
    # Parsing the command line arguments
    args = parseArgs()
    print( args.outFileName )
    args.outFileName = StripString( args.outFileName, 0, 1 )
    print( 'stripped' )
#    CreateXMLSystFromDataCard( '/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/h012/StatisticsChallenge/h013/inputs/datacard_ICHEP.txt', args.outFileName + '.xml' )    
    print( ConfigFile( args.outFileName + '.boost' ) )


# The program entrance
if __name__ == '__main__':
    main()
