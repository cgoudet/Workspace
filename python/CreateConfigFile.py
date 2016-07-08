categoriesNames = ["ggH_FwdLow", 'ggH_FwdHigh','ggH_CenLow', 'ggH_CenHigh', "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight",  "ttHhad", "ttHlep"]
#, "VHMET", "VHlep", "VHdilep"
#categoriesNames = ["ggH"]
coefNames=['a', 'b', 'c', 'd' ]
paramNames=['mean', 'sigma' ]
formNames=[ 'CB' ]
processes=['ggH', 'VBF', 'WH', 'ZH', 'ttH' ]

inputsFile='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/h012/StatisticsChallenge/h012/inputs/'
newProcesses = [ 'bbH', 'tHjb', 'tWH' ]

with open( "StatChallenge012_asimov.boost", 'w+' ) as configFile :
    configFile.write( "[General]\n" )
    configFile.write( 'catNames=' + ' '.join( categoriesNames ) + '\n' )
    configFile.write( 'process=' + ' '.join( processes + newProcesses) + '\n' )
    configFile.write( 'systFileName=' + inputsFile + 'datacard_2015_plus_DS1.txt\n' )
    configFile.write( 'outName=/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h012_asimov.root\n' )    


    for iCat in range( 0, len( categoriesNames ) ) : 
        configFile.write( '\n' )
        configFile.write( '[' + categoriesNames[iCat] + ']\n' )
        
        configFile.write( 'signalModel=0\n' )
        
        configFile.write( '\n'.join( [ "signal_" + proc 
                                       +'=' + inputsFile + 'ModelSignal/RAW/SigSimple_all_shape_categories_DBCB/Individual/SM/res_SM_DoubleCB_workspace.root sigPdf_SM_m125000_c' + str( iCat)
                                       for proc in processes+newProcesses  ] ) )
        configFile.write( '\n' )
        configFile.write( '\n'.join( [ "yield_" + proc 
                                       +'=' + inputsFile + 'ModelSignal/workspace_signal_yields_categories.root Yield_Signal_'+proc+'_SM_'+categoriesNames[iCat]
                                       for proc in processes+newProcesses  ] ) )
        
        configFile.write( '\n' )
        configFile.write( '\n'.join( [ param + form + '_' + proc 
                                       +'=' + (param if param != "mean" else 'mu')  + form + '_SM_m125000_c' + str( iCat )
                                       for param in paramNames for form in formNames for proc in processes+newProcesses  ] ) )
        
        configFile.write( '\n' )
        configFile.write( 'invMass=m_yy_m125000_c' + str(iCat) )
        configFile.write( '\n' )
        configFile.write( 'mHcomb=mResonance' )
        configFile.write( '\n' )
        
        configFile.write( '\n' )
        configFile.write( 'dataFileName ='+ inputsFile + '/PseudoData/ws_challenge_pseudo_data_' + categoriesNames[iCat] + '.root ws_challenge_pseudo_data_' + categoriesNames[iCat] + ' absdata_data_' + categoriesNames[iCat] + ' m_yy_' + categoriesNames[iCat] +'\n' )
        configFile.write( 'dataWeight=weight\n' )
        
