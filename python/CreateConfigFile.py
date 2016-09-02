categoriesNames = ["ggH_FwdLow", 'ggH_FwdHigh','ggH_CenLow', 'ggH_CenHigh', "VBFloose", "VBFtight", "VHMET", "VHlep", "VHdilep", "VHhad_loose", "VHhad_tight",  "ttHhad", "ttHlep"]
#
#categoriesNames = ["ggH"]
coefNames=['a', 'b', 'c', 'd' ]
paramNames=['mean', 'sigma' ]
formNames=[ 'CB' ]
processes=['ggH', 'VBF', 'WH', 'ZH', 'ttH' ]

inputsFile='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/h012/StatisticsChallenge/h013/inputs/'
newProcesses = [ 'bbH', 'tHjb', 'tWH' ]

with open( "StatChallenge013.boost", 'w+' ) as configFile :
    configFile.write( "[General]\n" )
    configFile.write( 'catNames=' + ' '.join( categoriesNames ) + '\n' )
    configFile.write( 'process=' + ' '.join( processes + newProcesses) + '\n' )
    configFile.write( 'systFileName=' + inputsFile + 'datacard_ICHEP.txt\n' )
    configFile.write( 'outName=/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h013.root\n' )    
    configFile.write( 'yieldScale=' + str( 10/13.27676 ) )

    for iCat in range( 0, len( categoriesNames ) ) : 
        configFile.write( '\n' )
        configFile.write( '[' + categoriesNames[iCat] + ']\n' )
        
        configFile.write( 'signalModel=0\n' )
        
        configFile.write( '\n'.join( [ "signal_" + proc 
                                       +'=' + inputsFile + 'ModelSignal/RAW/SigSimple_all_shape_categories_DBCB/Individual/SM/res_SM_DoubleCB_workspace.root sigPdf_SM_m125000_c' + str( iCat)
                                       for proc in processes+newProcesses  ] ) )
        configFile.write( '\n' )
        configFile.write( '\n'.join( [ "yield_" + proc 
                                       +'=' + inputsFile + 'workspace_signal_yields_categories.root Yield_Signal_'+proc+'_SM_'+categoriesNames[iCat]
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
        configFile.write( 'dataFileName ='+ inputsFile + 'PseudoData/ws_challenge_pseudo_data_' + categoriesNames[iCat] + '.root ws_challenge_pseudo_data_' + categoriesNames[iCat] + ' absdata_data_' + categoriesNames[iCat] + ' m_yy_' + categoriesNames[iCat] +'\n' )
        configFile.write( 'dataWeight=weight\n' )

        configFile.write( 'changeVarName=' + ' '.join( ['lumi_2015', 'lumi_2016' ] ) +'\n')
        configFile.write( 'changeVarVal=' + ' '.join( ['2.42', '7.58'] ) +'\n')

