categoriesNames = ["ggH", "VBF_low", "VBF_high", "VH_hadronic_low", "VH_hadronic_high", "VH_MET", "VH_leptonic", "VH_dileptons", "ttH_hadronic", "ttH_leptonic"]
#categoriesNames = ["ggH"]
coefNames=['a', 'b', 'c', 'd' ]
paramNames=['mean', 'sigma' ]
formNames=['GA', 'CB' ]
processes=['ggH', 'VBF', 'WH', 'ZH', 'ttH' ]

workspaceName='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/RAW_SignalModel/SigParam_%s_categories/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root'
workspaceYield='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/RAW_SignalModel/SigParam_all_shape_categories/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root'
datacard='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011_fix/datacard.txt'
model=3;

if model==0 :
    with open( "StatChallenge011.boost", 'w+' ) as configFile :
        configFile.write( "[General]\n" )
        configFile.write( 'catNames=' + ' '.join( categoriesNames ) + '\n' )
        configFile.write( 'process=' + ' '.join( processes ) + ' bbH \n' )
        configFile.write( 'outName=/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h011_fullreco.root\n' )
        configFile.write( 'systFileName=' + datacard + '\n' )
        
        
        for iCat in range( 0, len( categoriesNames ) ) : 
            configFile.write( '\n' )
            configFile.write( '[' + categoriesNames[iCat] + ']\n' )

            configFile.write( 'signalModel=1\n' )
            configFile.write( '\n'.join( [ coef + '_yieldVar_' + proc
                                           + '=' + workspaceName.replace('%s', proc) + ' yieldVar_' + coef + '_SM_c'+ str( iCat )
                                           for coef in coefNames for proc in processes ] ) )
            for proc in processes :        
                configFile.write( '\n'.join( [ coef + '_' + param + form + '_' + proc 
                                           +'=' + workspaceName.replace('%s', proc) + ' ' + coef + '_' + ( param if param != "mean" else "mu" ) + form + 'Nom_SM_c' + str( iCat )
                                               for coef in coefNames for param in paramNames for form in formNames  ] ) )
                configFile.write( '\n' )
                configFile.write( '\n'.join( [ coef + '_alphaCB_' + proc
                                               + '=' + workspaceName.replace('%s', proc ) + ' ' + coef + '_alphaCB_SM_c'+ str( iCat )
                                               for coef in coefNames  ] ) )
                configFile.write( '\n' )
                configFile.write( '\n'.join( [ 'nCB_'+ proc
                                               + '=' + workspaceName.replace( "%s", proc ) + ' nCB_SM_c'+ str( iCat )
                                               ] ) )
                configFile.write( '\n' )
                configFile.write( 'fCB_'+ proc + '=' + workspaceName.replace( "%s", proc ) + ' fracCB_SM_c'+ str( iCat ) )
                configFile.write( '\n' )
                
            configFile.write( '\n' )
            configFile.write( 'dataFileName = /sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/ws_challenge_pseudo_data.root ws_challenge_pseudo_data absdata_data m_yy category\n' )
            configFile.write( 'dataCut=category==category::Channel_' + categoriesNames[iCat] + '\n' )
                

elif model == 1 :
    with open( "StatChallenge011_pdf.boost", 'w+' ) as configFile :
        configFile.write( "[General]\n" )
        configFile.write( 'catNames=' + ' '.join( categoriesNames ) + '\n' )
        configFile.write( 'process=' + ' '.join( processes ) + ' bbH \n' )
        configFile.write( 'systFileName=' + datacard + '\n' )
        configFile.write( 'outName=/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h011_pdfReco.root\n' )    

        for iCat in range( 0, len( categoriesNames ) ) : 
            configFile.write( '\n' )
            configFile.write( '[' + categoriesNames[iCat] + ']\n' )

            configFile.write( 'signalModel=0\n' )
                   
            configFile.write( '\n'.join( [ "signal_" + proc 
                                           +'=' + workspaceYield + ' sigPdf_SM_c' + str( iCat)
                                           for proc in processes  ] ) )
            configFile.write( '\n' )
            configFile.write( '\n'.join( [ "yield_" + proc 
                                           +'=' + workspaceName.replace('%s', proc) + ' sigYield_SM_c' + str( iCat)
                                           for proc in processes  ] ) )
            configFile.write( '\n' )

            configFile.write( '\n'.join( [ param + form + '_' + proc 
                                           +'=' + (param if param != "mean" else 'mu')  + form + '_SM_c' + str( iCat )
                                           for param in paramNames for form in formNames for proc in processes  ] ) )

            configFile.write( '\n' )
            configFile.write( '\n'.join( [ 'nCB_'+ proc
                                           + '= nCB_SM_c'+ str( iCat )
                                           for proc in processes] ) )
            configFile.write( '\n' )
            configFile.write( '\n'.join( [ 'alphaCB_'+ proc
                                           + '= alphaCB_SM_c'+ str( iCat )
                                           for proc in processes] ) )
            configFile.write( '\n' )
            configFile.write( '\n'.join( ['fCB_'+ proc
                              + '= fracCB_SM_c'+ str( iCat ) 
                                          for proc in processes] ) )
            configFile.write( '\n' )
            configFile.write( 'invMass=m_yy' )
            configFile.write( '\n' )
            configFile.write( 'mHcomb=mResonance' )
            configFile.write( '\n' )

            configFile.write( '\n' )
            configFile.write( 'dataFileName = /sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/ws_challenge_pseudo_data.root ws_challenge_pseudo_data absdata_data m_yy category\n' )
            configFile.write( 'dataCut=category==category::Channel_' + categoriesNames[iCat] + '\n' )

elif model == 3 :
    workspaceName='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge/StatisticsChallenge/h011/inputs/ModelSignal/RAW/SigParam_%s_categories/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root'
    workspaceYield='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge/StatisticsChallenge/h011/inputs/ModelSignal/RAW/SigParam_all_shape_categories/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root'

    newProcesses = [ 'bbH', 'tHjb', 'tWH' ]
    with open( "StatChallenge011_asimov.boost", 'w+' ) as configFile :
        configFile.write( "[General]\n" )
        configFile.write( 'catNames=' + ' '.join( categoriesNames ) + '\n' )
        configFile.write( 'process=' + ' '.join( processes + newProcesses) + '\n' )
        configFile.write( 'systFileName=/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge/StatisticsChallenge/h011/inputs/datacard.txt\n' )
        configFile.write( 'outName=/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_asimov.root\n' )    


        for iCat in range( 0, len( categoriesNames ) ) : 
            configFile.write( '\n' )
            configFile.write( '[' + categoriesNames[iCat] + ']\n' )

            configFile.write( 'signalModel=0\n' )
                   
            configFile.write( '\n'.join( [ "signal_" + proc 
                                           +'=' + workspaceName.replace('%s', 'all_shape') + ' sigPdf_SM_c' + str( iCat)
                                           for proc in processes+newProcesses  ] ) )
            configFile.write( '\n' )
            configFile.write( '\n'.join( [ "yield_" + proc 
                                           +'=' + workspaceName.replace('%s', proc) + ' sigYield_SM_c' + str( iCat)
                                           for proc in processes  ] ) )
            configFile.write( '\n' )
            configFile.write( '\n'.join( [ "yield_" + proc 
                                           +'=/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge/StatisticsChallenge/h011/inputs/workspace_Yield_Signal_bbH_tHjb_tWH.root  Yield_Signal_tWH_SM_Per1fbMinus1_Channel_' + categoriesNames[iCat]
                                           for proc in newProcesses  ] ) )

            configFile.write( '\n' )

            configFile.write( '\n'.join( [ param + form + '_' + proc 
                                           +'=' + (param if param != "mean" else 'mu')  + form + '_SM_c' + str( iCat )
                                           for param in paramNames for form in formNames for proc in processes+newProcesses  ] ) )

            configFile.write( '\n' )
            configFile.write( '\n'.join( [ 'nCB_'+ proc
                                           + '= nCB_SM_c'+ str( iCat )
                                           for proc in processes+newProcesses] ) )
            configFile.write( '\n' )
            configFile.write( '\n'.join( [ 'alphaCB_'+ proc
                                           + '= alphaCB_SM_c'+ str( iCat )
                                           for proc in processes+newProcesses] ) )
            configFile.write( '\n' )
            configFile.write( '\n'.join( ['fCB_'+ proc
                              + '= fracCB_SM_c'+ str( iCat ) 
                                          for proc in processes+newProcesses] ) )
            configFile.write( '\n' )
            configFile.write( 'invMass=m_yy' )
            configFile.write( '\n' )
            configFile.write( 'mHcomb=mResonance' )
            configFile.write( '\n' )

            configFile.write( '\n' )
            configFile.write( 'dataFileName = /sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge/StatisticsChallenge/h011/inputs/ws_challenge_pseudo_data_Channel_' + categoriesNames[iCat] + '.root ws_challenge_pseudo_data_Channel_' + categoriesNames[iCat] + ' absdata_data_Channel_' + categoriesNames[iCat] + ' m_yy\n' )
            configFile.write( 'dataWeight=weight\n' )
