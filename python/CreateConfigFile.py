categoriesNames = ["ggH", "VBF_low", "VBF_high", "VH_hadronic_low", "VH_hadronic_high", "VH_MET", "VH_leptonic", "VH_dileptons", "ttH_hadronic", "ttH_leptonic"]
coefNames=['a', 'b', 'c', 'd' ]
paramNames=['mean', 'sigma' ]
formNames=['GA', 'CB' ]
processes=['ggH', 'VBF', 'WH', 'ZH', 'ttH' ]

workspaceName='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/RAW_SignalModel/SigParam_%s_categories/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root'
workspaceYield='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/RAW_SignalModel/SigParam_all_shape_categories/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root'
model=1;

if model==0 :
    with open( "StatChallenge011.boost", 'w+' ) as configFile :
        configFile.write( "[General]\n" )
        configFile.write( 'catNames=' + ' '.join( categoriesNames ) + '\n' )
        configFile.write( 'process=' + ' '.join( processes ) + '\n' )
        configFile.write( 'outName=/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h011_fullreco.root\n' )
        configFile.write( 'systFileName=/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/datacard.txt\n' )
        
        
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
        configFile.write( 'process=' + ' '.join( processes ) + '\n' )
        configFile.write( 'systFileName=/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/datacard.txt\n' )
    

        for iCat in range( 0, len( categoriesNames ) ) : 
            configFile.write( '\n' )
            configFile.write( '[' + categoriesNames[iCat] + ']\n' )

            configFile.write( 'signalModel=1\n' )
            # configFile.write( '\n'.join( [ coef + '_yieldVar_' + proc
            #                                + '=' + workspaceName.replace('%s', proc) + ' yieldVar_' + coef + '_SM_c'+ str( iCat )
            #                                for coef in coefNames for proc in processes ] ) )
            # configFile.write( '\n' )
                   
            configFile.write( '\n'.join( [ "signal_" + proc 
                                           +'=' + workspaceName.replace('%s', proc) + ' sigPdf_SM_c' + str( iCat)
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
