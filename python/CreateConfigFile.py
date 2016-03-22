categoriesNames = ["ggH", "VBF_low", "VBF_high", "VH_hadronic_low", "VH_hadronic_high", "VH_MET", "VH_leptonic", "VH_dileptons", "ttH_hadronic", "ttH_leptonic"]
coefNames=['a', 'b', 'c', 'd' ]
paramNames=['mean', 'sigma' ]
formNames=['GA', 'CB' ]
processes=['ggH', 'VBF', 'WH', 'ZH', 'ttH' ]

workspaceName='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/RAW_SignalModel/SigParam_%s_categories/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root'
workspaceYield='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/RAW_SignalModel/SigParam_all_shape_categories/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root'
with open( "StatChallenge011.boost", 'w+' ) as configFile :
    configFile.write( "[General]\n" )
    configFile.write( 'catNames=' + ' '.join( categoriesNames ) + '\n' )
    configFile.write( 'process=' + ' '.join( processes ) + '\n' )
    configFile.write( 'systFileName=/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/datacard.txt' )


    for iCat in range( 0, len( categoriesNames ) ) : 
        configFile.write( '\n' )
        configFile.write( '[' + categoriesNames[iCat] + ']\n' )
        configFile.write( '\n'.join( [ coef + '_' + param + form + '_' + proc 
                                       +'=' + workspaceName.replace('%s', proc) + ' ' + coef + '_' + ( param if param != "mean" else "mu" ) + form + 'Nom_SM_c' + str( iCat )
                                       for coef in coefNames for param in paramNames for form in formNames for proc in processes ] ) )
        configFile.write( '\n' )
        configFile.write( '\n'.join( [ coef + '_alphaCB_' + proc
                                       + '=' + workspaceName.replace('%s', proc ) + ' ' + coef + '_alphaCB_SM_c'+ str( iCat )
                                       for coef in coefNames for proc in processes ] ) )
        configFile.write( '\n' )
        configFile.write( '\n'.join( [ coef + '_yieldVar_' + proc
                                       + '=' + workspaceYield + ' ' + coef + '_yieldVar_SM_c'+ str( iCat )
                                       for coef in coefNames for proc in processes ] ) )
        configFile.write( '\n' )
        configFile.write( 'dataFileName = /sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/StatChallenge_h011/ws_challenge_pseudo_data.root ws_challenge_pseudo_data absdata_data m_yy \n' )
        configFile.write( 'dataCut=category==category::Channel_' + categoriesNames[iCat] + '\n' )
