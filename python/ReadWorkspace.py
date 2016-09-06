import os

mode=1


inFile=''
functions=[]
processes=['ggH', 'VBF', 'VH', 'ttH', 'bbH', 'tHjb', 'tWH' ]
categoriesNames = ["ggH_FwdLow", 'ggH_FwdHigh','ggH_CenLow', 'ggH_CenHigh', "VBFloose", "VBFtight", "VHMET", "VHlep", "VHdilep", "VHhad_loose", "VHhad_tight",  "ttHhad", "ttHlep"]
CBPars = ['nCBLo_SM_c', 'nCBHi_SM_c', 'alphaCBHi_SM_m125000_c', 'alphaCBLo_SM_m125000_c', 'muCBNom_SM_m125000_c', 'muCBNom_SM_m125000_c' ]

outName=''
if ( mode == 0 ) :
    inFile='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/h012/StatisticsChallenge/h013/WorkspaceForStudies/mu_i/file_workspace_CATEG_ICHEP_2015_2016.root'
    functions = [ 'Yield_Signal_' + proc + '_SM_'+cat for proc in processes for cat in categoriesNames ]
    functions += [ name + str(cat) for name in CBPars for cat in range(0, 13)]
    outName='marc'

elif ( mode == 1 ) :
    inFile='/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h013.root'
    functions = [ 'yield_' + proc + '_' + cat for proc in processes for cat in categoriesNames ]
    functions += [ name + str(cat) for name in CBPars for cat in range(0, 13)]
    outName='h013'

commandLine = 'PrintWorkspaceVariables --inFile ' + inFile + ' --outName ' + outName
if len( functions ) : commandLine += ' ' + ' '.join( [ '--function ' + name for name in functions ] )
print commandLine
os.system( commandLine )
