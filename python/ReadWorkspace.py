import os

#mode
#mode/10 : 0->Marc, 1->Christophe
#mode%10 : 0->PrintWorkspace, 1->correlationModel, 2->both
mode=12
inFile=''

functions=[]
processes=['ggH', 'VBF', 'VH', 'ttH', 'bbH', 'tHjb', 'tWH' ]
categoriesNames = ["ggH_FwdLow", 'ggH_FwdHigh','ggH_CenLow', 'ggH_CenHigh', "VBFloose", "VBFtight", "VHMET", "VHlep", "VHdilep", "VHhad_loose", "VHhad_tight",  "ttHhad", "ttHlep"]
CBPars = ['nCBLo_SM_c', 'nCBHi_SM_c', 'alphaCBHi_SM_m125000_c', 'alphaCBLo_SM_m125000_c', 'muCBNom_SM_m125000_c', 'muCBNom_SM_m125000_c' ]

outName='/sps/atlas/c/cgoudet/Plots/'

if ( mode /10 ==  0 ) :
    inFile='/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/h012/StatisticsChallenge/h013/WorkspaceForStudies/mu_i/file_workspace_CATEG_ICHEP_2015_2016.root'
    functions = [ 'Yield_Signal_' + proc + '_SM_'+cat for proc in processes for cat in categoriesNames ]
    functions += [ name + str(cat) for name in CBPars for cat in range(0, 13)]
    outName+='marc'

elif ( mode /10 == 1 ) :
    inFile='/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h013.root'
    functions = [ 'yield_' + proc + '_' + cat for proc in processes for cat in categoriesNames ]#
#    functions += [ name + str(cat) for name in CBPars for cat in range(0, 13)]
    functions += [ '_'.join(['asymParam_ATLAS_BR_gamgam',prod,cat])  for prod in processes for cat in categoriesNames ]
    print( functions )
    outName+='h013'

optionLine=' '
if mode%10==0 or mode%10==2 :
    optionLine += ' --mode printWS '
    optionLine += ' '.join( [ '--function ' + name for name in functions ] ) +' '
if mode%10==1 or mode%10==2 :
    optionLine += ' --mode corrModel '
    

commandLine = 'PrintWorkspaceVariables --inFile ' + inFile + ' --outName ' + outName + optionLine
#print commandLine
os.system( commandLine )
