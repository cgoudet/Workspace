import random
import ROOT
import os
import numpy as np
import subprocess as sub
import argparse

def FitRandomInit( nFits ) : 
    commandLine = 'LikelihoodProfile /sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge013.root  mu_XS_VBF 0.5 1.5 20    --data asimovData_1  --justMin 1  --strategy 1 --profiled mu_XS_ggH --profiled mu_XS_ttH --profiled mu_XS_VH  --scheme 1  '#--outfile /sps/atlas/c/cgoudet/Plots/FitMus'

    for iToy in range( 0, nFits ) :
        toyLine = commandLine + ' --outfile /sps/atlas/c/cgoudet/Hgam/TestInitValueFit/testInit6_' + str(iToy) + '_ '
        varInitToy = [random.random()*2  for i in range(0, 10 ) ]
        toWrite = [ str(iToy) ] + [ str(vVar) for vVar in varInitToy ] 
        toyLine +=  ' --varVal '.join( [''] + [ str(vVal) for vVal in varInitToy] )
        
        os.system( toyLine )

#=======================================
def PlotFits() :
    variables = [ "mu_XS_ggH", "mu_XS_VBF", "mu_XS_VH", "mu_XS_ttH" ]

    listFiles = sub.check_output(['ls /sps/atlas/c/cgoudet/Hgam/TestInitValueFit/*.root'], shell=1, stderr=sub.STDOUT).split()

    boostFileName = '/sps/atlas/c/cgoudet/Hgam/TestInitValueFit/TestInitValueFit.boost'
    boostFile = open( boostFileName, 'w' )
    boostFile.write( 'inputType=1\n' );
    boostFile.write( 'plotDirectory=/sps/atlas/c/cgoudet/Plots/ \n' )
    boostFile.write( '\n'.join( [ 'rootFileName=' + ' '.join( listFiles ) ]*len( variables ) ) + '\n' )
    boostFile.write( '\n'.join( [ 'objName=' + ' '.join( ['nll']*len(listFiles) ) ]*len( variables ) ) + '\n' )
    boostFile.write( '\n'.join( [ 'varName = ' + var  for var in variables ] ) + '\n')
    boostFile.write( 'varMin=0.99\n' )
    boostFile.write( 'varMax=1.01\n' )
    boostFile.write( '\n'.join( [ 'legend='+var.replace('mu_XS_', '' ) for var in variables ] ) + '\n' )
    boostFile.write( 'latex=minimizePrecision : 1e-4\nlatexOpt=0.16 0.9\n' )
    boostFile.close()
    
    os.system( 'CompareHist ' + boostFileName )

    return boostFileName

#=======================================
def parseArgs():
    parser = argparse.ArgumentParser(
        description="This text will be displayed when the -h or --help option are given"
                    "It should contain a short description of what this program is doing.")
    parser.add_argument('--doMode', help='Tag for workspace mode\n',default=1, type=int )
    parser.add_argument('--nFits', help='Tag for workspace mode\n',default=10, type=int )
    args = parser.parse_args()
    return args

#====================
def main():
    args = parseArgs()
    
    if args.doMode == 0 : FitRandomInit( args.nFits )
    elif args.doMode == 1 : print( PlotFits() )
    
# The program entrance
if __name__ == '__main__':
    main()


#     fileName = '/sps/atlas/c/cgoudet/Hgam/TestInitValueFit'
#     f = ROOT.TFile(fileName)
#     t = f.Get( "nll" )
    

#     for entry in t :
#         toWrite.append( str(entry.mu_XS_ggH) )
#         toWrite.append( str(entry.mu_XS_VBF) )
#         toWrite.append( str(entry.mu_XS_VH) )
#         toWrite.append( str(entry.mu_XS_ttH) )


#         outFile.write( ','.join( toWrite ) + '\n' )
# outFile.close()





