import os
import subprocess

input_file='/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h012_asimov.root'

variables = {}
variables["mu_XS_ggH"] = [ 0, 2, 50]    
variables["mu_XS_VBF"] = [ -0.5, 2.5, 50]    
variables["mu_XS_WH"]  = [ -3, 3, 50]    
variables["mu_XS_ZH"]  = [ -3, 3, 50]    
variables["mu_XS_ttH"] = [ -3, 3, 50]    
#variables["mu"] = [ 0, 2]    


fitperjob=50
options = {}
options["justMin"]    = 0
options['saveCsv']    = 1
options['save_np']    = 1
options['constraint'] = 0
options['strategy']   = 1
options['scheme']     = 4
options['data']       = 'obsData'
options['snapshot']   = ''

profiled = ''
#profiled = 'mu_XS_ggH mu_XS_VBF'
#====================================================
def ComputePointVal( varInfo, i ) :
    return varInfo[0] + float( varInfo[1]-varInfo[0] ) / ( varInfo[2]-1 ) * i

#=====================
def StripName( line, doPrefix = 1, doSuffix = 1 ) :
    if ( line.rfind( '.' ) != -1 and doSuffix ) : 
        line = line[0:line.rfind( '.' )]

    if ( line.rfind( '/' ) != -1 and doPrefix ) :
        line = line[line.rfind( '/' )+1:len( line )]

    return line

#==============================

optionLine = ' '.join( [ ( '--' + key + ' ' + str(options[key]) ) if str(options[key]) != '' else '' for key in options.keys() ] )
if profiled != '' : optionLine += ' ' + ' '.join( [ '--profiled ' + var for var in profiled.split(' ') ] ) 

#For 1D plot
commands = (
    [ [' '.join( [ var, str(ComputePointVal( variables[var], iPoint )), str(ComputePointVal( variables[var], min(iPoint+fitperjob, variables[var][-1]-1 ))), str(fitperjob) if iPoint+fitperjob<=variables[var][-1]-1  else str(variables[var][-1]-iPoint) ] ) 
       + ' ' +
       ' '.join( [ '--profiled ' + varP if varP!=var else '' for varP in variables.keys() ] )
       
       , var + '_' + str(int(ComputePointVal( variables[var], iPoint )*1e3))
       ] 
      for var in variables.keys()
      for iPoint in range( 0, variables[var][-1], fitperjob )
      ]
    if  not options["justMin"] else 
    [[ variables.keys()[0] + ' 1 1 1 ' + ' '.join( [ '--profiled ' + varP if varP != variables.keys()[0] else '' for varP in variables.keys() ] ), variables.keys()[0] ]]
    )
print commands

pathResults='/sps/atlas/c/cgoudet/Hgam/Couplages/JobsOutput/'
tag = subprocess.check_output(['date', '+%Y%m%d%H%M%S'])[0:-1]
directory=pathResults + tag + '/'
print directory
subprocess.check_output(['mkdir', directory])

for contents in commands :
    file     = directory + 'launcher_' + contents[1] + '.sh'
    with open( file, 'w' ) as bashFile:
        bashFile.write('server=`pwd`\n' 
                       + 'cd ${server} \n'
                       + 'ulimit -S -s 100000 \n'
                       + 'LD_LIBRARY_PATH=/afs/in2p3.fr/home/c/cgoudet/private/Couplings/RootCoreBin/lib:/afs/in2p3.fr/home/c/cgoudet/private/Couplings/RootCoreBin/bin:$LD_LIBRARY_PATH \n'
                       + 'cd /afs/in2p3.fr/home/c/cgoudet/private/Couplings/RootCoreBin/ \n'
                       + 'source local_setup.sh \n'
                       + 'cd ${server} \n'
                       + 'cp -v /afs/in2p3.fr/home/c/cgoudet/private/Couplings/RootCoreBin/obj/x86_64-slc6-gcc49-opt/PlotFunctions/bin/* . \n'
                       + 'cp ' + input_file + ' .\n'
                       )

        commandLine=[ 'LikelihoodProfile', StripName(input_file, 1, 0) , contents[0], optionLine ]
        bashFile.write( ' '.join( commandLine ) + '\n' )
        bashFile.write( 'rm '+ input_file[input_file.rfind('/')+1:] + '\n' )
        bashFile.write( 'cp *.root ' + directory + '\n')
        bashFile.write( 'cp -v *.csv ' + directory + '\n')
        bashFile.write( 'cp -v *.pdf ' + directory + '\n')

    qsubLine = '~/sub28.sh LP_' + ' '.join( [contents[1], directory + 'log_' + contents[1] + '.log', directory + 'logerror' + contents[1] + '.log', file ] )
    os.system(qsubLine)

