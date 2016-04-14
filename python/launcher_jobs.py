import os
import subprocess

path='/afs/in2p3.fr/home/c/cgoudet/private/Couplings/'
data='/sps/atlas/c/cgoudet/Hgam/Couplages/JobsOutput/'
tag = subprocess.check_output(['date', '+%Y%m%d%H%M%S'])[0:-1]
print tag

#MENU
#input_file='/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_Test12cat.root'
input_file='/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h011_test.root'
#input_file='/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h011_pdfReco.root'
#input_file='/sps/atlas/c/cgoudet/Hgam/Couplages/Outputs/StatChallenge_h011_fullreco.root'

dataset='--data obsData_G'
options='--saveCsv --save_np'

#variable1='mu'
variable1='mu_XS_ggH'
nx=50
xmin=1
xmax=1

#Set to '' to make true 1D profile
variable2=''
#Set to 1 to make 1D profile with variable 2 as free parameter 
ny=1 
ymin=0.6
ymax=1.4
#Put profiled variables
#variables=[]
variables=['mu_XS_VBF', 'mu_XS_WH', 'mu_XS_ZH', 'mu_XS_ttH']
#['mu_XS_WH', 'mu_XS_ZH', 'mu_XS_ggH', 'mu_XS_ttH', 'mHcomb'] 

strategy='--strategy 1'
fitperjob=5
justMin=1
#'-m 1' #Choose 0 if no specific changes to do
modif_scheme=' --scheme 0'
#'--Snapshot abcx'  Give the name of snapshot for asimov
snapshot=''

#########################################################################
########################################################################

#create code line for profiled variables
v3line=''
for v in variables:
 v3line+= ' --profiled ' + v

#Create the temporary file to store all .root files
directory=data + tag + '/'
print directory
subprocess.check_output(['mkdir', directory])

#create line for variable2
v2line = '' if variable2=='' else variable2 + ' ' + str(ymin) + ' ' + str(ymax) + ' ' + str(ny) 

if nx==1: valueVar1 = [xmin]
else: valueVar1=[ xmin + ( xmax-xmin ) / ( nx-1 ) * var1 for var1 in range( nx ) ]

valueVar2 = [0] if variable2=='' else ([ ymin ] if ny==1 else [ymin + ( ymax-ymin ) / ( ny-1 ) * var2 for var2 in range( ny ) ] )

if justMin==1 : fitperjob=1
if variable1=='' : ny=1

for step1 in range( -1, 0 if justMin==1 else nx , fitperjob ):
    for step2 in range( 0,  ny if step1 != -1 else 1 ):
        

        file     = directory + 'launcher_' + str(valueVar1[step1]) +  '_' + str(valueVar2[step2]) + '.sh'
        log      = directory + 'log_' + str(valueVar1[step1]) +  '_' +  str(valueVar2[step2]) + '.log'
        logerror = directory + 'logerror_' + str(valueVar1[step1]) +  '_' +  str(valueVar2[step2]) + '.log'
        
        v2line = '' if variable2 == '' else variable2 + ' ' + str(valueVar2[step2]) + ' ' + str(ymax) + ' ' + str(ny) + ' 1'
        v1line   = ' '.join( [ variable1, str(valueVar1[step1]), str(valueVar1[min(step1+fitperjob-1, len( valueVar1 )-1 )]) , str(fitperjob), v2line, v3line ] )
        
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

            commandLine=[ 'LikelihoodProfile', input_file[input_file.rfind('/')+1:] , dataset, options, modif_scheme, strategy, v1line, snapshot, '' if justMin!=1 else ' --justMin ' ]
            bashFile.write( ' '.join( commandLine ) + '\n' )
            bashFile.write( 'rm '+ input_file[input_file.rfind('/')+1:] + '\n' )
            bashFile.write( 'cp *.root ' + directory + '\n')
            bashFile.write( 'cp -v *.csv ' + directory + '\n')
            bashFile.write( 'cp -v *.pdf ' + directory + '\n')
            
        qsubLine = 'qsub -P P_atlas -m e -M goudet@lal.in2p3.fr -l ct=29:00:00 -l vmem=4G,fsize=10G -l sps=1 -N LP_' + str(step1) + '_' + str(nx) + '_' + str(step2) + '_' + str(ny) + ' -o ' + log + ' -e ' + logerror + ' ' + file
        os.system(qsubLine)

