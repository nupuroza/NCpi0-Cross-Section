import os         ## to check if response file exists
import argparse
import subprocess ## to make command line calls
import re         ## to pull out date from in_dir
import shutil     ## to copy response file to output directory

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################

parser = argparse.ArgumentParser(description='Script to run all steps of unfolding code up to runAll_unfold.py')
parser.add_argument('out_dir', help='Path to output directory', type=str,nargs='?')
parser.add_argument('server', help='Either gpvm or manannan', type=str,nargs='?')
parser.add_argument('user', help = 'The name of your user directory', type = str, nargs = '?')
parser.add_argument('reco_var_2g1p_input', help = '2g1p reco variable definition', type = str, nargs = '?')
parser.add_argument('reco_var_2g0p_input', help = '2g0p reco variable definition', type = str, nargs = '?')
parser.add_argument('true_var_input', help = 'True variable definition', type = str, nargs = '?')
parser.add_argument('--translateHiststest',help='Run in test mode using smaller number of syst universes (faster)',action='store_true')
parser.add_argument('--fakedata',help='Run with fake data',action='store_true')
parser.add_argument('--closureTest',help='Run as closure test',action='store_true')
parser.add_argument('--makeresponse', help = 'Remake Response Matrix', action = 'store_true')
parser.add_argument('--plotstest',help='Produce plots in test mode; useful for debugging',action='store_true')
p = parser.parse_args()

if p.out_dir < 0:
  print("ERROR: Output directory argument not provided")
  parser.print_help()
  exit(1)

## Read parser arguments
outFileDir = p.out_dir
if len(outFileDir) - 1 == outFileDir.rfind("/"):
  outFileDir = outFileDir[:-1] # Remove trailing forward slash.
server = p.server
user = p.user
translateHiststest = ""
fakedata = ""
reco_var_2g1p_input = ""
reco_var_2g0p_input = ""
true_var_input = ""
if p.reco_var_2g1p_input > 0:
  reco_var_2g1p_input = p.reco_var_2g1p_input
if p.reco_var_2g0p_input > 0:
  reco_var_2g0p_input = p.reco_var_2g0p_input
if p.true_var_input > 0:
  true_var_input = p.true_var_input
closureTest = ""
plotstest = ""
if p.translateHiststest > 0:
  translateHiststest = " --test"
if p.fakedata > 0:
  fakedata = " --fakedata"
if p.closureTest > 0:
  closureTest = " --closureTest"
if p.plotstest > 0:
  plotstest = " --test"

#############################################################################################################
### Run ResponseMaker.c if requested and get response matrix. ###############################################
#############################################################################################################

responseFileDir = outFileDir[:outFileDir.rfind("/") + 1] + "response_matrices"
if p.makeresponse > 0:
  command_string  = "root -l -q \"ResponseMaker.c(\\\"{0}\\\", \\\"{1}\\\", \\\"{2}\\\", \\\"{3}\\\", \\\"{4}\\\", \\\"{5}\\\")\"".format(responseFileDir, server, user, reco_var_2g1p_input, reco_var_2g0p_input, true_var_input)
  print("Running the following command: \"{0}\"".format(command_string))
  output = subprocess.check_output(command_string, shell = True, bufsize = 0)
  print(output)

if not os.path.isfile(responseFileDir + "/response_matrices_exclusive.root"):
  print "ERROR: Response file does not exist! Please rerun with option --makeresponse"
  print("responseFileDir is " + responseFileDir)
  parser.print_help()
  exit(1)

#############################################################################################################
### Run translateHists.py ###################################################################################
#############################################################################################################

command_string = "python translateHists.py " + outFileDir + " " + server + " " + user + translateHiststest + fakedata
print( "Running the following command: \"{0}\"".format(command_string))
output = subprocess.check_output(command_string, shell = True, bufsize = 0)
print(output)

#############################################################################################################
### Run runAll_unfold.py ####################################################################################
#############################################################################################################

command_string = "python runAll_unfold.py " + outFileDir + closureTest
print("Running the following command: \"{0}\"".format(command_string))
output = subprocess.check_output(command_string, shell = True, bufsize = 0)
print(output)

#############################################################################################################
### Run runAll_calculateXsection.py #########################################################################
#############################################################################################################

command_string = "python runAll_calculateXsection.py " + outFileDir + closureTest
print("Running the following command: \"{0}\"".format(command_string))
output = subprocess.check_output(command_string, shell = True, bufsize = 0)
print(output)

#############################################################################################################
### Run runAll_makeFakeDataStudyPlots.py ####################################################################
#############################################################################################################

command_string = "python runAll_makeFakeDataStudyPlots.py " + outFileDir + closureTest + plotstest
print("Running the following command: \"{0}\"".format(command_string))
output = subprocess.check_output(command_string, shell = True, bufsize = 0)
print(output)
