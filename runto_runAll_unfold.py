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
parser.add_argument('--test',help='Run in test mode using smaller number of syst universes (faster)',action='store_true')
parser.add_argument('--fakedata',help='Run with fake data',action='store_true')
parser.add_argument('--makeresponse', help = 'Remake Response Matrix', action = 'store_true')
p = parser.parse_args()

if p.out_dir < 0:
  print "ERROR: Output directory argument not provided"
  parser.print_help()
  exit(1)

## Read parser arguments
outFileDir = p.out_dir
test = ""
fakedata = ""
if p.test > 0:
  test = " --test"
fakedata = ""
if p.fakedata > 0:
  fakedata = " --fakedata"

#############################################################################################################
### Run translateHists.py ###################################################################################
#############################################################################################################

command_string = "python translateHists.py " + outFileDir + test + fakedata
print "Running the following command: \"{0}\"".format(command_string)
output = subprocess.check_output(command_string, shell = True, bufsize = 0)
print output

#############################################################################################################
### Run ResponseMaker.c if requested and get response matrix. ###############################################
#############################################################################################################

responseFileDir = outFileDir[:outFileDir.rfind("/") + 1] + "response_matrices"
if p.makeresponse > 0:
  command_string  = "root -l -q \"ResponseMaker.c(\\\"{0}\\\")\"".format(responseFileDir)
  print "Running the following command: \"{0}\"".format(command_string)
  output = subprocess.check_output(command_string, shell = True, bufsize = 0)
  print output

if not os.path.isfile(responseFileDir + "/response_matrices_exclusive.root"):
  print "ERROR: Response file does not exist! Please rerun with option --makeresponse"
  parser.print_help()
  exit(1)

#############################################################################################################
### Run runAll_unfold.py ####################################################################################
#############################################################################################################

command_string = "python runAll_unfold.py " + outFileDir
print "Running the following command: \"{0}\"".format(command_string)
output = subprocess.check_output(command_string, shell = True, bufsize = 0)
print output
