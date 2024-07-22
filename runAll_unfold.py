import argparse
import subprocess ## to make command line calls
import re         ## to pull out date from in_dir
import os

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################

parser = argparse.ArgumentParser(description='Script to run unfold.C for all configurations of unfolder')
parser.add_argument('in_dir', help='Path to input directory', type=str,nargs='?')
parser.add_argument('--closureTest',help='Input file corresponds to closure test',action='store_true')
p = parser.parse_args()

## If in_dir is not provided, exit
if p.in_dir < 0:
  print "ERROR: Input directory argument not provided"
  parser.print_help()
  exit(1)

inFileDir = p.in_dir

if inFileDir[-1] == '/':
  inFileDir = inFileDir[:-1]

is_closure_test = "true" if p.closureTest>0 else "false"

## Pull date out of in_dir
date_pattern = r"\d{4}-\d{2}-\d{2}"
match = re.search(date_pattern, inFileDir)  ## search for the pattern in the string
date_string = match.group()                ## extract the matched date

## Remove previous cumulative stats files if present
file_to_remove_2g1p = inFileDir + "/" + date_string + "_out_unfolded_stats_2g1p.txt"
file_to_remove_2g0p = inFileDir + "/" + date_string + "_out_unfolded_stats_2g0p.txt"
if os.path.exists(file_to_remove_2g1p):
  os.remove(file_to_remove_2g1p)
if os.path.exists(file_to_remove_2g0p):
  os.remove(file_to_remove_2g0p)

#############################################################################################################
### Run CalculateXsection.py ################################################################################
#############################################################################################################

configs_to_run = [
  ["WSVD-kIdentity","true","kIdentity"],
  ["WSVD-kFirstDerivative","true","kFirstDerivative"],
  ["WSVD-kSecondDerivative","true","kSecondDerivative"],
  ]

for i in range(5):
    configs_to_run.append(["DAgostini-{}-iteration".format(i+1), "false", "{}".format(i+1)])

for config,WSVD_switch,unfolder_option in configs_to_run:

  file_to_process = "{0}/{1}_out.root".format(inFileDir,date_string)
  command_string = "root -l -q \"unfold.C(\\\"{0}\\\",{1},\\\"{2}\\\", {3})\"".format(file_to_process,WSVD_switch,unfolder_option, is_closure_test)
  print "Running the following command: \"{0}\"".format(command_string)

  ## Run the command and capture the output
  output = subprocess.check_output(command_string, shell=True, bufsize=0)
  print(output)

