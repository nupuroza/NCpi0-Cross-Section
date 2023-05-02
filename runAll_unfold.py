import argparse
import subprocess ## to make command line calls
import re         ## to pull out date from in_dir

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################

parser = argparse.ArgumentParser(description='Script to run unfold.C for all configurations of unfolder')
parser.add_argument('in_dir', help='Path to input directory', type=str,nargs='?')
p = parser.parse_args()

## If in_dir is not provided, exit
if p.in_dir < 0:
  print "ERROR: Input directory argument not provided"
  parser.print_help()
  exit(1)

inFileDir = p.in_dir

## Pull date out of in_dir
date_pattern = r"\d{4}-\d{2}-\d{2}"
match = re.search(date_pattern, inFileDir)  ## search for the pattern in the string
date_string = match.group()                ## extract the matched date

#############################################################################################################
### Run CalculateXsection.py ################################################################################
#############################################################################################################

configs_to_run = [
  ["WSVD-kIdentity","true","kIdentity"],
  ["WSVD-kFirstDerivative","true","kFirstDerivative"],
  ["WSVD-kSecondDerivative","true","kSecondDerivative"],
  ["DAgostini-1-iteration","false","1"],
  ["DAgostini-2-iteration","false","2"],
  ["DAgostini-3-iteration","false","3"],
  ["DAgostini-4-iteration","false","4"],
  ["DAgostini-5-iteration","false","5"]
]

for config,WSVD_switch,unfolder_option in configs_to_run:

  file_to_process = "{0}/{1}_out.root".format(inFileDir,date_string)
  command_string = "root -l -q \"unfold.C(\\\"{0}\\\",{1},\\\"{2}\\\",false)\"".format(file_to_process,WSVD_switch,unfolder_option)
  print "Running the following command: \"{0}\"".format(command_string)

  ## Run the command and capture the output
  output = subprocess.check_output(command_string, shell=True, bufsize=0)
  print(output)

