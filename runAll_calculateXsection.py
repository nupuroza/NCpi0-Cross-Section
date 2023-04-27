import argparse
import subprocess ## to make command line calls
import re         ## to pull out date from in_dir

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################

parser = argparse.ArgumentParser(description='Script to perform the cross-section calculation for all configurations of unfolder')
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
  "WSVD-kIdentity",
  "WSVD-kFirstDerivative",
  "WSVD-kSecondDerivative",
  "DAgostini-1-iteration",
  "DAgostini-2-iteration",
  "DAgostini-3-iteration",
  "DAgostini-4-iteration",
  "DAgostini-5-iteration"
]

for config in configs_to_run:

  file_to_process = "{0}/{1}_out_unfolded_{2}.root".format(inFileDir,date_string,config)
  command_string = "python calculateXsection.py {0}".format(file_to_process)
  print "Running the following command: \"{0}\"".format(command_string)

  ## Run the command and capture the output
  output = subprocess.check_output(command_string, shell=True, bufsize=0)
  print(output)

