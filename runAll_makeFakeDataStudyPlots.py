import argparse
import subprocess ## to make command line calls
import re         ## to pull out date from in_dir

#############################################################################################################
### File Management #########################################################################################
#############################################################################################################

parser = argparse.ArgumentParser(description='script to make fake data study plots for all configurations of unfolder')
parser.add_argument('in_dir', help='Path to input directory', type=str,nargs='?')
parser.add_argument('--closureTest',help='Input file corresponds to closure test',action='store_true')
parser.add_argument('--test', help = 'Produce plots in test mode; useful for debugging', action = 'store_true')
p = parser.parse_args()

## If in_dir is not provided, exit
if p.in_dir < 0:
  print "ERROR: Input directory argument not provided"
  parser.print_help()
  exit(1)

inFileDir = p.in_dir

closureTest = "--closureTest" if p.closureTest>0 else ""
closureTest_ = "closureTest_" if p.closureTest > 0 else ""
test = "--test" if p.test > 0 else ""

## Pull date out of in_dir
date_pattern = r"\d{4}-\d{2}-\d{2}"
match = re.search(date_pattern, inFileDir)  ## search for the pattern in the string
date_string = match.group()                ## extract the matched date

#############################################################################################################
### Run makeFakeDataStudyPlots.py ###########################################################################
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

  file_to_process = "{0}/{1}_out_unfolded_{2}_{3}xsec-extracted.root".format(inFileDir,date_string,config, closureTest_)
  plotDir = "{0}/{1}_xsec-plots_{2}/".format(inFileDir,date_string,config)
  command_string = "python makeFakeDataStudyPlots.py {0} {1} {2} {3}".format(file_to_process,plotDir, closureTest, test)
  print "Running the following command: \"{0}\"".format(command_string)

  ## Run the command and capture the output
  output = subprocess.check_output(command_string, shell=True, bufsize=0)
  print(output)

