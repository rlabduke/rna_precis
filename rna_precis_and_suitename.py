from __future__ import absolute_import, division, print_function
import os, sys
try:
  from libtbx import easy_run
  import iotbx.phil
except:
  sys.stderr.write("""
  These scripts must be run in a Phenix environment
  run 'source PATH_to_phenix/builds/setpaths.sh'
  then try again""")
  sys.exit()
try:
  from rna_scripts import rna_precis_parameters, rna_suite_parameters
except:
  print("""
  RNA processing scripts could not be found.
  Try running with rna_precis_params as the working dir.
  Or contact Christopher.""", file=sys.stderr)
  sys.exit()

## Commandline parsing ##
def get_master_phil():
  return iotbx.phil.parse(input_string="""
    elbow = True
      .type = bool
      .help = run readyset to generate restraint files for non-standard residues
    complete_only = False
      .type = bool
      .help = only print output for residues with complete parameterization
    model = None
      .type = path
      .help = the input file, in pdb or mmcif format
    cleanup_after_readyset=True
      .type = bool
      .help = remove readyset_temp dir and contents, set False for debugging
""", process_includes=True)

usage_string = """This script generates suitename parameters, suitename analysis, and RNAprecis parameters"""

args = sys.argv[1:]

cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=get_master_phil(),
    pdb_file_def="model",
    usage_string=usage_string)
params = cmdline.work.extract()
pdbfile = params.model

if not pdbfile.endswith(".pdb") and not pdbfile.endswith(".cif") and not pdbfile.endswith(".ent"):
  usage()
  sys.exit()
## End Commandline parsing ##

## Run scripts to generate data ##
print("Running Readyset and suite parameters, this may take a while . . .", file=sys.stderr)
suite_params = rna_suite_parameters.run(pdbfile, do_elbow=params.elbow, do_cleanup=params.cleanup_after_readyset)
print("Running new rna parameters . . .", file=sys.stderr)
new_rna_params = rna_precis_parameters.run(pdbfile)
## Scripts done ##

## match outputs with each other
residues = {}

suite_header = suite_params[0]
for line in suite_params[1:]: #slice drops the header line. Can't believe it took me so long to figure that out.
  #blank:1:chainID,resnum,i_code,altloc,resname,alpha,beta,gamma,delta,epsilon,zeta
  # :1: A:   6: : :  U:__?__:157.567:55.308:81.262:-150.997:-77.008
  # :1: A:   7: : :  U:-59.498:-178.720:62.506:137.257:-105.934:-51.960
  if not line.strip(): continue
  x = line.split(":")
  idstr = ":".join(x[2:7])
  data = ":".join(x[7:])
  residues[idstr] = {"suite":None,"new":None}
  residues[idstr]["suite"] = data

new_rna_header = new_rna_params[0]
for line in new_rna_params[1:]:
  #blank:1:chainID:resnum:i_code:altloc:resname:Pperp1:d1'1':ang1:ang2:dhNN:pto1'1':Pdown:dhPN:Pperp2
  # :1: A:   7: : :  U:4.508:5.039:96.237:115.345:9.002:4.722:2.381:54.401:1.414
  # :1: A:   8: : :  U:1.414:9.429:43.978:88.650:94.335:3.223:4.262:84.523:4.762
  if not line.strip(): continue
  x = line.split(":")
  idstr = ":".join(x[2:7])
  data = ":".join(x[7:])
  if idstr not in residues:
    residues[idstr] = {"suite":None,"new":None}
  residues[idstr]["new"] = data
## matchup done ##

## output the combined results
print("chainID:resnum:i_code:altloc:resname:alpha:beta:gamma:delta:epsilon:zeta:bin:suitename:suiteness:triage:Pperp1:d1'1':ang1:ang2:dhNN:pto1'1':Pdown:dhPN:Pperp2")
reskeys = residues.keys()
reskeys.sort()
for idstr in reskeys:
  suite = residues[idstr]["suite"]
  new = residues[idstr]["new"]
  if params.complete_only:
    if suite is None or new is None:
      continue
    if "__?__" in suite or "__?__" in new:
      continue
  if suite is None:
    suite = "__?__:__?__:__?__:__?__:__?__:__?__"
  if new is None:
    new = "__?__:__?__:__?__:__?__:__?__:__?__:__?__:__?__:__?__"
  print(":".join([idstr,suite,new]))
## output done ##

