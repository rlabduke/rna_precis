from __future__ import absolute_import, division, print_function
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.validation import restraints
from libtbx import easy_run
import os, sys

def run_readyset(pdbfile):
  tempdir = "readyset_temp"
  #libdir = "readyset_library"
  if not os.path.isdir(tempdir):
    os.mkdir(tempdir)
  #if not os.path.isdir(libdir):
  #  os.mkdir(libdir)
  os.chdir(tempdir)
  ready_set_command = "phenix.ready_set ../%s" % (pdbfile)
  ready_set_out = easy_run.fully_buffered(ready_set_command)

  cif_fnames = os.listdir('.')
  cif_objects = prepare_restaints(cif_fnames)
  return cif_objects

def cleanup_readyset():
  tempdir = "readyset_temp"
  #Move new restraints files to libdir, delete all files in tempdir
  filelist = os.listdir(tempdir)
  #liblist = os.listdir(libdir)
  for filename in filelist:
    filepath = os.path.join(tempdir,filename)
    os.remove(filepath)
  os.removedirs(tempdir)
    #if len(filename) <= 7 and filename.endswith('.cif') and filename not in liblist:
    #  os.rename(filepath,os.path.join(libdir,filename))
    #else:
    #  os.remove(filepath)
  #could also remove the tempdir, I think I'll leave it as evidence for now.
  #os.remove(tempdir)

def prepare_restaints(cif_filenames):
  cif_objects = list()
  for cif_file_name in cif_filenames:
    if len(cif_file_name) > 7: continue
    if cif_file_name.endswith('.cif'):
      cif_object = iotbx.cif.reader(file_path=cif_file_name).model()
      cif_objects.append((cif_file_name, cif_object))
  return cif_objects

def find_atom_xyz(rg,atomname):
  #check the '' and 'A' alternates for an atoms of a given name
  for ag in rg.atom_groups():
    if ag.altloc not in ['','A']: continue
    for atom in ag.atoms():
      if atom.name == atomname:
        return atom.xyz
  #if rg contains no matching atom
  return None

def gamma_from_atoms(pdb_hierarchy,chainID,resnum,i_code):
  #Many modified NAs do not have gamma/delta specified in their restraints
  #given residues identifiers, search hierarchy for the residue
  #Search for required atoms (gamma = ["O5'","C5'","C4'","C3'"])
  #calculate gamma
  #return gamma = '__?__' on any failure
  from cctbx import geometry_restraints
  for chain in pdb_hierarchy.chains():
    if chain.id.strip() != chainID.strip(): continue
    for rg in chain.residue_groups():
      if rg.resseq == resnum  and rg.icode == i_code:
        gamma_atoms = [find_atom_xyz(rg," O5'"),
                       find_atom_xyz(rg," C5'"),
                       find_atom_xyz(rg," C4'"),
                       find_atom_xyz(rg," C3'")]
        if None in gamma_atoms: return '__?__'
        gamma = geometry_restraints.dihedral(sites=gamma_atoms,angle_ideal=180, weight=1).angle_model
        return "%.3f" % gamma
  #if no matching residue found
  return '__?__'

def delta_from_atoms(pdb_hierarchy,chainID,resnum,i_code):
  #Many modified NAs do not have gamma/delta specified in their restraints
  #given residues identifiers, search hierarchy for the residue
  #Search for required atoms (delta = ["C5'","C4'","C3'","O3'"])
  #calculate delta
  #return delta = '__?__' on any failure
  from cctbx import geometry_restraints
  for chain in pdb_hierarchy.chains():
    if chain.id.strip() != chainID.strip(): continue
    for rg in chain.residue_groups():
      if rg.resseq == resnum  and rg.icode == i_code:
        delta_atoms = [find_atom_xyz(rg," C5'"),
                       find_atom_xyz(rg," C4'"),
                       find_atom_xyz(rg," C3'"),
                       find_atom_xyz(rg," O3'")]
        #print(delta_atoms)
        if None in delta_atoms: return '__?__'
        delta = geometry_restraints.dihedral(sites=delta_atoms,angle_ideal=180, weight=1).angle_model
        return "%.3f" % delta
  #if no matching residue found
  return '__?__'    

##Copied from mmtbx/validation/utils, then edited for local inputs
#This makes suitename input-formatted rna dihedral params
#alpha, beta, gamma, delta, epsilon, zeta
# :1: A:  14: : :  A:83.353:-158.279:-114.591:92.038:-125.519:-57.251
# :1: A:  15: : :  G:-55.064:162.460:51.892:79.826:-136.303:-143.947
# :1: A:  16: : :H2U:-6.131:91.184:__?__:__?__:-61.796:-131.179
################################################################################
def get_inverted_atoms(atoms, improper=False):
  temp = []
  if not improper:
    temp.append(atoms[3])
    temp.append(atoms[2])
    temp.append(atoms[1])
    temp.append(atoms[0])
  else:
    temp.append(atoms[3])
    temp.append(atoms[1])
    temp.append(atoms[2])
    temp.append(atoms[0])
  return temp

def match_dihedral_to_name(atoms):
  name = None
  alpha = ["O3'","P","O5'","C5'"]
  beta = ["P","O5'","C5'","C4'"]
  gamma = ["O5'","C5'","C4'","C3'"]
  delta = ["C5'","C4'","C3'","O3'"]
  epsilon = ["C4'","C3'","O3'","P"]
  zeta = ["C3'","O3'","P","O5'"]
  if atoms == alpha:
    name = "alpha"
  elif atoms == beta:
    name = "beta"
  elif atoms == gamma:
    name = "gamma"
  elif atoms == delta:
    name = "delta"
  elif atoms == epsilon:
    name = "epsilon"
  elif atoms == zeta:
    name = "zeta"
  return name

def build_name_hash(pdb_hierarchy):
  i_seq_name_hash = dict()
  for atom in pdb_hierarchy.atoms():
    i_seq_name_hash[atom.i_seq]=atom.pdb_label_columns()
  return i_seq_name_hash

def get_rna_backbone_dihedrals(pdb_hierarchy=None, geometry=None):
  # at present, this will only return measurements for angles arising from
  # atoms with altloc ' ' or altloc 'A'
  # TO-DO: extend to more alternates JJH 140108
  from cctbx import geometry_restraints
  from collections import defaultdict
  bb_dihedrals = defaultdict(dict)
  formatted_out = []
  alt_tracker = {}
  ##if (processed_pdb_file is not None):
  ##  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart
  ##  geometry = processed_pdb_file.geometry_restraints_manager()
  ##  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  ##else :
  ##  assert (not None in [geometry, pdb_hierarchy])
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  dihedral_proxies = geometry.dihedral_proxies
  i_seq_name_hash = build_name_hash(pdb_hierarchy=pdb_hierarchy)

  def is_blank_or_alt_a(proxy):
    for i in proxy.i_seqs:
       alt = i_seq_name_hash[i][4:5]
       if alt not in [' ', 'A']:
         return False
    return True

  for dp in dihedral_proxies:
    atoms = []
    debug_key = ""
    invert_sign = False
    dp.sort_i_seqs()
    for i in dp.i_seqs:
      atoms.append(i_seq_name_hash[i][0:4].strip())
      debug_key = debug_key+i_seq_name_hash[i]
    if len(atoms) != 4:
      continue
    name = match_dihedral_to_name(atoms=atoms)
    #handle dihedral equivalences
    if name == None:
      inverted_atoms = get_inverted_atoms(atoms=atoms, improper=False)
      name = match_dihedral_to_name(atoms=inverted_atoms)
      if name == None:
        inverted_atoms = get_inverted_atoms(atoms=atoms, improper=True)
        name = match_dihedral_to_name(atoms=inverted_atoms)
        if name is not None:
          invert_sign = True
    if (name is not None) and (is_blank_or_alt_a(dp)):
      restraint = geometry_restraints.dihedral(
                                               sites_cart=sites_cart,
                                               proxy=dp)
      key = i_seq_name_hash[dp.i_seqs[1]][4:]
      if alt_tracker.get(key[1:]) is None:
        alt_tracker[key[1:]] = []
      if key[0:1] not in alt_tracker[key[1:]]:
        alt_tracker[key[1:]].append(key[0:1])
      bb_dihedrals[key][name] = restraint.angle_model
      if invert_sign:
        bb_dihedrals[key][name] = bb_dihedrals[key][name] * -1.0
  for key in list(bb_dihedrals.keys()):
    altloc = key[0:1]
    resname = key[1:4]
    chainID = key[4:6]
    resnum = key[6:10]
    i_code = key[10:]
    if 'A' in alt_tracker[key[1:]]:
      if altloc != 'A':
        continue
    if bb_dihedrals[key].get('alpha') is not None:
      alpha = "%.3f" % bb_dihedrals[key]['alpha']
    # FIXME will the lookup below ever actually work?
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('alpha') is not None:
      alpha = "%.3f" % bb_dihedrals[' '+key[1:]]['alpha']
    else:
      alpha = '__?__'
    if bb_dihedrals[key].get('beta') is not None:
      beta = "%.3f" % bb_dihedrals[key]['beta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('beta') is not None:
      beta = "%.3f" % bb_dihedrals[' '+key[1:]]['beta']
    else:
      beta = '__?__'
    if bb_dihedrals[key].get('gamma') is not None:
      gamma = "%.3f" % bb_dihedrals[key]['gamma']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('gamma') is not None:
      gamma = "%.3f" % bb_dihedrals[' '+key[1:]]['gamma']
    else:
      ###Added by CJW, attempt to calculate gamma from atoms rather than proxies
      gamma = gamma_from_atoms(pdb_hierarchy,chainID,resnum,i_code)
      #gamma = '__?__'
    if bb_dihedrals[key].get('delta'):
      delta = "%.3f" % bb_dihedrals[key]['delta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('delta') is not None:
      delta = "%.3f" % bb_dihedrals[' '+key[1:]]['delta']
    else:
    ###Added by CJW, attempt to calculate delta from atoms rather than proxies
      delta = delta_from_atoms(pdb_hierarchy,chainID,resnum,i_code)
      #delta = '__?__'
    if bb_dihedrals[key].get('epsilon'):
      epsilon = "%.3f" % bb_dihedrals[key]['epsilon']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('epsilon') is not None:
      epsilon = "%.3f" % bb_dihedrals[' '+key[1:]]['epsilon']
    else:
      epsilon = '__?__'
    if bb_dihedrals[key].get('zeta'):
      zeta = "%.3f" % bb_dihedrals[key]['zeta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('zeta') is not None:
      zeta = "%.3f" % bb_dihedrals[' '+key[1:]]['zeta']
    else:
      zeta = '__?__'
    eval = "%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s" \
           % (" ",
              "1",
              chainID,
              resnum,
              i_code,
              altloc,
              resname,
              alpha,
              beta,
              gamma,
              delta,
              epsilon,
              zeta)
    formatted_out.append(eval)
  formatted_out.sort()
  backbone_dihedrals = ""
  for line in formatted_out:
    backbone_dihedrals += line+'\n'
  return backbone_dihedrals
################################################################################

def parse_suitename_output(suite_lines):
#exec("phenix.suitename -report -pointIDfields 7 -altIDfield 6 < $midfile > $outfile");
# :1: A:  13: : :  C 33 t 1c 0.824
# :1: A:  14: : :  A trig !! 0.000 gamma
# :1: A:  15: : :  G 33 p 1a 0.484
# :1: A:  16: : :H2U inc  __ 0.000
#...
#Found 52 complete suites derived from 76 entries
#1 suites were triaged, leaving 51 assigned to bins
#For all 51 suites: average suiteness== 0.750 (power==3.00)
#     2 suites are  outliers
#     0 suites have suiteness == 0   
  suite_results = {}
  for line in suite_lines:
    if not line.startswith(" :"):
      continue
    x = line.split(":")
    resid = ":".join(x[0:6])
    s = x[6]
    resname = s[0:3]
    cat = s[4:8]
    suite = s[9:11]
    suiteness = s[12:17]
    triage = s[18:].strip()
    suite_results[resid] = ":".join([cat,suite,suiteness,triage])
  return suite_results


def run(pdbfile, do_elbow=True, do_cleanup=True):
  pdb_inp = iotbx.pdb.input(file_name = pdbfile)

  if do_elbow:
    cif_objects = run_readyset(pdbfile)
  else:
    cif_objects=None


  model = mmtbx.model.manager(
    model_input = pdb_inp,
    build_grm   = True,
    restraint_objects = cif_objects,
    stop_for_unknowns = True,
    log = sys.stderr)
  pdb_hierarchy = model.get_hierarchy()
  grm = model.get_restraints_manager().geometry
  xray_structure = model.get_xray_structure()

  rna_bb = get_rna_backbone_dihedrals(pdb_hierarchy=pdb_hierarchy, geometry=grm)
  #This is just one big text string

  suitename_output = easy_run.fully_buffered("phenix.suitename -report -pointIDfields 7 -altIDfield 6 -", stdin_lines=rna_bb)
  suitename_results = parse_suitename_output(suitename_output.stdout_lines)

  result = []
  result.append("blank:1:chainID:resnum:i_code:altloc:resname:alpha:beta:gamma:delta:epsilon:zeta:bin:suitename:suiteness:triage")
  for line in rna_bb.split("\n"):
    if not line: continue
    resid = ":".join(line.split(":")[0:6])
    try:
      suite = suitename_results[resid]
    except KeyError:
      suite = "NULL"
    result.append(line+":"+suite)

  if do_elbow:
    os.chdir("..")
    if do_cleanup:
      cleanup_readyset()

  return result


if __name__ == "__main__":
#input validation
  #Parse cmdline
  #if len(sys.argv) > 3:
  #  print("Needs one or two arguments: a pdb/cif structure file and an optional restraints dir")
  #  sys.exit()
  pdbfile = sys.argv[1].strip()
  if len(os.path.basename(pdbfile)) <= 7:
    print("Short filenames will interfere with ligand restraint generation.")
    print("Please use a filename at least 8 total characters (4 before file extention)")
    print("  Example: 1234.pdb or longer")
    sys.exit()
  if len(sys.argv) != 2:
    print("This script accepts a single pdb or mmcif formatted file")
    sys.exit()
  if sys.argv[1][-4:] not in [".pdb",".cif"]:
    print("This script accepts a single pdb or mmcif formatted file")
    sys.exit()
#end input validation

  infilepath = sys.argv[1]
  result = run(infilepath)
  for line in result:
    print(line)
