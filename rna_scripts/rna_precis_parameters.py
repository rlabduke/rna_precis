from __future__ import absolute_import, division, print_function
#from mmtbx.conformation_dependent_library import generate_dna_rna_fragments
from cctbx import geometry_restraints
from libtbx import easy_run
import iotbx.pdb
import mmtbx.model
import os, sys

#return geometry_restraints.dihedral(sites=[CA1.xyz,C1.xyz,N2.xyz,CA2.xyz],
#  angle_ideal=180, weight=1).angle_model

class rna_params():
  def __init__(self, three):
    self.three = three
    self.p_im1, self.n_im1, self.c_im1 = self.find_pnc(three[0])
    self.p_i, self.n_i, self.c_i = self.find_pnc(three[1])
    self.p_ip1, self.n_ip1, self.c_ip1 = self.find_pnc(three[2])
    if None in [self.c_im1, self.n_im1,
                self.p_i, self.c_i, self.n_i,
                self.p_ip1]:
      self.complete = False
    else:
      self.complete = True
    if self.complete: self.do_params()

  def find_pnc(self,res):
    p = find_atom(res, " P  ")
    c = find_atom(res, " C1'")
    #need N9 from A/G, need N1 from C/U, need ? from modified base
    #check N9 then N1 to see if one is within bonding distance of C1' 
    n = find_atom(res, " N9 ")
    if not are_glyco_bonded(n,c):
      n = find_atom(res, " N1 ")
      if not are_glyco_bonded(n,c):
        n = None
    return p,n,c

  def do_pperp(self, p,n,c):
    intersect_point = perptersect(n.xyz,c.xyz,p.xyz)
    d = distance(intersect_point,p.xyz)
    return d

  def do_params(self):
    #perp distance from P(n) to the N(-1) to C1'(-1) line
    self.pperp1 = self.do_pperp(self.p_i, self.n_im1, self.c_im1)
    #distance from C1'(-1) to C1'(n)
    self.d11 = distance(self.c_im1.xyz, self.c_i.xyz)
    self.ang1 = geometry_restraints.angle(sites=[self.n_im1.xyz, self.c_im1.xyz, self.c_i.xyz],
      angle_ideal=120, weight=1).angle_model
    self.ang2 = geometry_restraints.angle(sites=[self.c_im1.xyz, self.c_i.xyz, self.n_i.xyz],
      angle_ideal=120, weight=1).angle_model
    self.dhNN = geometry_restraints.dihedral(sites=[self.n_im1.xyz, self.c_im1.xyz, self.c_i.xyz, self.n_i.xyz],
      angle_ideal=180, weight=1).angle_model
    #distance from P(n) to perp intersect x along C1'(-1) to C1'(n) line
    x = perptersect(self.c_im1.xyz, self.c_i.xyz, self.p_i.xyz)
    self.pto11 = distance(x, self.p_i.xyz)
    #distance from that intersect x, to C1'(n)
    self.pdown = distance(x, self.c_i.xyz)
    self.dhPN = geometry_restraints.dihedral(sites=[self.p_i.xyz, x, self.c_i.xyz, self.n_i.xyz],
      angle_ideal=180, weight=1).angle_model
    self.pperp2 = self.do_pperp(self.p_ip1, self.n_i, self.c_i)
    #NB: The geometry restraints module from cctbx is a convenient way to get geometry calculations
    #  angle_ideal and weight are required arguments for these functions
    #  I believe the angle_ideal allows you to get sigmas standard deviation

  def print_params(self,pdbid="pdbid"):
    #format from suitename input
    #"blank:1:chainID,resnum,i_code,altloc,resname,alpha,beta,gamma,delta,epsilon,zeta"
    return "%s:%s:%2s:%s:%s:%1s:%s:%.3f:%.3f:%.3f:%.3f:%.3f:%.3f:%.3f:%.3f:%.3f" % (
           " ",
           "1",
           self.three[1].parent().parent().id,
           self.three[1].resseq,
           self.three[1].icode,
           self.three[1].parent().altloc,
           self.three[1].resname,
           self.pperp1,
           self.d11,
           self.ang1,
           self.ang2,
           self.dhNN,
           self.pto11,
           self.pdown,
           self.dhPN,
           self.pperp2)

def find_atom(residue, atomname):
  for atom in residue.atoms():
    if atom.name == atomname:
      return atom
  return None

def distance(a,b):
  return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)**0.5

def perptersect(a1, a2, b1):
  #Finds the line from a1 to a2, drops a perpendicular to it from b1, and returns
  #  the point of intersection.
  A = [a2[0]-a1[0], a2[1]-a1[1], a2[2]-a1[2]]
  #Find the slope of line A in each direction, A is in vector notation
  t = (A[0]*(b1[0]-a1[0]) + A[1]*(b1[1]-a1[1]) + A[2]*(b1[2]-a1[2])) / ((A[0]**2)+(A[1]**2)+(A[2]**2))
  #Solve the parametric equations (dot of perpendiculars=0). . .
  b2 = [a1[0]+A[0]*t, a1[1]+A[1]*t, a1[2]+A[2]*t]
  # . . . and use the result to find the new point b2 on the line
  return b2

def are_linked(prevres, res):
  target_dist = 1.673 #1.593 + 0.02*4
  #from the geostd cif for A:
  #A  P    O5*   coval  1.593  0.020  1.593
  #the needed O3'-P distance should be the same
  #set the cutoff at geometry outlier distance of 4 sigma
  o = find_atom(prevres, " O3'")
  p = find_atom(res," P  ")
  if None in [o,p]: return False
  dist = distance(o.xyz,p.xyz)
  if dist > target_dist: return False
  return True

def are_glyco_bonded(n,c):
  #n and c should be 1.4 to 1.5A apart (generous)
  #this check is used to see if the atom search has found the correct n atom
  #I don't trust naming in non-standard or modified bases
  if None in [n,c]: return False 
  d = distance(c.xyz,n.xyz)
  if d < 1.4 or d > 1.5: return False
  return True

def three_linked(three):
  if are_linked(three[0],three[1]) and are_linked(three[1], three[2]):
    return True
  else:
    return False

def do_params(three):
  pass

def csv_header():
  #format from suitename input
  #"blank:1:chainID,resnum,i_code,altloc,resname,alpha,beta,gamma,delta,epsilon,zeta"
  headerlist = ["blank",
                "1",
                "chainID",
                "resnum",
                "i_code",
                "altloc",
                "resname",
                "Pperp1",
                "d1'1'",
                "ang1",
                "ang2",
                "dhNN",
                "pto1'1'",
                "Pdown",
                "dhPN",
                "Pperp2"]
  return ":".join(headerlist)

def run_mpgeo(pdbfilepath):
  mpgeo_command = "mmtbx.mp_geo rna_backbone=True pdb=%s" % pdbfilepath
  mpgeo_out = easy_run.fully_buffered(mpgeo_command)
  #if mpgeo_out.return_code != 0:
  #  raise RuntimeError("mmtbx.mp_geo crashed, check file formatting - dumping stderr:\n%s" %
  #    "\n".join(mpgeo_out.stderr_lines))
  return mpgeo_out.stdout_lines
  #"mmtbx.mp_geo rna_backbone=True pdb=$infile > $midfile",$arg_list_filler,$mpgeo_return_code

def store_mpgeo_lines(lines):
  #for formatting see mmtbx/validations/utils line~280
  #blank:1:chain:resseq:icode:altloc:resname:alpha,beta,gamma,delta,epsion,zeta
  # :1: A:  57: : :  G:-65.658:167.070:57.463:81.657:__?__:__?__
  # :1: A:  59: : :  U:__?__:-158.840:63.678:84.575:-148.820:-53.724
  # :1: A:  60: : :  C:-72.162:179.454:66.018:148.304:-97.069:-66.356
  for line in lines:
    x = line.split(":")
    chain = x[2].strip()
    resseq = x[3]
    i_code = x[4]
    altloc = x[5].strip()
    resname = x[6]
    alpha = x[7]
    beta = x[8]
    gamma = x[9]
    delta = x[10]
    epsilon = x[11]
    zeta = x[12].strip()
  reskey = ":".join([chain,resseq,icode,altloc])


#infilepath = sys.argv[1]
#filebase = os.path.basename(infilepath).split('.')[0]

def run(infilepath):
  filebase = os.path.basename(infilepath).split('.')[0]
  pdb_io = iotbx.pdb.input(infilepath)
  pdb_hierarchy = pdb_io.construct_hierarchy()

  result = []
  result.append(csv_header())
  for chain in pdb_hierarchy.chains():
    if not chain.is_na(): continue
    for conf in chain.conformers():
      if conf.altloc not in ['','A']: continue
      three = [None,None,None]
      for res in conf.residues():
        three.append(res)
        popped = three.pop(0)
        if None in three: continue
        if not three_linked(three): continue
        rna = rna_params(three)
        if rna.complete:
          result.append(rna.print_params(filebase))
  return result

if __name__ == "__main__":
#input validation
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
