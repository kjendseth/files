# =============================================================================
# 
# Script that converts an Amber formatted topology file and a coordinate file to a Orca force field parameter file and a Orca-readable xyz file.
# Written by Åsmund Røhr Kjendseth (2019)
# 
# Input:  
# parm7 topology file named: prmtop  
# coordinate file named: inpcrd
# These files are written by the AmberTools program tleap using the command: saveamberparm mol prmtop inpcrd
#
# More input/output options will be included.
#
# =============================================================================

import parmed as pmd

# Read Amber topology and coordinate files
parm = pmd.load_file('prmtop', xyz='inpcrd')


# Non bonded LJ 1-4 scaling factor
scnb = 2.0


# Colomb conversion factors (not used below) but may be the cause of slight differences in ee-terms.
CC_AMBER = 332.0522173 
CC_CHARMM = 332.054
CC_NAMD = 332.0636



# Write Amber formatted frcmod file (force constants). This file can be useful but is not used further below.
pmd.tools.writeFrcmod(parm, 'amber.frcmod').execute()





# List of elements according to atomic number
elements = [None,
         "H", "He",
         "Li", "Be",
         "B", "C", "N", "O", "F", "Ne",
         "Na", "Mg",
         "Al", "Si", "P", "S", "Cl", "Ar",
         "K", "Ca",
         "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
         "Ga", "Ge", "As", "Se", "Br", "Kr",
         "Rb", "Sr",
         "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
         "In", "Sn", "Sb", "Te", "I", "Xe",
         "Cs", "Ba",
         "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
         "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
         "Tl", "Pb", "Bi", "Po", "At", "Rn",
         "Fr", "Ra",
         "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
         "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Uub"]


# Add atom numbers
for i in range(0, len(parm.atoms)):
    parm.atoms[i].number = i
    

# Write Orca formatted force field file.
ff = open("ORCAFF.prms", "w")

# Write atom numbers, type, charge and LJ parameters. Note the scaling of 1-4 parameters are performed here and that the LJ parameteres are chosen to fit the charmm energy expression (sigma to rmin and epsilon * -1)
ff.write("$atoms\n")
ff.write("{0:1d} {1:1d} {2:1d}\n".format(len(parm.atoms), 1, 4))
for i in range(0, len(parm.atoms)):
    ff.write("{0:6d}   {1:2s}   {2:10.6f}   {3:10.6f}   {4:10.6f}   {5:10.6f}   {6:10.6f}\n".format(i + 1, elements[parm.atoms[i].atomic_number], parm.atoms[i].charge, -1 * parm.atoms[i].epsilon, parm.atoms[i].rmin *2, -1 * parm.atoms[i].epsilon / scnb , parm.atoms[i].rmin * 2))

# Amber and Charmm format identical
ff.write("$bonds\n")
ff.write("{0:1d} {1:1d} {1:1d} \n".format(len(parm.bonds), 2))
for i in range(0, len(parm.bonds)):
    ff.write("{0:6d} {1:10d} {2:10.6f} {3:10.6f}\n".format(parm.bonds[i].atom1.number + 1, parm.bonds[i].atom2.number + 1, parm.bonds[i].type.req, parm.bonds[i].type.k))

# Amber and Charmm format identical
ff.write("$angles\n")
ff.write("{0:1d} {1:1d} {2:1d}\n".format(len(parm.angles), 3, 2))
for i in range(0, len(parm.angles)):
    ff.write("{0:6d}   {1:6d}   {2:6d}   {3:10.6f}   {4:10.6f}\n".format(parm.angles[i].atom1.number + 1, parm.angles[i].atom2.number + 1, parm.angles[i].atom3.number + 1, parm.angles[i].type.theteq, parm.angles[i].type.k))

# Amber and Charmm format identical. Note that Amber impropers are included here.
ff.write("$dihedrals\n")
ff.write("{0:1d} {1:1d} {2:1d}\n".format(len(parm.dihedrals), 4, 3))
for i in range(0, len(parm.dihedrals)):
    ff.write("{0:6d}   {1:6d}   {2:6d}   {3:6d}   {4:10.6f}   {5:10.6f}   {6:6d}\n".format(parm.dihedrals[i].atom1.number + 1, parm.dihedrals[i].atom2.number + 1, parm.dihedrals[i].atom3.number + 1, parm.dihedrals[i].atom4.number + 1, parm.dihedrals[i].type.phase, parm.dihedrals[i].type.phi_k, parm.dihedrals[i].type.per))

ff.close()

# Write a formatted xyz file for ORCA
xyz = open("MM.xyz", "w")
xyz.write("{0:1d}\n".format(len(parm.atoms)))
xyz.write("Coordinates for Orca MM job\n")
for i in range(0, len(parm.atoms)):
    xyz.write("{0:2s}  {1:10.6f}  {2:10.6f}  {3:10.6f}\n".format(elements[parm.atoms[i].atomic_number], parm.atoms[i].xx, parm.atoms[i].xy, parm.atoms[i].xz))
xyz.close()

