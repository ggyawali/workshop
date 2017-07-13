# # # # # # # #
# Molecular Monte Carlo Library
# Steven W. Rick and Gaurav Gyawali
# Email: ggyawali@uno.edu
# # # # # # # #

'''This library provides an environment to run molecular mechanical calculations. It's been especially written to run Monte Carlo minimization for
Summer Worshop on Computational Sciences'''

from math import sqrt, exp
from random import *
from numpy import linalg as la
import numpy as np
     
class atom:
      """ Class which builds an atom object that stores information about an atom"""
      count = 0

      def __init__(self, type, x, y, z):
            """Attributes:
            id = index number
            type = atom type
            xyz = coordinates of the form (x,y,z)
            x = x coordinate
            y = y coordinate
            z = z coordinate
            """
      
            atom.count = atom.count + 1
            self.id = atom.count
            self.type = type
            self.xyz = np.array((x,y,z))
            self.x   = x
            self.y   = y
            self.z   = z
            
      def info(self):
            """ Prints out the id, type as well as coordinates of an atom"""
            print self.id, self.type, self.xyz.get_xyz

      def translate ( self, dx, dy, dz):
            """Translates the coordinates from (x,y,z) to (x + dx, y+ dy, z + dy) by using the given translation
            vector (dx, dy, dz)"""
            self.x = self.x + dx
            self.y = self.y + dy
            self.z = self.z + dz
            self.xyz = np.array((self.x, self.y, self.z))

class box:
      """Class that defines a cubical box of given length consisting of the atom objects"""
      count = 0

      def __init__(self, length, alist, epsilon = 0.2385, sigma = 3.4):
            """Attributes:
             id = box id
             atoms = list of atoms
             natoms = # of atoms
             length = lenght of the box
             epsilon = epsilon parameter of Lennard-Jones potential
             sigma = sigma parameter of Lennard-Jones potential """
      
            box.count = box.count + 1
            self.id    = box.count
            self.atoms = alist
            self.natoms = len(alist)
            self.length = length
            self.epsilon = epsilon
            self.sigma = sigma
            
      def pair_energy(self,e, s, r):
            """Calculates the pairwise Lennard-Jones potential using given epsilon, sigma and radial distance """
            return 4.0*e*((s/r)**12-(s/r)**6)

      def pbc_correction(self, r):
            """Finds the nearest periodic boundary neighbor of given coordinate"""
            return (  self.length * int(round(r[0]/self.length)), self.length * int(round(r[1]/self.length)), self.length * int(round(r[2]/self.length)) )

      def energy(self):
            """Calculates the total energy of the box by summing up all the pairwise potential energies"""
            sum_energy = 0.0
            for i in range(0,self.natoms-1):
                for j in range(i+1,self.natoms):
                    rij = (self.atoms[i].xyz - self.atoms[j].xyz)
                    rij = rij - self.pbc_correction(rij)
                    mag_rij = la.norm(rij)
                    sum_energy = sum_energy + self.pair_energy(self.epsilon, self.sigma, mag_rij)                    
            return sum_energy

      def writexyz(self,fname):
            """Writes out the current atomic coordinates to the trajectory file in xyz format for visualization purposes """
            xyzfile = open(fname + ".xyz","a+")
            xyzfile.write(str(self.natoms) + "\n\n")
            for a in self.atoms:
                  	cxyz = a.xyz - np.array(self.pbc_correction(a.xyz))
			xyzfile.write(str(a.type) +  "\t" + str(cxyz[0]) + "\t" + str(cxyz[1]) + "\t" + str(cxyz[2]) + "\n")
            xyzfile.close()                 

      def writepdb(self,fname):
            """Writes out the "pdb" file for the given box"""
            pdbfile = open(fname + ".pdb", "w")
            for a in self.atoms:
                  pdbfile.write(str(a.type) + "\t" + str(a.x) + "\t" + str(a.y) + "\t" + str(a.z) + "\n")
            pdbfile.close()
                  

class pdbreader:
      
      """Class that reads in the atomic coordinates in pdb format
      Attributes:
      atoms: list of all the atoms read"""

      count = 0
      
      def __init__(self, pdbfile):
            self.atoms = []
            for line in open(pdbfile, "r"):
                  list = line.split()
                  if len(list)>0:
                        self.atoms.append(atom(list[0], float(list[1]), float(list[2]), float(list[3])))                        

class grcalc():
      """Class that calculates the self correlation function of a given box g(r), normalizes it and saves
      it to an output file"
      Attributes:
      r_max=  Maximum distance for the calculation of g(r)
      dr=     bin size
      r=      list of radii
      gr=     list of g(r)
      natoms= # of atoms
      nmax  = # of bins
      n_steps = # of steps for the calculation of g(r)
      length  = length of the box"""

      def __init__(self, natoms, r_max, dr):
            self.r_max = r_max
            self.dr    = dr
            self.r = []
            self.gr = []
            self.natoms = natoms
            self.n_max = int(r_max/dr)
            self.n_steps = 0
            self.length = 0
            for i in range(self.n_max):
                  self.r.append(i*dr)
                  self.gr.append(0)

      def pbc_correction(self, r):
            return (  self.length * int(round(r[0]/self.length)), self.length * int(round(r[1]/self.length)),  self.length * int(round(r[2]/self.length)) ) 

      def calculate(self, b):
            """Calculate the g(r) for a given box"""
            self.n_steps = self.n_steps + 1
            self.length  = b.length
            self.natoms = b.natoms
            for i in range(0,self.natoms-1):
                for j in range(i+1,self.natoms):
                    rij = (b.atoms[i].xyz - b.atoms[j].xyz)
                    rij = rij - self.pbc_correction(rij)
                    mag_rij = la.norm(rij)
                    bin_no  = int(round(mag_rij/self.dr))
                    if bin_no <= self.n_max:
                      self.gr[bin_no] = self.gr[bin_no] + 1

      def normalize(self):
            """Normalizes the g(r)"""
            dens  = float(self.natoms)/(self.length**3)
            for i in range(1,self.n_max):
                  volel = 4 * 3.1415 * i**2 * self.dr**3
                  norm  = 1/(dens * float(self.natoms-1)*self.n_steps*volel*0.5)
                  self.gr[i] = float(self.gr[i]) * norm

      def write(self,grfile):
            """Writes out the g(r) to an output file"""
            grfile = open(grfile + ".gr","w")
                      
            for i in range(1,self.n_max):
                  grfile.write(str(self.r[i]) + "\t" + str(self.gr[i]) + "\n")

            grfile.close()
                      
                      
                  
            
      
