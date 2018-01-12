#
#
#
#
#
from MMC import *
import copy
import math

pdb1 = pdbreader("argon.pdb")
b_curr = box(18.35, pdb1.atoms)
k = 0.0019872041    #Boltzmann constant in kcals/(mol.K)
T = 91.8
nsteps = 500
N = b_curr.natoms
nacc = 0

for i in range(nsteps):
    e_curr = b_curr.energy()
    b_trial = copy.deepcopy(b_curr)
    r1 = randint(0,N-1)
    dx = random() - 0.5
    dy = random() - 0.5
    dz = random() - 0.5
    b_trial.atoms[r1].translate(dx, dy, dz)
    e_trial = b_trial.energy()
    e_move = e_trial - e_curr
    if e_move <= 0:
        b_curr = copy.deepcopy(b_trial)
        nacc +=1
        print "Accepted by 1"
    else:
        boltz = math.exp(-e_move/(k * T))
        r2 = random()
        if boltz > r2:
            b_curr = copy.deepcopy(b_trial)
            nacc +=1
            print "Accepted by 2"
        else:
            print "Rejected"
    b_curr.writexyz("run1")
       





