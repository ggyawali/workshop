import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties
title_font = {'fontname':'Times New Roman', 'size':'15', 'color':'black', 'weight':'normal','verticalalignment':'bottom'} # Bottom vertical alignment for more space
lable_font = {'fontname':'Times New Roman', 'size':'15'}
fontP = FontProperties()
fontP.set_size('small')
plt.ylim([-0.3,0.2])
plt.title("Lennard-Jones Potential")
plt.xlabel("r(A)")
plt.ylabel("V(kcals/mol)")
plt.legend(loc='upper right',fancybox=True,prop=fontP)
plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)

def energy(e, s, r):
            """Calculates the pairwise Lennard-Jones potential using given epsilon, sigma and radial distance """
            return 4.0*e*((s/r)**12-(s/r)**6)


def plot_lj(e,s):
	
	rmax = 10
	dr = 0.01
	N = int(rmax/dr)
	rlist = []
	ljlist = []
	for i in range(1,N):
		rlist.append(i*dr)
		ljlist.append(energy(e,s,i*dr))
	plt.plot(rlist, ljlist, label= "e = " + str(e) + " s = " + str(s), linewidth = 2)


 """"This is the part where you make the changes""""
plot_lj(0.24, 3.4)
plot_lj(0.1,3.4)
""" syntax: plot_lj(epsilon, sigma)"""

plt.show()
