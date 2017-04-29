# This python script reads from the command line the filename (its root) and
# the largest exponent of 10. This defines the number of mesh points.
# This script calls then an executable from a c++ that solves
# an eigenvalue problem for an electron and two electrons in a harmonic oscillator well
# It makes in turn the various plots as pdf files and finally sets up the basis for
# a report and its pertinent latex file

import sys, os
from  matplotlib import pyplot as plt
import numpy as np

# Command line arguments using sys.argv[]
try:
    filename = sys.argv[1]
    dimension = int(sys.argv[2])
except:
    print ("Usage of this script", sys.argv[0], "infile", sys.argv[1], "Dimension", sys.argv[2], sys.exit(1))

# Define command line text string 
cmdline  = './main.x '+filename +' ' + str(dimension)
# Now run code, here c++ code  which has been compiled and linked
cmd = cmdline
failure = os.system(cmd)
if failure:
   print ('running project2 failed', sys.exit(1))

figfile_psi = filename+'psi.pdf'
figfile_psi_square = filename+'psi_square.pdf'
count = 0
colors = ['r','m','g','b']
maxdim = [50,10,5,5]
j = 0

for i in [0.01,0.5,1.0,5.0]:
    fileread = filename+str(count)
    data = np.loadtxt(fileread)
    r = data[:,0]
    psi = data[:,1] #r*np.exp(-E[count]*r**2/2)
    plt.hold(True)
    plt.plot(r,psi,colors[1]+':.',linewidth = 2.0,label = 'non-interacting')
    count+=1

    fileread2 = filename+str(count)
    data2 = np.loadtxt(fileread2)
    r2 = data2[:,0]
    psi2 = data2[:,1] #r*np.exp(-E[count]*r**2/2)
    plt.plot(r2,psi2,colors[2]+':.',linewidth = 2.0,label = 'interacting')
    #plt.yscale('log')
    plt.autoscale(enable=True, axis = 'y')
    plt.xlabel(r'$r$',fontsize=16)
    plt.ylabel(r'$rPsi$',fontsize=16)
    plt.xlim((0,maxdim[j]))
    plt.title('omega ='+str(i))
    plt.legend()
    plt.savefig(filename+'psi'+str(count-1))
    plt.clf()
    count+=1

    plt.hold(True)
    plt.plot(r,psi**2,colors[1]+':.',linewidth = 2.0,label = 'non-interacting')
   
    plt.plot(r2,psi2**2,colors[2]+':.',linewidth = 2.0,label = 'interacting')
    #plt.yscale('log')
    plt.autoscale(enable=True, axis = 'y')
    plt.xlabel(r'$r$',fontsize=16)
    plt.ylabel(r'$r^2Psi^2$',fontsize=16)
    plt.xlim((0,maxdim[j]))
    plt.title('omega ='+str(i))
    plt.legend()
    plt.savefig(filename+'psisquare'+str(count-2))
    plt.clf()
    j+=1




# Now prepare latex file, r in front avoids backslashes being treated
# as control chars in strings. What follows are plain  latex commands
preamb = r"""\documentclass[10pt,showpacs,preprintnumbers,footinbib,amsmath,amssymb,aps,prl,twocolumn,groupedaddress,superscriptaddress,showkeys]{revtex4-1}
\usepackage{graphicx}
\usepackage{dcolumn}
\usepackage{bm}
\usepackage[colorlinks=true,urlcolor=blue,citecolor=blue]{hyperref}
\usepackage{color}
\begin{document}
\title{Project 1}
\author{A.~N.~Author}
\affiliation{Department of Something, University of Somewhere, Outer Space}
\begin{abstract}
We present our Ferrari algorithm for solving linear equations. Our best algorithm runs as $4n$ FLOPS with $n$ the dimensionality of the matrix.
\end{abstract}
\maketitle

"""

figure = r"""\begin{figure}[hbtp]
\includegraphics[scale=0.4]{test1.pdf}
\caption{Exact and numerial solutions for $n=10$ mesh points.} 
\label{fig:n10points}
\end{figure}

"""


introduction = r"""\section{Introduction}

"""

theory = r"""\section{Theory, algorithms and methods}

"""

results = r"""\section{Results and discussions}

"""

conclusions = r"""\section{Conclusions}

"""

references = r"""\begin{thebibliography}{99}
\bibitem{miller2006} G.~A.~Miller, A.~K.~Opper, and E.~J.~Stephenson, Annu.~Rev.~Nucl.~Sci.~{\bf 56}, 253 (2006).
\end{thebibliography}

"""

# Dump to file:
filename = 'ReportProject2'
f = open(filename + '.tex', "w")
f.write(preamb)
f.write(introduction)
f.write(theory)
f.write(results)
f.write(figure)
f.write(conclusions)
f.write(references)
f.write("""\end{document}""")
f.close()

