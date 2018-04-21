#! /usr/bin/environment python
import scipy
from scipy import interpolate
import numpy as np
from numpy import random
import math
import re
def main():
   
    DecayRate = 5*390000.0

    peakFile = open("10Mpeak.dat","r")
    totFile = open("10Mtot.dat","r")

    x1var = list()
    x2var = list()
    peakvar = list()
    totvar = list()

    for line in peakFile:
        x1var.append(float(line.split()[0]))
        peakvar.append(float(line.split()[1]))
    peakFile.close()

    for line in totFile:
        x2var.append(float(line.split()[0]))
        totvar.append(float(line.split()[1]))
    totFile.close()

    def GetPeakEfficiency(InputEnergy):
        cs = interpolate.CubicSpline(x1var,peakvar)
        return 0.950*cs(InputEnergy)

    def GetTotalEfficiency(InputEnergy):
        ctot = interpolate.CubicSpline(x2var,totvar)
        return 0.995*ctot(InputEnergy)

     
    measuredExp     = np.zeros(shape=(100,100))
    measuredExpErr  = np.zeros(shape=(100,100))
    branchingRat    = np.zeros(shape=(100,100))
    branchingRatErr = np.zeros(shape=(100,100))

    with open("27Al.dat","r") as inputfile:
        content = inputfile.readlines()
    content   = [x.strip() for x in content]
    

    """ determine what line everything is on """
    branchingline = 0
    measuredline  = 0


    """ First, read header"""
    numLevels = int(content[1].split()[0])+1

    for i in range(0,len(content)):
        if re.search("Energy-Levels",content[i]):
           energies  = np.array(map(lambda v: float(v[0]),map(lambda b: b.split() ,content[i+1:i+1+numLevels])))
           f         = np.array(map(lambda v: float(v[1]),map(lambda b: b.split() ,content[i+1:i+1+numLevels])))

    for i in range(0,len(content)):
        if re.search("Measured",content[i]):
            measuredline = i
            temp = np.array(map(lambda v: (int(v[0]),int(v[1]),float(v[2]),float(v[3])) , map(lambda b: b.split(),content[i+1:])))
            print temp
            for j in temp:
                measuredExp[int(j[0]),int(j[1])]    = j[2]
                measuredExpErr[int(j[0]),int(j[1])] = j[3]

    for i in range(0,len(content)):
        if re.search("B-Values",content[i]):
            branchingline = i
            temp = np.array(map(lambda v: (int(v[0]),int(v[1]),float(v[2]),float(v[3])) , map(lambda b: b.split(),content[i+1:measuredline-1])))
            print temp
            for j in temp:
                branchingRat[int(j[0]),int(j[1])]    = j[2]
                branchingRatErr[int(j[0]),int(j[1])] = j[3]
    

    x = np.mat([[0]*numLevels]*numLevels,dtype=float)

    measured = np.mat([[0]*numLevels]*numLevels,dtype=float)

    """
    measuredExp[1][0] = 55515.0
    measuredExp[3][0] = 4424.0
    measuredExp[6][3] = 5633.0
    measuredExp[6][1] = 3977.0
    measuredExp[6][0] = 3383.0
    measuredExp[13][6] = 6996.7
    measuredExp[13][4] = 11848.0
    measuredExp[13][1] = 2478.0

    measuredExpErr[1][0] = 320.0
    measuredExpErr[3][0] = 102.0
    measuredExpErr[6][3] = 113.0
    measuredExpErr[6][1] = 104.0
    measuredExpErr[6][0] = 84.0
    measuredExpErr[13][6] = 109.0
    measuredExpErr[13][4] = 138.0
    measuredExpErr[13][1] = 86.0
    """
    
    index = 0
    for j in range(0,numLevels):
        for i in range(0,numLevels):
            if measuredExpErr[j][i] != 0.0:
                index += 1.0
    table=[list() for element in range(0,int(index))]

    for uncIter in range(0,100):
        print uncIter
        """define Branching ratio matrix xji"""
        for j in range(1,numLevels):
            for i in range(0,numLevels):
                x[j,i] = branchingRat[j][i] + random.normal(0.0,branchingRatErr[j][i])

        """define measured (Sji)"""
        for j in range(1,numLevels):
            for i in range(0,numLevels):
                if measuredExpErr[j][i] != 0.0:
                    measured[j,i] = measuredExp[j][i]+random.normal(0.0,measuredExpErr[j][i])

        """Normalize Branching ratio matrix xji"""
        for j in range(1,numLevels):
            rowSum = 0.0
            for i in range(0,numLevels):
                rowSum += x[j,i]
            for i in range(0,numLevels):
                temp = 0.0
                temp += x[j,i]
                x[j,i] =  float(temp/rowSum)
            
        """define the "c" matrix"""
        c = np.mat([[0]*numLevels]*numLevels,dtype=float)
        for j in range(0,numLevels):
            for i in range(0,numLevels):
                c[j,i] = x[j,i]/(1.0)
                
        """define the "e_p" matrix"""
        e_p = np.mat([[0]*numLevels]*numLevels,dtype=float)
        for i in range(0,numLevels-1):
            for j in range(i+1,numLevels):
                e_p[j,i] = GetPeakEfficiency( math.fabs(energies[j]-energies[i]) )

        """define the "e_t" matrix"""
        e_t = np.mat([[0]*numLevels]*numLevels,dtype=float)
        for i in range(0,numLevels-1):
            for j in range(i+1,numLevels):
                e_t[j,i] = GetTotalEfficiency( math.fabs(energies[j]-energies[i]) )

        for iteration in range(0,10):

            """define the "a" matrix"""
            a = np.mat([[0]*numLevels]*numLevels,dtype=float)
            for i in range(0,numLevels-1):
                for j in range(i+1,numLevels):
                    a[j,i] = c[j,i]*e_p[j,i]

            """define the "e" matrix"""
            e = np.mat([[0]*numLevels]*numLevels,dtype=float)
            for i in range(0,numLevels-1):
                for j in range(i+1,numLevels):
                    e[j,i] = c[j,i]*e_t[j,i]

            """define the "b" matrix"""
            b = np.mat([[0]*numLevels]*numLevels,dtype=float)
            for i in range(0,numLevels-1):
                for j in range(i+1,numLevels):
                    b[j,i] = x[j,i] - e[j,i]

            """define the "A" matrix"""
            A = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            for i in range(1,numLevels):
                A += a**int(i)

            """define the "A0" matrix"""
            A0 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            for j in range(0,numLevels):
                for i in range(0,numLevels):
                    A0[j,i] = a[j,i]

            """define the "A1" matrix"""
            A1 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            A1 = a*a

            """define the "E" matrix"""
            E = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            for i in range(0,numLevels):
                E[i,i] = 1.0

            """define the "B" matrix"""
            B = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            B = E.copy()
            for i in range(1,numLevels):
                B += b**int(i)

            """define the "B0" matrix"""
            B0 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            B0 = E.copy()
            for i in range(1,numLevels):
                B0 += x**int(i)

            """define the "B1" matrix"""
            B1 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            B1 = E.copy()
            for k in range(1,numLevels):
                for l in range(0,(k-1)+1):
                    B1 -= (x**l)*e*(x**(k-l-1))
            
            """define the "N" matrix"""
            N = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            for i in range(0,numLevels):
                    N[i,i]=(f*B)[0,i]

            """define the "N0" matrix"""
            N0 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            for i in range(0,numLevels):
                    N0[i,i]=(f*B0)[0,i]

            """define the "N1" matrix"""
            N1 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            for i in range(0,numLevels):
                    N1[i,i]=(f*B1)[0,i]

            """define the "M" matrix"""
            M = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            for i in range(0,numLevels):
                    M[i,i]=B[i,0]
            
            """define the "M0" matrix"""
            M0 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            M0 = E.copy()

            """define the "M1" matrix"""
            M1 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            for i in range(0,numLevels):
                    M1[i,i]=B1[i,0]

            """define the "D0" matrix"""
            D0 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)

            """define the "D1" matrix"""
            D1 = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            D1 = N1*A0+N0*A1+N0*A0*M1

            """define the "D" matrix"""
            D = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            D = N0*A0*(M-M0)+N0*(A-A0)*M+(N-N0)*A*M
            
            S_nc = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            S_nc = DecayRate*(N0*A0)

            S_c = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            S_c = DecayRate*N*A*M

            I  = np.mat([ [0]*numLevels ]*numLevels,dtype=float)
            for j in range(0,numLevels):
                for i in range(0,numLevels):
                    I[j,i] = DecayRate*(N0*c)[j,i]
            #print "{0} Order Efficiencies".format(iteration)

            index = -1 
            for j in range(0,numLevels):
                for i in range(0,numLevels):
                    if measuredExpErr[j][i] != False:
                        index += 1.0
                        Ej=energies[j]
                        Ei=energies[i]
                        if iteration==0:
                            e_p[j,i] = measured[j,i]/I[j,i]
                        elif iteration==9:
                            table[int(index)].append(e_p.getA()[j][i])
                        else:
                            e_p[j,i] = (measured[j,i]/I[j,i]) - (D[j,i]/(N0*c)[j,i])
    index = -1 
    for j in range(0,numLevels):
        for i in range(0,numLevels):
            if measuredExpErr[j][i] != False:
                index += 1
                Ej=energies[j]
                Ei=energies[i]
                print "{0:1.4f}\t\t{1:1.5f}\t\t{2:1.5f}\t\t{3:1.5f}".format(Ej-Ei,np.percentile(table[index],16), np.percentile(table[index],50),np.percentile(table[index],84))
    """
    print "{0}\t\t\t\t\t{1}\t{2}".format("Energy","--no summing effects--"," --with coincidence summing--")
    for j in range(0,numLevels):
        for i in range(0,j):
            if S_c.getA()[j][i] != 0.0:
                Ej=energies[j]
                Ei=energies[i]
                no_summing = S_nc.getA()[j][i]
                with_summing = S_c.getA()[j][i]
                print "{0:6.1f} -> {1:6.1f} ({4:4.1f} keV)\t\t{2:10}\t\t{3:10}\t\t{5:2.2f}".format(Ej,Ei,no_summing,with_summing,(Ej-Ei),100.0*(no_summing/with_summing))
    for j in range(0,numLevels):
        for i in range(0,j):
            if S_c.getA()[j][i] != 0.0:
                Ej=energies.getA1()[j]
                Ei=energies.getA1()[i]
                print "{0:6.1f} -> {1:6.1f} ({2:4.1f} keV)\t\t{3:10}\t\t{4:10}".format(Ej,Ei,(Ej-Ei), e_p.getA()[j][i], S_c.getA()[j][i]/I.getA()[j][i])

    for j in range(0,numLevels):
        for i in range(0,j):
            if S_c.getA()[j][i] != 0.0:
                Ej=energies.getA1()[j]
                Ei=energies.getA1()[i]
                print "{0:6.1f} -> {1:6.1f} ({2:4.1f} keV)\t\t{3:10}".format(Ej,Ei,(Ej-Ei), I.getA()[j][i])
    """
if __name__ == '__main__':
    main()
