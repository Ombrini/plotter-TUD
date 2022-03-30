import numpy as np
from scipy import io, integrate, linalg, signal
from scipy.sparse.linalg import eigs
from matplotlib import pyplot as plt
import daetools
from daetools.pyDAE import *
from pyUnits import *
from pyCore import *

#setting a fiew constants
e = 1.6e-19 
eps0 = 8.85e-12
epsR = 78.5 #dielectic permittivity of water
eps = eps0*epsR
kB = 1.38e-23
c0 = 1000 # initial conc of anions and cations in mol/m**3
psi0V = 0.025  #value of potential in Volt
Temp = 298
Na = 6.02e23 
D = 1e-9 #diff coeff in m**2/s

#setting times and speces to get clear results
LambdaD = ((eps*kB*Temp)/(2*(e**2)*Na*c0))**(1/2)
L = LambdaD*8 
RC = (L*LambdaD)/D
time = RC*8

C = eps/LambdaD
Q = C*psi0V
print(Q)


#solving the linearinzed solution 

class PoissonEq(daeModel):
    def __init__(self, Name, Parent = None, Description = ""):
        daeModel.__init__(self,Name,Parent,Description)

        self.x = daeDomain("x", self,m, "x axis domain")

        self.D = daeParameter("D", (m**2)/s, self, "Diffusion coefficinet")
        self.LambdaD = daeParameter("LambdaD", m, self, "debye lenght")
        self.psi0 = daeParameter("psi0", m/m , self, "potential at the wall")

        self.rho = daeVariable("rho", no_t, self)
        self.rho.DistributeOnDomain(self.x)
        self.rho.Description = "charge density"

        self.c = daeVariable("c", no_t, self)
        self.c.DistributeOnDomain(self.x)
        self.c.Description = "neutral density"

        self.psi = daeVariable("psi", no_t, self)
        self.psi.DistributeOnDomain(self.x)
        self.psi.Description = "potential"

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        eq = self.CreateEquation("chargeEq", "Charge density conservation")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = dt(self.rho(x)) - d(self.D()*(d(self.rho(x),self.x,eCFDM)+self.c(x)*d(self.psi(x),self.x,eCFDM)), self.x,eCFDM)

        eq = self.CreateEquation("neutralEq", "netrual density conservation")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = dt(self.c(x)) - d(self.D()*(d(self.c(x),self.x,eCFDM)+self.rho(x)*d(self.psi(x),self.x,eCFDM)), self.x,eCFDM)

        eq = self.CreateEquation("PoissonEq", "Poisson equation")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = d2(self.psi(x),self.x,eCFDM)+ ((1/self.LambdaD())**2)*self.rho(x)

        #Boundary Conditions


        # Potential BC
        eq = self.CreateEquation("DirichBC1", "0 potential at x=0")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.psi(x)

        eq = self.CreateEquation("DirichBC2", "psi0 potential at x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = self.psi(x) + self.psi0()


        #Flux BC
        eq = self.CreateEquation("fluxCondRho", "rho result from the 0 flux condition at x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = d(self.rho(x),self.x, eCFDM) + self.c(x)*d(self.psi(x),self.x, eCFDM)

        eq = self.CreateEquation("fluxCondC", "c result from the 0 flux condition at x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = d(self.c(x),self.x, eCFDM) + self.rho(x)*d(self.psi(x),self.x, eCFDM)

        #Bulk BC
        eq = self.CreateEquation("chargeinBulk0", "0 charge in bulk elyte condition")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.rho(x)

        eq = self.CreateEquation("entrualnBulk1", "1 neutral in bulk elyte condition")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.c(x) - 1

class simPoissonEq(daeSimulation):
    def __init__(self):
        daeSimulation.__init__(self)
        self.m = PoissonEq("name")

    def SetUpParametersAndDomains(self):
        
        self.m.x.CreateStructuredGrid(200,0,L)

        self.m.psi0.SetValue(psi0V*e/(kB*Temp))
        self.m.D.SetValue(D *((m**2)/s))
        self.m.LambdaD.SetValue(LambdaD * m)
        
    def SetUpVariables(self):
        for x in range(1,self.m.x.NumberOfPoints -1):
            self.m.rho.SetInitialCondition(x,0)
            self.m.c.SetInitialGuess(x,1)
            self.m.psi.SetInitialCondition(x,0)
            

log = daePythonStdOutLog()
daesolver = daeIDAS()
datareporter = daeTCPIPDataReporter()
simulation = simPoissonEq()


simulation.m.SetReportingOn(True)

simulation.ReportingInterval = time/50
simulation.TimeHorizon = time

datareporter.Connect("", "PoissonEq")

simulation.Initialize(daesolver,datareporter,log)
simulation.SolveInitial()
simulation.Run()
simulation.Finalize()