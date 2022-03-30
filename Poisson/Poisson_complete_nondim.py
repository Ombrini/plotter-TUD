import numpy as np
from scipy import io, integrate, linalg, signal
from scipy.sparse.linalg import eigs
from matplotlib import pyplot as plt
import daetools
from daetools.pyDAE import *
from pyUnits import *
from pyCore import *

elec_pot_t = daeVariableType(name = "elec_pot_t", units = V, lowerBound = -1e20, upperBound = 1e20, initialGuess = 0, absTolerance = 1e-5)

#setting a fiew constants
e = 1.6e-19 
eps0 = 8.85e-12
epsR = 78.5 #dielectic permittivity of water
eps = eps0*epsR
kB = 1.38e-23
R = 8.3144598 #[J/(mol*K)]

 
Temp = 298
Na = 6.02e23 

psi0 = 1 
DA = 1.1e-9
DC = 0.9e-9
c0 = 1000 #bulk conc in mol/m3
cA0 = 1 #Aninon norm conc
cC0 = cA0 #Cation norm conc

#setting times and speces to get clear results
LambdaD = ((eps*kB*Temp)/(2*(e**2)*Na*c0))**(1/2)
L = LambdaD*8 
RC = (L*LambdaD)/DA
time = RC*15

#Number of grid points
Npoints = 200

#solving the linearinzed solution 

class PoissonEq(daeModel):
    def __init__(self, Name, Parent = None, Description = ""):
        daeModel.__init__(self,Name,Parent,Description)

        self.x = daeDomain("x", self,m, "x axis domain")

        self.DA = daeParameter("DA", (m**2)/s, self, "Diffusion coefficinet of anions")
        self.DC = daeParameter("DC", (m**2)/s, self, "Diffusion coefficinet of cations")
        self.LambdaD = daeParameter("LambdaD", m, self, "debye lenght")

        self.psi0 = daeParameter("psi0", V/V , self, "potential at the wall")
        self.cA0 = daeParameter("cA0", mol/mol, self, "initial anion concentration")
        self.cC0 = daeParameter("cC0", mol/mol, self, "initial cation concentration")

        self.cA = daeVariable("cA", no_t, self )
        self.cA.DistributeOnDomain(self.x)
        self.cA.Description = "Anion conc"

        self.cC = daeVariable("cC", no_t, self )
        self.cC.DistributeOnDomain(self.x)
        self.cC.Description = "Cation conc"

        self.psi = daeVariable("psi", no_t, self)
        self.psi.DistributeOnDomain(self.x)
        self.psi.Description = "potential"

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        ################## fluxes
        # -self.DA()*(d(self.cA(x),self.x,eCFDM) + self.cA(x)*d(self.psi(x),self.x,eCFDM))
        # -self.DC()*(d(self.cC(x),self.x,eCFDM) - self.cC(x)*d(self.psi(x),self.x,eCFDM))
        ###############

        eq = self.CreateEquation("AnionCons", "Anion Conservation eq")
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = dt(self.cA(x)) + self.DA()*d(-d(self.cA(x),self.x,eCFDM) + self.cA(x)*d(self.psi(x),self.x,eCFDM), self.x, eCFDM)

        eq = self.CreateEquation("CationCons", "Cation Conservation eq")
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = dt(self.cC(x)) + self.DC()*d(-d(self.cC(x),self.x,eCFDM) - self.cC(x)*d(self.psi(x),self.x,eCFDM), self.x, eCFDM)

        eq = self.CreateEquation("PoissonEq", "Poisson equation")
    
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = d2(self.psi(x),self.x,eCFDM) + 0.5*(1/(self.LambdaD()**2))*(self.cC(x)-self.cA(x))


        #Boundary Conditions for the ions

        eq = self.CreateEquation("ZeroBulkchargeA", "bulk anion con")

        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.cA(x) - self.cA0()
        
        eq = self.CreateEquation("ZeroBulkchargeC", "bulk cation conc")
    
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.cC(x) - self.cC0()


        eq = self.CreateEquation("fluxCond0Anions", "no anion flux x=L")

        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = (d(self.cA(x),self.x,eCFDM) + self.cA(x)*d(self.psi(x),self.x,eCFDM))

        eq = self.CreateEquation("fluxCond0Cations", "no cation flux x=L")

        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = (d(self.cC(x),self.x,eCFDM) - self.cC(x)*d(self.psi(x),self.x,eCFDM))
           #Boundary Conditions for psi

        eq = self.CreateEquation("DirichBC1", "0 potential at x=0")
 
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.psi(x)

        eq = self.CreateEquation("DirichBC2", "psi0 potential at x=L")
    
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = self.psi(x) + self.psi0()

        
class simPoissonEq(daeSimulation):
    def __init__(self):
        daeSimulation.__init__(self)
        self.m = PoissonEq("name")

    def SetUpParametersAndDomains(self):
        

        self.m.x.CreateStructuredGrid(Npoints,0,LambdaD*8)

        # pointsArray = np.zeros(1)
        # for i in range(Npoints):
        #     pointsArray = np.append(pointsArray, pointsArray[i]+(LambdaD*8-pointsArray[i])/100)

        # pointsArray = np.zeros(1)
        # for i in range(Npoints):
        #     pointsArray = np.append(pointsArray, np.log(i+2))
        # pointsArray = (pointsArray/np.amax(pointsArray))*Npoints
        # print(pointsArray)

        # new_grid = pointsArray
        # self.m.x.Points = new_grid

        self.m.psi0.SetValue(psi0)
        self.m.DA.SetValue(DA *((m**2)/s))
        self.m.DC.SetValue(DC *((m**2)/s))
        # self.m.e.SetValue(e * A*s)
        # self.m.eps.SetValue(eps * F/m)
        self.m.cA0.SetValue(cA0)
        self.m.cC0.SetValue(cC0)
        # self.m.Na.SetValue(Na * 1/(mol))
        # self.m.T.SetValue(Temp *K)
        # self.m.R.SetValue(R * J/(mol*K))
        # self.m.F.SetValue(e*Na *A*s/mol)
        self.m.LambdaD.SetValue(LambdaD *m)


        # self.m.LambdaD.SetValue(LambdaD * m)
        
    def SetUpVariables(self):
        for x in range(1,self.m.x.NumberOfPoints -1):
            self.m.cC.SetInitialCondition(x,1)
            self.m.cA.SetInitialCondition(x,1)
            # self.m.psi.SetInitialGuess(x,0)

log = daePythonStdOutLog()
daesolver = daeIDAS()
datareporter = daeTCPIPDataReporter()
simulation = simPoissonEq()


simulation.m.SetReportingOn(True)

simulation.ReportingInterval = time/500
simulation.TimeHorizon = time

datareporter.Connect("", "PoissonEq")

simulation.Initialize(daesolver,datareporter,log)
simulation.SolveInitial()
simulation.Run()
simulation.Finalize()