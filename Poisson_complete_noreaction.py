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
epsR = 30 #dielectic permittivity of elyte
eps = eps0*epsR
kB = 1.38e-23
R = 8.3144598 #[J/(mol*K)]


# psi0V = 0.001  #value of potential in electron-Volt
Temp = 298
Na = 6.02e23 
Fara = e*Na

psiStern = 0.01 #stern layer drop in volt
xStern = 5e-10 #stern layer thickness




DA = 2.2e-10
DC = 2.9e-10
# DA = 0
cA0 = 1000 #Aninon conc
cC0 = cA0 #Cation conc

#setting times and speces to get clear results
LambdaD = ((eps*kB*Temp)/(2*(e**2)*Na*cA0))**(1/2)
L = LambdaD*8 
RC = (L*LambdaD)/((DA+DC)*0.5)
time = RC*15

psiDL = (2*R*Temp/Fara)*np.arcsinh(np.sqrt(eps/(8*R*Temp*cA0))/xStern*psiStern) #works decently for similar D, ok for liquid elyte

#Number of grid points
Npoints = 500

#solving the linearinzed solution 

class PoissonEq(daeModel):
    def __init__(self, Name, Parent = None, Description = ""):
        daeModel.__init__(self,Name,Parent,Description)

        self.x = daeDomain("x", self,m, "x axis domain")

        self.DA = daeParameter("DA", (m**2)/s, self, "Diffusion coefficinet of anions")
        self.DC = daeParameter("DC", (m**2)/s, self, "Diffusion coefficinet of cations")

        # self.psi0 = daeParameter("psi0", V , self, "potential at the wall")
        self.psiStern = daeParameter("psiStern", V, self, "potential drop in the stern layer")
        self.xStern = daeParameter("xStern", m , self, "thickness of the stern layer")
     
        self.eps = daeParameter("eps", F/m, self, "dielectric constant")
        self.cA0 = daeParameter("cA0", mol*(m**(-3)), self, "initial anion concentration")
        self.cC0 = daeParameter("cC0", mol*(m**(-3)), self, "initial cation concentration")
     
        self.T = daeParameter("T", K, self, "temperature")
        self.R = daeParameter("R", J/(mol*K), self, "gas constant")
        self.F = daeParameter("F", (A*s)/mol, self, "Faraday constant")

        self.cA = daeVariable("cA", molar_concentration_t, self )
        self.cA.DistributeOnDomain(self.x)
        self.cA.Description = "Anion conc"

        self.cC = daeVariable("cC", molar_concentration_t, self )
        self.cC.DistributeOnDomain(self.x)
        self.cC.Description = "Cation conc"

        self.psi = daeVariable("psi", elec_pot_t, self)
        self.psi.DistributeOnDomain(self.x)
        self.psi.Description = "potential"

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        ################## fluxes
        # -self.DA()*d(self.cA(x),self.x,eCFDM) + (self.DA()/(self.R()*self.T()))*self.F()*self.cA(x)*d(self.psi(x),self.x,eCFDM)
        # -self.DC()*d(self.cC(x),self.x,eCFDM) - (self.DC()/(self.R()*self.T()))*self.F()*self.cC(x)*d(self.psi(x),self.x,eCFDM)
        ###############

        eq = self.CreateEquation("AnionCons", "Anion Conservation eq")
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = dt(self.cA(x)) + d(-self.DA()*d(self.cA(x),self.x,eCFDM) + (self.DA()/(self.R()*self.T()))*self.F()*self.cA(x)*d(self.psi(x),self.x,eCFDM), self.x, eCFDM)

        eq = self.CreateEquation("CationCons", "Cation Conservation eq")
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = dt(self.cC(x)) + d(-self.DC()*d(self.cC(x),self.x,eCFDM) - (self.DC()/(self.R()*self.T()))*self.F()*self.cC(x)*d(self.psi(x),self.x,eCFDM), self.x, eCFDM)

        eq = self.CreateEquation("PoissonEq", "Poisson equation")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eOpenOpen)
        eq.Residual = self.eps()*d2(self.psi(x),self.x,eCFDM) + self.F()*(self.cC(x)-self.cA(x))
         
         
         #Boundary Conditions for psi

        eq = self.CreateEquation("DirichBC1", "0 potential at x=0")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.psi(x)

        #BC for imposed stern layer potentia drop
        eq = self.CreateEquation("NewmannBC", "derivative at x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual =d(self.psi(x),self.x,eCFDM) - self.psiStern()/self.xStern()


        #Boundary Conditions for the ions

        eq = self.CreateEquation("ZeroBulkchargeA", "0 charge in bulk elyte condition anions")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.cA(x) - self.cA0()
        
        eq = self.CreateEquation("ZeroBulkchargeC", "0 charge in bulk elyte condition cations")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.cC(x) - self.cC0()


        eq = self.CreateEquation("fluxCond0Anions", "no anion flux x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = -d(self.cA(x),self.x,eCFDM) + (1/(self.R()*self.T()))*self.F()*self.cA(x)*d(self.psi(x),self.x,eCFDM)

        eq = self.CreateEquation("fluxCond0Cations", "no cation flux x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = -d(self.cC(x),self.x,eCFDM) - (1/(self.R()*self.T()))*self.F()*self.cC(x)*d(self.psi(x),self.x,eCFDM)

  


        
class simPoissonEq(daeSimulation):
    def __init__(self):
        daeSimulation.__init__(self)
        self.m = PoissonEq("name")

    def SetUpParametersAndDomains(self):
        
        DomainLenght = LambdaD*10
        DL = DomainLenght
        self.m.x.CreateStructuredGrid(Npoints,0,DL)

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

        # self.m.psi0.SetValue(psi0V * V)
        self.m.psiStern.SetValue(psiStern * V)
        self.m.xStern.SetValue(xStern * m)

        self.m.DA.SetValue(DA *((m**2)/s))
        self.m.DC.SetValue(DC *((m**2)/s))
        # self.m.e.SetValue(e * A*s)
        self.m.eps.SetValue(eps * F/m)
        self.m.cA0.SetValue(cA0 * mol/(m**3))
        self.m.cC0.SetValue(cC0 * mol/(m**3))
        # self.m.Na.SetValue(Na * 1/(mol))
        self.m.T.SetValue(Temp *K)
        self.m.R.SetValue(R * J/(mol*K))
        self.m.F.SetValue(e*Na *A*s/mol)

        


        # self.m.LambdaD.SetValue(LambdaD * m)
        
    def SetUpVariables(self):
        DomainLenght = LambdaD*8
        DL = DomainLenght
        for x in range(1,self.m.x.NumberOfPoints -1):
            self.m.cC.SetInitialCondition(x,cC0 * mol*(m**(-3)))
            self.m.cA.SetInitialCondition(x,cA0 * mol*(m**(-3)))
            self.m.psi.SetInitialGuess(x,0*V)

log = daePythonStdOutLog()
daesolver = daeIDAS()
datareporter = daeTCPIPDataReporter()
simulation = simPoissonEq()


simulation.m.SetReportingOn(True)

simulation.ReportingInterval = RC
simulation.TimeHorizon = RC*300

datareporter.Connect("", "PoissonEq")

simulation.Initialize(daesolver,datareporter,log)
simulation.SolveInitial()
simulation.Run()
simulation.Finalize()

print("analytical calculation for psiDL is: " , psiDL)