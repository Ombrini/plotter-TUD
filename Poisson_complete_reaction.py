import numpy as np
from scipy import io, integrate, linalg, signal
from scipy.sparse.linalg import eigs
from matplotlib import pyplot as plt
import daetools
from daetools.pyDAE import *
import daetools.pyDAE as dae
from pyUnits import *
from pyCore import *

dae.daeGetConfig().SetString("daetools.IDAS.MaxNumItersIC","1000")
dae.daeGetConfig().SetString("daetools.IDAS.MaxNumSteps","10000")
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

psiStern = -0.01 #stern layer drop in volt
xStern = 5e-10 #stern layer thickness
potential1C = 0.3 #V

psiStern = -potential1C

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

jO = (4*DC*cC0/30e-4)*100 #J oxidation in mol/m2*s
jR = jO/cA0 # k reduction in m/s

k0 = 0.01 # exchange current density [A/m^2]
k0_mole = k0/Fara #exchange flux density in mol/m^2*s
current = 20 # normal 1C current [A/m^2]
current_mol = current/Fara #flux in mol/m^2*s

jO = k0_mole
jR = jO

psiDL = (2*R*Temp/Fara)*np.arcsinh(np.sqrt(eps/(8*R*Temp*cA0))/xStern*psiStern) #works decently for similar D, ok for liquid elyte

#Number of grid points
Npoints = 5000
DomainLenght = LambdaD*20
DL = DomainLenght

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


        self.jR = daeParameter("jR", mol/((m**2)*s), self, "reduction rate")
        self.jO = daeParameter("jO", mol/((m**2)*s), self, "oxidation flux")

        self.cA = daeVariable("cA", molar_concentration_t, self )
        self.cA.DistributeOnDomain(self.x)
        self.cA.Description = "Anion conc"

        self.cC = daeVariable("cC", molar_concentration_t, self )
        self.cC.DistributeOnDomain(self.x)
        self.cC.Description = "Cation conc"

        self.psi = daeVariable("psi", elec_pot_t, self)
        self.psi.DistributeOnDomain(self.x)
        self.psi.Description = "potential"

    def FluxC(self):
        JC = -self.DC()*d(self.cC(x),self.x,eCFDM) # - (self.DC()/(self.R()*self.T()))*self.F()*self.cC(x)*d(self.psi(x),self.x,eCFDM)
        return JC

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        ################## fluxes
        # -self.DA()*d(self.cA(x),self.x,eCFDM) + (self.DA()/(self.R()*self.T()))*self.F()*self.cA(x)*d(self.psi(x),self.x,eCFDM)
        # -self.DC()*d(self.cC(x),self.x,eCFDM) - (self.DC()/(self.R()*self.T()))*self.F()*self.cC(x)*d(self.psi(x),self.x,eCFDM)
        ###############
        # JA = -self.DC()*d(self.cC(x),self.x,eCFDM) - (self.DC()/(self.R()*self.T()))*self.F()*self.cC(x)*d(self.psi(x),self.x,eCFDM)

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

        #####################################################
        #Boundary Conditions for psi
        ######################################################

        eq = self.CreateEquation("DirichBC1", "0 potential at x=0")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.psi(x)

        # eq = self.CreateEquation("NewmannBC1", "derivative at x=0")
        # # eq.CheckUnitsConsistency = False
        # x = eq.DistributeOnDomain(self.x, eLowerBound)
        # eq.Residual =d(self.psi(x),self.x,eCFDM) - self.psiStern()/self.xStern()

        #BC for imposed stern layer potentia drop
        eq = self.CreateEquation("NewmannBC2", "derivative at x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual =d(self.psi(x),self.x,eCFDM) - self.psiStern()/self.xStern()


        #########################################
        #Boundary Conditions for the ions
        ##########################################


        # Left side BC

        eq = self.CreateEquation("ZeroBulkchargeA", "0 charge in bulk elyte condition anions")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.cA(x) - self.cA0()
        
        eq = self.CreateEquation("ZeroBulkchargeC", "0 charge in bulk elyte condition cations")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.cC(x) - self.cC0()

        # eq = self.CreateEquation("fluxCond0Anions", "no cation flux x=0")
        # # eq.CheckUnitsConsistency = False
        # x = eq.DistributeOnDomain(self.x, eLowerBound)
        # eq.Residual = d(self.cC(x),self.x,eCFDM) + (1/(self.R()*self.T()))*self.F()*self.cC(x)*d(self.psi(x),self.x,eCFDM)

        # self.IF(Time() < Constant(RC*100*s))
        # eq = self.CreateEquation("fluxCond0Cations", "no anion flux x=0")
        # # eq.CheckUnitsConsistency = False
        # x = eq.DistributeOnDomain(self.x, eLowerBound)
        # eq.Residual = d(self.cA(x),self.x,eCFDM) - (1/(self.R()*self.T()))*self.F()*self.cA(x)*d(self.psi(x),self.x,eCFDM)

        # self.ELSE()
        # eq = self.CreateEquation("fluxReduction", "anion flux x=0")
        # # eq.CheckUnitsConsistency = False
        # x = eq.DistributeOnDomain(self.x, eLowerBound)
        # eq.Residual = -self.DA()*(d(self.cA(x),self.x,eCFDM) - (1/(self.R()*self.T()))*self.F()*self.cA(x)*d(self.psi(x),self.x,eCFDM)) \
        # + (self.jR()*self.cA(x)*np.exp(-0.5*self.F()*self.psiStern()/(self.R()*self.T())) -
        #  self.jO()*np.exp(+0.5*self.F()*self.psiStern()/(self.R()*self.T())))

        # self.END_IF()

        #Right side BC

        eq = self.CreateEquation("fluxCond0Anions", "no anion flux x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = -d(self.cA(x),self.x,eCFDM) + (1/(self.R()*self.T()))*self.F()*self.cA(x)*d(self.psi(x),self.x,eCFDM) 
         
        # eq.Residual =self.JA

        self.IF(Time() < Constant(RC*300*s))
        eq = self.CreateEquation("fluxCond0Cations", "no cation flux x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = -d(self.cC(x),self.x,eCFDM) - (1/(self.R()*self.T()))*self.F()*self.cC(x)*d(self.psi(x),self.x,eCFDM)

        self.ELSE()
        eq = self.CreateEquation("fluxReduction", "cation flux x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = -self.DC()*(d(self.cC(x),self.x,eCFDM) - (1/(self.R()*self.T()))*self.F()*self.cC(x)*d(self.psi(x),self.x,eCFDM)) \
        - ( self.jO()*np.exp(+0.5*self.F()*self.psiStern()/(self.R()*self.T())) -
         self.jR()*(self.cC(x)/self.cC0())*np.exp(-0.5*self.F()*self.psiStern()/(self.R()*self.T())) )

        self.END_IF()


        
class simPoissonEq(daeSimulation):
    def __init__(self):
        daeSimulation.__init__(self)
        self.m = PoissonEq("name")

    def SetUpParametersAndDomains(self):
        
        DomainLenght = LambdaD*10
        DL = DomainLenght
        self.m.x.CreateStructuredGrid(Npoints,0,DL)

        pointsArray = np.empty(Npoints+1)
        pointsArray[:]=DL*(np.log(np.linspace(1,100,Npoints+1))/np.log(100))
        new_grid = pointsArray
        # self.m.x.Points = new_grid

        # self.m.psi0.SetValue(psi0V * V)
        self.m.psiStern.SetValue(psiStern * V)
        self.m.xStern.SetValue(xStern * m)
        
        self.m.jR.SetValue(jR * mol /(m**2*s))
        self.m.jO.SetValue(jO * mol /(m**2*s))

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

        for x in range(1,self.m.x.NumberOfPoints -1):
            self.m.cC.SetInitialCondition(x,cC0 * mol*(m**(-3)))
            self.m.cA.SetInitialCondition(x,cA0 * mol*(m**(-3)))
            self.m.psi.SetInitialGuess(x,0*V)

log = daePythonStdOutLog()
daesolver = daeIDAS()
datareporter = daeTCPIPDataReporter()
simulation = simPoissonEq()


simulation.m.SetReportingOn(True)

simulation.ReportingInterval = RC/10
simulation.TimeHorizon = RC*600

datareporter.Connect("", "PoissonEq")

simulation.Initialize(daesolver,datareporter,log)
simulation.SolveInitial()
simulation.Run()
simulation.Finalize()

print("analytical calculation for psiDL is: " , psiDL)
print("jR " , jR)
print("jO " , jO)