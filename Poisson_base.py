import numpy as np
from scipy import io, integrate, linalg, signal
from scipy.sparse.linalg import eigs
from matplotlib import pyplot as plt
import daetools
from daetools.pyDAE import *
from pyUnits import *
from pyCore import *

elec_pot_t = daeVariableType(name = "elec_pot_t", units = V, lowerBound = -1e20, upperBound = 1e20, initialGuess = 0, absTolerance = 1e-5)

# there are some discrepancies on the unit of maesures, but we can fix it easily, for the moment
# the problem is the initialization

class PoissonEq(daeModel):
    def __init__(self, Name, Parent = None, Description = ""):
        daeModel.__init__(self,Name,Parent,Description)

        self.x = daeDomain("x", self,m, "x axis domain")

        self.D = daeParameter("D", (m**2)/s, self, "Diffusion coefficinet")
        self.R = daeParameter("RT", (J/(mol*K)),self, "Gas constant")
        self.e = daeParameter("e", A*s,self, "electron charge" )
        self.eps = daeParameter("eps", F/m, self, "dielectric constant")
        self.c0 = daeParameter("c0", mol*(m**(-3)), self, "initial concentration")
        self.phi0 = daeParameter("phi0", V, self, "potential at the wall")
        self.Na = daeParameter("Na", mol**(-1), self, "avogadro's number")
        self.T = daeParameter("T", K, self, "temperature")

        self.rho = daeVariable("rho", molar_concentration_t, self)
        self.rho.DistributeOnDomain(self.x)
        self.rho.Description = "charge density"

        self.psi = daeVariable("psi", elec_pot_t, self)
        self.psi.DistributeOnDomain(self.x)
        self.psi.Description = "potential"

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        eq = self.CreateEquation("chargeEq", "Charge density equation")
        eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eOpenOpen)

        eq.Residual = dt(self.rho(x)) - self.D()*(d2(self.rho(x),self.x, eCFDM) + ((2*(self.e()**2)*self.c0()*self.Na()**2)/(self.eps()*self.R()*self.T()))*self.rho(x))

        eq = self.CreateEquation("PoissonEq", "Poisson equation")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eOpenOpen)

        eq.Residual = d2(self.psi(x),self.x, eCFDM) + (self.e()*self.Na()/self.eps())*self.rho(x)

#         boundary conditions
        eq = self.CreateEquation("DirichBC1", "0 potential at x=0")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.psi(x)

        eq = self.CreateEquation("DirichBC2", "fixed potential at x=L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = self.psi(x) + self.phi0()

        eq = self.CreateEquation("fluxCond0", "result from the 0 flux condition at x+L")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eUpperBound)
        eq.Residual = d(self.rho(x), self.x, eCFDM) + ((self.e()*self.c0()*self.Na())/(self.R()*self.T()))*d(self.psi(x),self.x, eCFDM)

        eq = self.CreateEquation("chargeinBulk0", "0 charge in bulk elyte condition")
        # eq.CheckUnitsConsistency = False
        x = eq.DistributeOnDomain(self.x, eLowerBound)
        eq.Residual = self.rho(x)

class simPoissonEq(daeSimulation):
    def __init__(self):
        daeSimulation.__init__(self)
        self.m = PoissonEq("name")


    def SetUpParametersAndDomains(self):
        self.m.x.CreateStructuredGrid(30,0,1)
    
        # self.m.x.CreateDistributed(eCFDM,2,n,0,0.1)
        # self.m.x.CreateArray(10)
        # self.m.x.Points = [x.value for x in range(30)]


        self.m.e.SetValue(1.6e-19 *A*s)
        self.m.eps.SetValue(78.5*8.85e-12 *(F/m))
        self.m.D.SetValue( 1e-9 *((m**2)/s))
        self.m.R.SetValue(8.314 *J/(mol*K))
        self.m.c0.SetValue(0.1* mol*(m**(-3)))
        self.m.phi0.SetValue(-0.01 *V)
        self.m.Na.SetValue(6.02e23 * mol**(-1))
        self.m.T.SetValue(298 * K)

    def SetUpVariables(self):
        for x in range(1,self.m.x.NumberOfPoints -1):
            self.m.rho.SetInitialGuess(x,0 * mol*(m**(-3)))
            self.m.psi.SetInitialCondition(x,0*V)



log = daePythonStdOutLog()
daesolver = daeIDAS()
datareporter = daeTCPIPDataReporter()
simulation = simPoissonEq()


simulation.m.SetReportingOn(True)

simulation.ReportingInterval = 30
simulation.TimeHorizon = 1000

datareporter.Connect("", "PoissonEq")

simulation.Initialize(daesolver,datareporter,log)
simulation.SolveInitial()
simulation.Run()
simulation.Finalize()

