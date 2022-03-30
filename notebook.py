import numpy as np
from scipy import io, integrate, linalg, signal
from scipy.sparse.linalg import eigs
from matplotlib import pyplot as plt
import daetools
from daetools.pyDAE import *
from pyUnits import *
from pyCore import *

# pointsArray = np.zeros(1)
# for i in range(100):
#      pointsArray = np.append(pointsArray, pointsArray[i]+(100-pointsArray[i])/50)
# print(pointsArray)

number = 100
pointsArray = np.zeros(1)
for i in range(number):
     pointsArray = np.append(pointsArray, np.log(i+2))
pointsArray = (pointsArray/np.amax(pointsArray))*number
print(pointsArray)