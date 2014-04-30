import numpy as nm
import string as st
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *
import pylab as pl
import copy as cp
import minuit

class vars:
    def __init__(self):
        self.x = 10.0
        self.y = 15.0

def f(x,y):
    return ((x-2) / 3)**2 + y**2 + y**4

x=10
y=15
m = minuit.Minuit(f,x=10,y=15)



