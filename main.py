import numpy as np
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
import os, math
from scipy.interpolate import interp2d
from scipy import ndimage, misc
from calculator import spiralArmMotionCalculator


calculatorA = spiralArmMotionCalculator('x1a.fits')
poa = calculatorA.converting(-50,0,150,360,40,2)
calculatorB = spiralArmMotionCalculator('x1b.fits')
pob = calculatorB.converting(-50,0,150,360,40,2)

Z1,K1,P1=calculatorA.fitForCenter(poa,0,100,50,85,width=10,theta=None,rad=None,prob=None)
Z2,K2,P2=calculatorB.fitForCenter(pob,0,100,50,85,width=10,theta=None,rad=None,prob=None)

popt, pcov=spiralArmMotionCalculator.fitting(Z1,Z2,K1,K2,P1,P2,est=-0.1,withRn=True)
print(popt)
