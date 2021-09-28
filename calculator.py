import numpy as np
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit
import os, math
from scipy.interpolate import interp2d
from scipy import ndimage, misc
from utils import tand, sind, cosd, gaussianDistributionFunction


class spiralArmMotionCalculator:
    def __init__(self, imgPath):
        self.imgPath = imgPath
        self.img = fits.open(imgPath)[0].data


    def converting(self, PA=0,inclination=0,final_radius=100,pw=360,r0=40,times=2,center=None):
        theta_ , R_ = np.meshgrid(np.linspace(0, 2*np.pi, pw),
                                np.arange(0, final_radius))
        if center == None:
            center = [len(self.img)/2, len(self.img)/2]

        x = np.linspace(0, len(self.img)-1, len(self.img), dtype=int)
        y = np.linspace(0, len(self.img)-1, len(self.img), dtype=int)
        f = interp2d(x,y,self.img,kind='cubic')
        z2 = [[[]for i in range(pw)] for i in range(final_radius)]
        theta = 0
        while theta < 90:
            r = 1
            while r < final_radius:
                rr = (r * ((sind(theta)*cosd(inclination))**2 + (cosd(theta))**2)**0.5)
                x = rr * cosd(np.arctan(tand(theta)*cosd(inclination))*180/math.pi + PA) + center[0]
                y = rr * sind(np.arctan(tand(theta)*cosd(inclination))*180/math.pi + PA) + center[1]
                z2[r][theta] = float(f(x,y))*((rr/r0)**times)
                r += 1
            theta += 1
        theta = 90
        r = 1
        while r < final_radius:
            rr = (r * ((sind(theta)*cosd(inclination))**2 + (cosd(theta))**2)**0.5)
            x = (r * (cosd(inclination))*cosd(90+PA)) + center[0]
            y = (r * (cosd(inclination))*sind(90+PA)) + center[1]
            z2[r][theta] = float(f(x,y))*((rr/r0)**times)
            r += 1


        theta = 91
        while theta <= 180:
            r = 1
            while r < final_radius:
                rr = (r * ((sind(theta)*cosd(inclination))**2+(cosd(theta))**2)**0.5)
                x = rr * cosd(np.arctan(tand(theta)*cosd(inclination))*180/math.pi+180+PA) + center[0]
                y = rr * sind(np.arctan(tand(theta)*cosd(inclination))*180/math.pi+180+PA) + center[1]
                z2[r][theta] = float(f(x,y))*((rr/r0)**times)

                r += 1
            theta += 1

        while theta < 270:
            r = 1

            while r < final_radius:
                rr = (r * ((sind(theta)*cosd(inclination))**2+(cosd(theta))**2)**0.5)
                x = rr * cosd(np.arctan(tand(theta)*cosd(inclination))*180/math.pi+180+PA) + center[0]
                y = rr * sind(np.arctan(tand(theta)*cosd(inclination))*180/math.pi+180+PA) + center[1]
                z2[r][theta] = float(f(x,y))*((rr/r0)**times)

                r += 1
            theta += 1

        theta = 270
        r = 1
        while r < final_radius:
            rr = (r * ((sind(theta)*cosd(inclination))**2+(cosd(theta))**2)**0.5)
            x =- (r * (cosd(inclination))*cosd(90+PA)) + center[0]
            y =- (r * (cosd(inclination))*sind(90+PA)) + center[1]
            z2[r][theta] = float(f(x,y))*((rr/r0)**times)
            r += 1
        theta = 271
        while theta < 360:
            r = 1
            while r < final_radius:
                rr = (r * ((sind(theta)*cosd(inclination))**2+(cosd(theta))**2)**0.5)
                x = (r * ((sind(theta)*cosd(inclination))**2+(cosd(theta))**2)**0.5) * cosd(np.arctan(tand(theta)*cosd(inclination))*180/math.pi+PA) + center[0]
                y = (r * ((sind(theta)*cosd(inclination))**2+(cosd(theta))**2)**0.5) * sind(np.arctan(tand(theta)*cosd(inclination))*180/math.pi+PA) + center[1]
                z2[r][theta] = float(f(x,y))*((rr/r0)**times)
                r += 1
            theta += 1

        r = 0
        theta = 0
        while theta < 360:
            z2[r][theta] = float(f(center[0],center[1]))
            theta += 1
        plt.imshow(z2, cmap= 'gray',origin='lower')
        return z2


    @staticmethod
    def fitForCenter(img_p,x_init,x_fin,y_init,y_fin,width=10,mini=None,
                   theta=None,rad=None,prob=None):
        if theta == None and rad == None and prob == None:
            theta = []
            rad = []
            prob = []
        elif theta == None or rad == None or prob == None:
            raise ValueError('please input theta, radius and probability together')
        w = width
        j = x_init
        while (j <= x_fin):
            X1 = []
            Y1 = []
            i = y_init
            max1 = 0
            a = 0
            while i < y_fin:
                if img_p[i][j] >= 0:
                    y = img_p[i][j]
                    X1.append(i)
                    Y1.append(y)
                    if img_p[i][j] > max1:
                        max1 = img_p[i][j]
                        a = i          #a is the location of the maximum
                i += 1
            xdata = X1
            ydata = Y1
            b = img_p[a][j]
            c = w/2
            d1 = int(min(X1))
            cen = Y1.index(max(Y1))
            w = width
            popt, pcov = curve_fit(gaussianDistributionFunction,xdata[cen-w:cen+w],ydata[cen-w:cen+w],[b-d1,c,a,d1],maxfev = 5000000)

            if (popt[2] < max(X1)) & (popt[2] > min(X1)):
                rad.append(popt[2])
                theta.append(j)
                p = np.sqrt(pcov[2][2])
                prob.append(p)
            j += 1
        return theta, rad, prob


    @staticmethod
    def fitting(theta1,theta2,rad1,rad2,prob1,prob2,est=-0.1,withRn=True):
        def f_withRn(X,rp,n,a0,a1,a2,a3,a4,a5):
            z,y,k = X
            x = z + y*rp*(k/60)**n
            return a0+a1*x**1+a2*(x)**2+a3*(x)**3+a4*x**4+a5*x**5

        def f_withoutRn(X,rp,n,a0,a1,a2,a3,a4,a5):
            z,y = X
            x = z + y*rp
            return a0+a1*x**1+a2*(x)**2+a3*(x)**3+a4*x**4+a5*x**5

        theta = np.hstack((theta1,theta2))
        theta1b = np.zeros(len(theta1))
        theta2b = np.ones(len(theta2))

        thetab = np.hstack((theta1b,theta2b))
        rad = np.hstack((rad1,rad2))
        P = np.hstack((prob1,prob2))

        if withRn == True:
            thetaf = np.vstack((theta,thetab,rad))
            popt, pcov = curve_fit(f_withRn,thetaf,rad,[est]*8,sigma=P,absolute_sigma=True,maxfev=5000000)

        else:
            thetaf = np.vstack((theta,thetab))
            popt, pcov = curve_fit(f_withoutRn,thetaf,rad,[est]*7,sigma=P,absolute_sigma=True,maxfev=5000000)

        return popt, pcov
