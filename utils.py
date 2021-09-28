import math
def tand(x):
    return math.tan(x * math.pi / 180)


def sind(x):
    return math.sin(x * math.pi / 180)


def cosd(x):
    return math.cos(x * math.pi / 180)


def gaussianDistributionFunction(x,A1,sigma1,miu1,A0):#gaussian distribution function
    return A1*(1/((2**0.5)*sigma1))*math.e**((-(x-miu1)**2)/(2*sigma1**2)) + A0
