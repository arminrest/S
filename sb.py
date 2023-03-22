#!/usr/bin/env python
 
import re,sys,string,math,os,types,shutil
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/LESPEC')
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
from texttable import txttableclass
from pipeclasses import paramfileclass
import optparse
#import matplotlib
#matplotlib.use('GTKAgg')
#matplotlib.use('Agg')
from lightechoprocs import *
import pylab
from tools import rmfile,makepath4file
import numpy as np
from scipy import interpolate


# LMC/SMC extinction law
# MCextinctionA = <A(lambda)/A(V)> * R_V * E(B-V)
# flux = flux0 * 10^(-0.4*CCMextinctionA)
# from Gordon, CLayton, Misselt, Landoldt, Wolff 2003
# <A(lambda)/A(V)> is equation 5
# possible predefined dusttypes in 'setmodelparams' are:
# SMCbar, SMCwing, LMC2, LMCave, see Table 3
class MCextinctionclass:
    def __init__(self,dustmodel=None):
        if dustmodel==None:
            dustmodel='LMCave'
        self.setmodelparams(dustmodel)
    def setmodelparams(self,dustmodel):
        if  dustmodel == 'SMCbar':
            # Table 3
            self.c1=-4.959
            self.c2= 2.264
            self.c3= 0.389
            self.c4= 0.461
            self.x0= 4.6
            self.gamma= 1.0
            # Table 2
            self.Rv = 2.74
        elif  dustmodel == 'SMCwing':
            self.c1= -0.856
            self.c2= 1.038
            self.c3= 3.215
            self.c4= 0.107
            self.x0= 4.703
            self.gamma= 1.212
            # Table 2
            self.Rv = 2.05
        elif dustmodel == 'LMC2':
            self.c1=-1.475
            self.c2= 1.132
            self.c3= 1.463
            self.c4= 0.294
            self.x0= 4.558
            self.gamma= 0.945
            # Table 2
            self.Rv = 2.76
        elif dustmodel == 'LMCave':
            self.c1=-0.890
            self.c2= 0.998
            self.c3= 2.719
            self.c4= 0.400
            self.x0= 4.579
            self.gamma= 0.934
            # Table 2
            self.Rv = 3.41
        elif dustmodel == 'Sk-69 213':
            self.c1=-2.791
            self.c2= 1.594
            self.c3= 1.816
            self.c4= 0.527
            self.x0= 4.564
            self.gamma= 0.735
            # Table 2
            self.Rv = 3.96
        else:
            raise RuntimeError('dust model %s is not defined!' % dustmodel)
            
    def Ax_AV(self,x):
        # Gordon, CLayton, Misselt, Landoldt, Wolff 2003, equation 5
        # x is in micron^-1
        x2 = x*x
        d2 = x2 - self.x0*self.x0
        D = x2/(d2*d2 + x2*self.gamma*self.gamma)
        if x>=5.9:
            x59 = (x-5.9)
            F = 0.5392 * x59*x59 + 0.05644 * x59*x59*x59
        else: 
            F = 0.0
        Ax_AV_value = 1/self.Rv*(self.c1 + self.Rv + self.c2*x + self.c3*D + self.c4*F)
        return(Ax_AV_value)
        
    def MCextinctionA(self,lambdaval,EBmV):
        #return(0.0)
        # lmabdaval is in Angstrom
        inverselambda = 1.0/(lambdaval*1E-4) # convert lambda(angstroem) into lambda^-1 in microns: inverselambda = 1.0/(lambda/10000.0)
        return(self.Ax_AV(inverselambda) * EBmV * self.Rv)

# Galactic extinction law
# CCMextinctionA = <A(lambda)/A(V)> * R_V * E(B-V)
# flux = flux0 * 10^(-0.4*CCMextinctionA)
# from Cardelli, Clayton, Mathis 89
def CCMextinctionA(lambdaval,Rv,EBmV):
    """
    flux = flux0 * 10^(-0.4*CCMextinctionA)
    CCMextinctionA = (a + b/Rv) * EBmV* Rv
    from Cardelli, Clayton, Mathis 89
    """
    inverselambda = 1.0/(lambdaval*1E-4) # convert lambda(angstroem) into lambda^-1 in microns: inverselambda = 1.0/(lambda/10000.0)
    
    if( (inverselambda >= 0.3) and (inverselambda < 1.1)):
        a = 0.574* math.pow(inverselambda,1.61)
        b = -0.527* math.pow(inverselambda,1.61)
    elif ((inverselambda>=1.1) and (inverselambda < 3.3)):
        y = inverselambda-1.82
        y2 = y*y
        y3 = y2*y
        a =  1.0 + 0.17699*y - 0.50447*y2 - 0.02427*y3 + 0.72085*y2*y2 + 0.01979*y2*y3 - 0.77530*y3*y3 + 0.32999*y3*y3*y;
        b =        1.41338*y + 2.28305*y2 + 1.07233*y3 - 5.38434*y2*y2 - 0.62251*y2*y3 + 5.30260*y3*y3 - 2.09002*y3*y3*y;
    elif((inverselambda>=3.3) and (inverselambda < 8.0) ):
        if(inverselambda < 5.9):
            a = 1.752 -0.316*inverselambda - 0.104/(math.pow(inverselambda-4.67,2) + 0.341)
            b = -3.090 + 1.825*inverselambda + 1.206/(math.pow(inverselambda-4.62,2) + 0.263)
        else:
            a = 1.752 - 0.316*inverselambda - 0.104/(math.pow(inverselambda-4.67,2) + 0.341) - 0.04473*math.pow(inverselambda-5.9,2) - 0.009779*math.pow(inverselambda-5.9,3)
            b = -3.090 + 1.825*inverselambda + 1.206/(math.pow(inverselambda-4.62,2) + 0.263) + 0.2130*math.pow(inverselambda-5.9,2)  + 0.1207*math.pow(inverselambda-5.9,3)
    else:
        raise RuntimeError('inverse lambda out of range: %f (lambda = %f)' % (inverselambda,lambdaval))
    aratio = a + b/Rv;
    return(EBmV * Rv * aratio)
    
class S_tableclass:
    def __init__(self):
        self.S=None
        self.lambda_A=None
        self.theta_deg=None

        self.normfactor=1.0
        
    def loadtable(self,filename):
        print('Loading ',filename)
        lines = open(filename,'r').readlines()
        
        lambdarange = None
        thetarange = None
        
        self.lambda_A = None
        self.S = []
        self.theta_deg = []
        
        
        thetacounter = 0
        for line in lines:
            if not re.match('\#',line):
                if lambdarange == None or thetarange == None:
                    print('ERROR: Could not determine lambdarange or thetarange!')
                    sys.exit(0)
                if self.lambda_A == None:
                    print('ERROR: Could not determine lambda!')
                    sys.exit(0)
                linedummy = re.sub('^\s+|\s+$','',line,count=2)
                Stmp = re.split('\s+',linedummy)
                Stmp = [float(x) for x in Stmp] # convert to float
                self.theta_deg.append(Stmp[0]) # the first column is theta
                self.S.append(Stmp[1:])        # the rest is S
                if len(self.S[thetacounter]) != len(self.lambda_A):
                    print('Bug! inconsistent number of S vals (%d!=%d) in data line %d' % (len(self.S[thetacounter]),len(self.lambda_A),thetacounter))
                    sys.exit(0)                                                    
                thetacounter+=1
            elif re.match('^\#lambdarange:\s*',line):
                linedummy = re.sub('^\#lambdarange:\s*|\s+$','',line,count=2)
                lambdarange = re.split('\s+',linedummy)
                if len(lambdarange) != 3:
                    print('ERROR! expecting three values for lambdarange!',lambdarange)
                    sys.exit(0)                                                    
                self.lambdamin  = float(lambdarange[0])
                self.lambdamax  = float(lambdarange[1])
                self.lambdastep = float(lambdarange[2])
            elif re.match('^\#thetarange:\s*',line):
                linedummy = re.sub('^\#thetarange:\s*|\s+$','',line,count=2)
                thetarange = re.split('\s+',linedummy)
                if len(thetarange) != 3:
                    print('ERROR! expecting three values for thetarange!',thetarange)
                    sys.exit(0)                                                    
                self.thetamin  = float(thetarange[0])
                self.thetamax  = float(thetarange[1])
                self.thetastep = float(thetarange[2])
            elif re.match('^\#theta/lambda\s*',line):
                linedummy = re.sub('^\#theta/lambda\s*|\s+$','',line,count=2)
                self.lambda_A = re.split('\s+',linedummy)
                self.lambda_A = [float(x) for x in self.lambda_A] # convert to float
               
        if len(self.lambda_A) != int((self.lambdamax-self.lambdamin)/self.lambdastep+1.0):
            print('Bug! number of lambdas inconsistent with range! %d != %d\n' % (len(self.lambda_A),int((self.lambdamax-self.lambdamin)/self.lambdastep+1.0)))
            sys.exit(0)

    def norm1(self,thetaval,lambdaval):
        Sval = self.S_cm2(thetaval,lambdaval,normfactor=1.0)* 1E22
        self.normfactor = 1.0/Sval
        print('Setting normalization factor of scattering S to %f' % (self.normfactor))

    def S_cm2(self,thetaval,lambdaval,normfactor=None):
        thetaindex  = int((thetaval - self.thetamin)/self.thetastep)
        if thetaindex<0 or thetaindex>=len(self.theta_deg)-1:
            print('ERROR! theta %f out of range! index = %d' % (thetaval,thetaindex))
            sys.exit(0)                                                    
        
        if (thetaval<self.theta_deg[thetaindex] or thetaval>=self.theta_deg[thetaindex+1]):
            print('ERROR! theta %f not between theta %f and %f' % (thetaval,self.theta_deg[thetaindex],self.theta_deg[thetaindex+1]))
            sys.exit(0)

        lambdaindex = int((lambdaval - self.lambdamin)/self.lambdastep)
        if lambdaindex<0 or lambdaindex>=len(self.lambda_A)-1:
            print('ERROR! lambda %f out of range! index = %d' % (lambdaval,lambdaindex))
            sys.exit(0)

        if (lambdaval<self.lambda_A[lambdaindex] or lambdaval>=self.lambda_A[lambdaindex+1]):
            print('ERROR! lambda %f not between lambda %f and %f' % (lambdaval,self.lambda_A[lambdaindex],self.lambda_A[lambdaindex+1]))
            sys.exit(0)

        # calculated the weighted mean of the 4 corners. The weight is 
        weight_lambda = (lambdaval-self.lambda_A[lambdaindex])/(self.lambda_A[lambdaindex+1]-self.lambda_A[lambdaindex])
        weight_theta = (thetaval-self.theta_deg[thetaindex])/(self.theta_deg[thetaindex+1]-self.theta_deg[thetaindex])

        Sval  = self.S[thetaindex  ][lambdaindex  ] * (1.0-weight_lambda) * (1.0-weight_theta)  
        Sval += self.S[thetaindex  ][lambdaindex+1] * (    weight_lambda) * (1.0-weight_theta)  
        Sval += self.S[thetaindex+1][lambdaindex+1] * (    weight_lambda) * (    weight_theta)  
        Sval += self.S[thetaindex+1][lambdaindex  ] * (1.0-weight_lambda) * (    weight_theta)  

        #print '%d: theta: %f' % (thetaindex,self.theta_deg[thetaindex])
        #print '%d: lambda: %f' % (lambdaindex,self.lambda_A[lambdaindex])
        if normfactor==None:
            Sval *= self.normfactor
        else:
            Sval *= normfactor            
        
        return(Sval)

    def getdustfilename(self,dust):
        dustfilename  = '%s/' % (os.environ['DUST_ROOTDIR'])
        if dust == 'LMCavg':
            dustfilename += 'S.LMCavg.bC2.0.dat'
        elif dust == 'LMC2':
            dustfilename += 'S.LMC2.bC1.0.dat'
        elif dust == 'SMCbar':
            dustfilename += 'S.SMC.bC0.0.dat'
        elif dust == 'MWG':
            dustfilename += 'S.MWG_A.bC5.6RV3.1.dat'
        else:
            raise RuntimeError('dust %s not defined!' % (dust))
            
        if not os.path.isfile(dustfilename):
            raise RuntimeError('Could not find dust file %s' % (dustfilename))
        return(dustfilename)
if __name__=='__main__':

    S_table          = S_tableclass()

    dustfilename = S_table.getdustfilename('MWG')
    print('Loading dust properties...')
    S_table.loadtable(dustfilename)

    tout = txttableclass()
    tout.configcols(['theta'],'f','%.0f',visible=1)
    tout.configcols(['rho_t','z_t','r_t','1/(rho_t*r_t)','1/(r_t)','1/(r_t*r_t)','1/(r_t*r_t*r_t)','z_over_rho','inv_tan_theta','tan_theta'],'f','%.3f',visible=1)
    tout.configcols(['w4000','w6000','w8000'],'f','%.3e',visible=1)
    tout.configcols(['w4000n','w6000n','w8000n'],'f','%.3f',visible=1)    
    tout.configcols(['w4000n_rhor','w6000n_rhor','w8000n_rhor'],'f','%.3f',visible=1)    
    tout.configcols(['w4000n_r','w6000n_r','w8000n_r'],'f','%.3f',visible=1)    
    tout.configcols(['w4000n_rr','w6000n_rr','w8000n_rr'],'f','%.3f',visible=1)    
    tout.configcols(['w4000n_rrr','w6000n_rrr','w8000n_rrr'],'f','%.3f',visible=1)    

    w4000_90deg = S_table.S_cm2(90,4000.0)
    w6000_90deg = S_table.S_cm2(90,6000.0)
    w8000_90deg = S_table.S_cm2(90,8000.0)
    for theta in xrange(1,180):
        tan_theta = math.tan(theta*deg2rad)
        inverse_tan_theta = 1.0/tan_theta
        #rho_t = 0.5*inverse_tan_theta + math.sqrt((inverse_tan_theta**2)/4.0 +0.5)
        rho_t = inverse_tan_theta + math.sqrt(inverse_tan_theta**2 + 1)
        z_t = 0.5*(rho_t*rho_t - 1.0)
        r_t = math.sqrt(rho_t*rho_t + z_t*z_t)
        inv_rho_t_r_t = 1.0/(rho_t*r_t)
        inv_r = 1.0/(r_t)
        inv_rr = 1.0/(r_t*r_t)
        inv_rrr = 1.0/(r_t*r_t*r_t)
        w4000 = S_table.S_cm2(theta,4000.0)
        w6000 = S_table.S_cm2(theta,6000.0)
        w8000 = S_table.S_cm2(theta,8000.0)
        
        key = tout.newrow({'theta':theta,
                           'rho_t':rho_t,
                           'z_t':z_t,
                           'r_t':r_t,
                           '1/(rho_t*r_t)':inv_rho_t_r_t,
                           '1/(r_t)':inv_rr,
                           '1/(r_t*r_t)':inv_rr,
                           '1/(r_t*r_t*r_t)':inv_rrr,
                           'z_over_rho':z_t/rho_t,
                           'tan_theta':tan_theta,
                           'inv_tan_theta':inverse_tan_theta,
                           'w4000':w4000,
                           'w6000':w6000,
                           'w8000':w8000,
                           'w4000n':w4000/w4000_90deg,
                           'w6000n':w6000/w6000_90deg,
                           'w8000n':w8000/w8000_90deg,
                           'w4000n_rhor':w4000/w4000_90deg*inv_rho_t_r_t,
                           'w6000n_rhor':w6000/w6000_90deg*inv_rho_t_r_t,
                           'w8000n_rhor':w8000/w8000_90deg*inv_rho_t_r_t,
                           'w4000n_r':w4000/w4000_90deg*inv_r,
                           'w6000n_r':w6000/w6000_90deg*inv_r,
                           'w8000n_r':w8000/w8000_90deg*inv_r,
                           'w4000n_rr':w4000/w4000_90deg*inv_rr,
                           'w6000n_rr':w6000/w6000_90deg*inv_rr,
                           'w8000n_rr':w8000/w8000_90deg*inv_rr,
                           'w4000n_rrr':w4000/w4000_90deg*inv_rrr,
                           'w6000n_rrr':w6000/w6000_90deg*inv_rrr,
                           'w8000n_rrr':w8000/w8000_90deg*inv_rrr
                           })
        #rho2_t = 0.5*inverse_tan_theta - math.sqrt((inverse_tan_theta**2)/4 +0.5)
        #print '%.0f %.3f %.3f' % (theta,rho_t)
    tout.save2file('S_vs_theta')
