#!/usr/bin/env python
 
import re,sys,string,math,os,types,shutil
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/LESPEC')
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
from texttable import txttableclass
from lightechoprocs import *
from sb import S_tableclass

if __name__=='__main__':


    dusttypes = ['LMCavg','LMC2','SMCbar','MWG']
    ws = [4000.0,6000.0,7999.0]

    tout = txttableclass()
    tout.configcols(['theta'],'f','%.0f',visible=1)
    tout.configcols(['rho_t','z_t','r_t'],'f','%.3f',visible=1)

    S = {}
    for dust in dusttypes:
        S[dust]=S_tableclass()
        dustfilename =  S[dust].getdustfilename(dust)
        print(f'Loading dust properties from {dustfilename}...')
        S[dust].loadtable(dustfilename)
        for w in ws:
            tout.configcols(['S_%s_%.0f' % (dust,w)],'f','%.3e',visible=1)    
        for w in ws:
            tout.configcols(['Sn_%s_%.0f' % (dust,w)],'f','%.3f',visible=1)    

    for theta in range(1,180):
        tan_theta = math.tan(theta*deg2rad)
        inverse_tan_theta = 1.0/tan_theta
        rho_t = inverse_tan_theta + math.sqrt(inverse_tan_theta**2 + 1)
        z_t = 0.5*(rho_t*rho_t - 1.0)
        r_t = math.sqrt(rho_t*rho_t + z_t*z_t)
        key = tout.newrow({'theta':theta,
                           'rho_t':rho_t,
                           'z_t':z_t,
                           'r_t':r_t})
        
        for dust in dusttypes:
            for w in ws:
                Sval = S[dust].S_cm2(theta,w)
                tout.setentry(key,'S_%s_%.0f' % (dust,w),Sval)
                tout.setentry(key,'Sn_%s_%.0f' % (dust,w),Sval/S[dust].S_cm2(90.0,w))

        #rho2_t = 0.5*inverse_tan_theta - math.sqrt((inverse_tan_theta**2)/4 +0.5)
        #print '%.0f %.3f %.3f' % (theta,rho_t)
    outfilename = 'S_vs_theta.txt'
    print(f'Saving {outfilename}')
    tout.save2file(outfilename)
