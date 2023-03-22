#!/usr/bin/env python
 
import re,sys,string,math,os,types,shutil
#sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
#sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/LESPEC')
#sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
from texttable import txttableclass
from lightechoprocs import *
from sb import S_tableclass,CCMextinctionA

if __name__=='__main__':

    dusttypes = ['LMCavg','LMC2','SMCbar','MWG']
    thetas = [15.0,20.0,40.0,60.0,90.0,130.0,150.0,179.0]

    tout = txttableclass()
    tout.configcols(['wave'],'f','%.0f',visible=1)

    S = {}
    wmin={'LMCavg':4000.0,'LMC2':4000.0,'SMCbar':4000.0,'MWG':3500.0}
    wmax={'LMCavg':7999.0,'LMC2':7999.0,'SMCbar':7999.0,'MWG':9999.0}

    for dust in dusttypes:
        S[dust]=S_tableclass()
        dustfilename =  S[dust].getdustfilename(dust)
        print(f'Loading dust properties from {dustfilename}...')
        S[dust].loadtable(dustfilename)
        for theta in thetas:
            tout.configcols(['S_%s_%.0f' % (dust,theta)],'f','%.3e',visible=1)    
        for theta in thetas:
            tout.configcols(['Sn_%s_%.0f' % (dust,theta)],'f','%.3f',visible=1)    

    for w in range(3500,10000,5):
        key = tout.newrow({'wave':w})        
        for dust in dusttypes:
            if w<wmin[dust] or w>wmax[dust]: continue
            for theta in thetas:
                Sval = S[dust].S_cm2(theta,w)
                tout.setentry(key,'S_%s_%.0f' % (dust,theta),Sval)
                tout.setentry(key,'Sn_%s_%.0f' % (dust,theta),Sval/S[dust].S_cm2(theta,7999.0))

        #rho2_t = 0.5*inverse_tan_theta - math.sqrt((inverse_tan_theta**2)/4 +0.5)
        #print '%.0f %.3f %.3f' % (theta,rho_t)
    outfilename = 'S_vs_lambda.txt'
    print(f'Saving {outfilename}')
    tout.save2file(outfilename)
