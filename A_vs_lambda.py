#!/usr/bin/env python
 
import re,sys,string,math,os,types,shutil
#sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
#sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/LESPEC')
#sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
from texttable import txttableclass
from lightechoprocs import *
from sb import CCMextinctionA

if __name__=='__main__':
    Rvs = [2.0,3.1,5.5]
    EBmV = 0.1
    tout = txttableclass()
    tout.configcols(['wave'],'f','%.0f',visible=1)
    for Rv in Rvs:
        tout.configcols(['A_Rv%.1f' % (Rv)],'f','%.3e',visible=1)    
    for Rv in Rvs:
        tout.configcols(['An_Rv%.1f' % (Rv)],'f','%.3f',visible=1)    
   
    for w in range(3500,10000,5):
        key = tout.newrow({'wave':w})        
        for Rv in Rvs:
            A_lambda = CCMextinctionA(w,Rv,EBmV)
            An_lambda = A_lambda/CCMextinctionA(8000.0,Rv,EBmV)
            tout.setentry(key,'A_Rv%.1f' % (Rv),A_lambda)
            tout.setentry(key,'An_Rv%.1f' % (Rv),An_lambda)
    outfilename = 'A_vs_lambda.txt'
    print(f'Saving {outfilename}')
    tout.save2file(outfilename)
