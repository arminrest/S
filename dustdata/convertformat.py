#!/usr/bin/env python
 
import re,sys,string,math,os,types
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
import texttable

if __name__=='__main__':

    if len(sys.argv)<2:
        print 'Usage: convertformat.py datfile'
        print 'NOTE: multiplies lambda with 10000.0 to get Angstroems!!!!'
        print 'output format: lambda albedo      g      Cext         K'
        sys.exit(0)

    t = texttable.txttableclass()
    t.loadfile(sys.argv[1])
    t.configcols(['Cext','K'],'f','%.3e')
    t.configcols(['albedo','g'],'f','%.4f')
    t.configcols(['lambda'],'f','%.1f')
    for key in t.rowkeys():
        t.setentry(key,'lambda',t.getentry(key,'lambda')*10000.0)
    print 'Saving converted file to temp/%s' % sys.argv[1]
    t.save2file('temp/%s' % sys.argv[1])
