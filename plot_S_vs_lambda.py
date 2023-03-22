#!/usr/bin/env python
import re,sys,string,math,os,types,exceptions,time,fcntl,shutil
# put the tools directory into the path
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/LESPEC')
import matplotlib
#matplotlib.use('GTKAgg')
#matplotlib.use('Agg')
import pylab as matlib
from   matplotlib.ticker import FormatStrFormatter,MultipleLocator

from texttable import txttableclass 
#http://matplotlib.sourceforge.net/
import optparse
import copy

from plotit import mkplot,mkplot_options,createLegend


def mkxypanel(spxy,dustwidth,trho,trhofilename,options,legendflag=False,diagonalflag=False):

    matlib.setp(spxy.get_xticklabels(), visible=False)

    if options.getlimitsforequalaxis:
        xcolname = 'z_arcsec'
        #options.xmin = None
        #options.xmax = None
        options.axisequal = True
        options.major_ticks_x = None
    else:
        xcolname = 'z'
    options.x =xcolname
    
    plots4legend = []
    # make the dust
    if diagonalflag:
        delta_y =  0.165/math.sqrt(2.0)
        x = [0.3393,-0.3393]
        y1 = [-7+delta_y,7+delta_y]
        y2 = [-7-delta_y,7-delta_y]
        spxy.fill_between(x,y1,y2,color=dustcolor,alpha=dustalpha)
        spxy.plot(x, y1, 'k-',lw=0.2) 
        spxy.plot(x, y2, 'k-',lw=0.2) 
        # make a dummy axhspan for the event pulse legend
        p = spxy.axhspan(-6.0,-6.0,fill=True,facecolor=dustcolor,alpha=dustalpha)
        plots4legend.append(p)
    else:    
        p = spxy.axvspan(-0.5*dustwidth,+0.5*dustwidth,fill=True,facecolor=dustcolor,alpha=dustalpha)
        plots4legend.append(p)


    # make the event
    t0 = 300.0
    dt = [0,30,-100]
    tp = []
    for i in xrange(len(dt)):
        tp.append(t0+dt[i]/365.25)

    x = trho.col_asfloat_list(xcolname)
    y1 = trho.col_asfloat_list('rho(%.3f)' % tp[2])
    y2 = trho.col_asfloat_list('rho(%.3f)' % tp[1])
    
    spxy.fill_between(x,y1,y2,color=eventcolor,alpha=eventalpha)

    yFormatter   = FormatStrFormatter('%3.1f')
    spxy.yaxis.set_major_formatter( yFormatter )

        
    # make a dummy axhspan for the event pulse legend
    p = spxy.axhspan(-6.0,-6.0,fill=True,facecolor=eventcolor,alpha=eventalpha)
    plots4legend.append(p)

    # 
    options.y = 'rho(%.3f),rho(%.3f),rho(%.3f)' % (tp[0],tp[1],tp[2])
    (spxy,dummy1,dummy2) = mkplot(spxy,[trhofilename],options)
    #sp.text(0.125,301.8,'A',size=16)
    #spxy.grid(True)

    spxy.set_ylim((options.ymin, options.ymax))

    if legendflag:
        names = ['Dust','Light Pulse']
        createLegend(spxy,plots4legend,names,loc='upper left', fontsize=options.labelfontsize,labelspacing=options.legendlabelspacing,borderpad=options.legendborderpad,handletextpad=options.legendhandletextpad)

    if options.getlimitsforequalaxis:
        matlib.show()
        distance_ly= 10000.0
        print 'THESE ARE THE AXIS VALUES TO GET EQUAL AREA ASSUMING A DISTANCE OF %f light years!!!!' % distance_ly
        arcsec2ly = distance_ly * math.tan(1.0/3600.0 * deg2rad)
        xmin,xmax = spxy.get_xlim()
        xmin_ly = xmin*arcsec2ly
        xmax_ly = xmax*arcsec2ly
        print 'xxx',xmin,xmax,' arcsec'
        print 'xxx',xmin_ly,xmax_ly,' ly, add --xmin %.4f --xmax %.4f' % (xmin_ly,xmax_ly)
        print 'yyy',spxy.get_ylim(),' arcsec'

        
        print 0.008,'ly is ',0.008/arcsec2ly,'"'
        print 7,'" is ',7*arcsec2ly,'ly'

        sys.exit(0)



def mkfluxpanel(spflux,filelist,phasemin_day,phasemax_day,options,ylabel=False,legendflag=False,stretch=None,psfsize=None):
    
    matlib.setp(spflux.get_yticklabels(), visible=False)
    matlib.setp(spflux.get_xticklabels(), visible=False)

    ax2 = spflux.twinx()
    ax2.figure.canvas.draw()
    matlib.setp(ax2.get_xticklabels(), visible=False)


    ax2options = copy.deepcopy(options)
    if stretch!=None:
        #ax2options.mult2y = '%f' % stretch        
        phasemin_day /= stretch
        phasemax_day /= stretch
    ax2options.x = 'flux'
    ax2options.y = 'dt'
    ax2options.mult2x = '97'
    ax2options.format = 'r-'
    ax2options.color = 'CBred'
    ax2options.nolegend = True
    ax2options.ymin=phasemax_day
    ax2options.ymax=phasemin_day
    ax2options.major_ticks_y = 50
    ax2options.minor_ticks_y = 5    


    (ax2,lcplot,dummy2)=mkplot(ax2,['sn1993j.rest.interp.lc.verysmallsteps'],ax2options)

    #parser1 = mkplot_options()
    #lcoptions,  dummy3 = parser1.parse_args(args=re.split('\s+','-x flux -y dt --mult2x 100  -f r: --nolegend'))
    #(ax2,lcplot,dummy2)=mkplot(ax2,['sn1993j.rest.interp.lc.verysmallsteps'],lcoptions)
    #ax2.set_ylim(phasemax_day,phasemin_day)
    

    if ylabel:
        ax2.set_ylabel('phase in days',rotation=-90.0)
    else:
        ax2.set_ylabel(' ',rotation=-90.0)
        matlib.setp(ax2.get_yticklabels(), visible=False)
        
    (spflux,dummy1,dummy2) = mkplot(spflux,filelist,options)
    blacklineplot = dummy1[0]

    if psfsize!=None:
        spflux.text(80,-3,'PSF: %.2f"'% psfsize,size=16)    

    if legendflag:
        names = ['LE flux profile','SN 1993J lightcurve']
        createLegend(spflux,[blacklineplot,lcplot],names,loc='upper right', fontsize=options.labelfontsize,labelspacing=options.legendlabelspacing,borderpad=options.legendborderpad,handletextpad=options.legendhandletextpad)
    

if __name__=='__main__':

    parser = mkplot_options()
    options,  filelist = parser.parse_args()

    matlib.clf()

    if options.figsize != None:
        print 'Changing figure size!'
        m = re.split(',',options.figsize)
        if len(m) == 1:
            xfigsize = yfigsize = float(m[0])
        elif  len(m) == 2:
            xfigsize = float(m[0])
            yfigsize = float(m[1])
        else:
            RuntimeError,"Error: %s must be a comma-separated list with maximum two entries!" % options.figsize
        matlib.figure(figsize=(xfigsize,yfigsize))
     
    matlib.subplots_adjust(wspace=0.0)
    matlib.subplots_adjust(hspace=0.0)

    spS = matlib.subplot(2,1,1)
    spSnorm = matlib.subplot(2,1,2,sharex=spS)

    (spSnorm,dummy1,dummy2) = mkplot(spSnorm,filelist,options)

    options.y = 'S_MWG_90,S_LMCavg_90,S_SMCbar_90'
    options.ylabel = '$S(\\lambda,\\theta)$'
    options.legend = 'MWG\,$\ \\theta=90^o$,LMC avg\,$\ \\theta=90^o$,SMC bar\,$\ \\theta=90^o$'
    options.colors = 'black'
    options.ymin=1E-24
    options.ymax=None

    options.major_ticks_y=1E-22
    yFormatter   = FormatStrFormatter('%.0e')
    spS.yaxis.set_major_formatter( yFormatter )
    
    (spS,dummy1,dummy2) = mkplot(spS,filelist,options)


    matlib.setp(spS.get_xticklabels(), visible=False)

    if options.savefile!='':
        print 'Saving file ',options.savefile
        matlib.savefig(options.savefile)
    matlib.show()

