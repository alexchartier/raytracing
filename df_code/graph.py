# module for graph settings in matplotlib 
import matplotlib as mpl
#_______________________________________________________________________________
def setGraphParameters(fontSize=17,typeFace='Helvetica',marker='.',markerSize=8,
                      lineWidth=2,lineStyle='-',cycleEnable=False):
   # set standard configuration for plot appearance
   # marker reduction for markers that are not a circle
   markerRed   = 6
   markerSize2 = markerSize - markerRed
   # line settings
   mpl.rcParams['lines.marker']     = marker
   mpl.rcParams['lines.markersize'] = markerSize
   mpl.rcParams['lines.linewidth']  = lineWidth
   mpl.rcParams['lines.linestyle']  = lineStyle
   if cycleEnable:
      mpl.rcParams['axes.prop_cycle']  = cycler( marker=['.','s','o','x','*','v','<','>'],
                                                 markersize=[markerSize,markerSize2,markerSize2,markerSize2,markerSize2,markerSize2,markerSize2,markerSize2],
                                                 color=['black','blue','red','green','black','blue','red','green'] )
   # font settings
   font = {'family': 'sans-serif', 'sans-serif': [typeFace], 'size': fontSize }
   mpl.rc('font',**font)
   # to avoid cell block limit...
   mpl.rcParams['agg.path.chunksize'] = 10000
   return 0
