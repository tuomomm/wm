from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection
from os.path import exists
import scipy.stats
import calcconds

from scipy.stats import linregress

addition = ''

def boxoff(ax,whichxoff='top'):
    ax.spines[whichxoff].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


def sci_notation(num, decimal_digits=1):
    """Format a float into LaTeX scientific notation string."""
    exponent = int(np.floor(np.log10(abs(num))))
    coeff = round(num / 10**exponent, decimal_digits)
    return f"{coeff}$\\times10^{{{exponent}}}$"
  
def drawdiscontinuity(ax,y,yoffset,x=0,xoffset=0.1,lw=2.0,lw2=1.0,whitecol='#FFFFFF',blackcol='#000000'):
  thisline = ax.plot([x-xoffset,x+xoffset],[y-yoffset,y],'-',linewidth=lw2,color=blackcol)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y,y+yoffset],'-',linewidth=lw2,color=blackcol)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-1.5*xoffset,x+1.5*xoffset],[y-0.75*yoffset,y+0.75*yoffset],'-',color=whitecol,zorder=100,linewidth=lw)
  thisline[0].set_clip_on(False)


blockeds_PFC = [['Calbin,CalbinC','0.985,0.985'], ['CK','0.959'], ['DAGK','0.998'], ['Gi','0.971'], ['Gqabg','1.063'], ['Ng','0.962'], ['PDE4','1.012'], ['PKA','1.048'], ['PKC','0.996'], ['comb','PFC']] 
blockeds_ACC = [['DAGK','1.039'], ['Gi','0.984'], ['NCX','0.999'], ['Ng','0.944'], ['PDE4','0.977'], ['PKA','1.039'], ['PKC','1.037'], ['PLA2','0.933'], ['PP1','0.966'], ['comb','ACC']] 

toPlot = ['PP1','PKA','NCX','PLA2','PDE4','PKC','CK','Ng','Gqabg','Calbin,CalbinC','DAGK','Gi','comb']
toPlotNames = ['PP1','PKA','NCX','PLA2','PDE4','PKC','CaMKII','Ng','Gq','Calbin','DAGK','Gi','Comb.']

f,axs = subplots(3,5)
axarr = axs.reshape(prod(axs.shape),).tolist()
axnew = []
axarr[13].set_visible(False)
axarr[14].set_visible(False)
for iax in range(0,13):
  iaxx = iax%4
  iaxy = int(iax/4)
  axarr[iax].set_position([0.56+0.115*iaxx,0.76-0.22*iaxy,0.08,0.155])
  axarr[iax].tick_params(axis='both', which='major', labelsize=4, length=2)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)
  boxoff(axarr[iax])

for iax in range(0,6):
  axnew.append(f.add_axes([0.09,0.8-0.15*iax,0.40,0.13]))
  axnew[iax].tick_params(axis='both', which='major', labelsize=4, length=2)
  for axis in ['top','bottom','left','right']:      
    axnew[iax].spines[axis].set_linewidth(0.3)
  boxoff(axnew[iax])
  for ipol in range(0,len(toPlot)):
    polygon = Polygon(array([[-1.5+5*ipol,-1.5+5*ipol,3.5+5*ipol,3.5+5*ipol],[0,100,100,0]]).T)
    p = PatchCollection([polygon], cmap=matplotlib.cm.jet, zorder=1)
    p.set_facecolor('#EEEEEE' if ipol%2 == 1 else '#FFFFFF')
    p.set_edgecolor(None)
    axnew[iax].add_collection(p)
  axnew[iax].set_xlim([-1.5,5*len(toPlot)-1.5])
  axnew[iax].set_xticks([])


def mybar(ax,x,y,facecolor=[],linewidth=0.3,w=0.4):
  qs = quantile(y, [0,0.25,0.5,0.75,1])
  polygon = Polygon(array([[x-w,x+w,x+w,x-w],[qs[1],qs[1],qs[3],qs[3]]]).T)
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  if type(facecolor) is not list or len(facecolor) > 0:
    p.set_facecolor(facecolor)
  p.set_edgecolor('#000000')
  p.set_linewidth(0.3)
  ax.add_collection(p)
  if isinf(qs[4]) or isnan(qs[4]):
    a2 = []
    a2.append(ax.plot([x-w,x+w,x,x,x-w,x+w,x,x],[qs[0],qs[0],qs[0],qs[2],qs[2],qs[2],qs[2],qs[3]+(qs[3]-qs[1])*0.5],'k-',lw=linewidth))
    a2.append(ax.plot([x,x],[qs[3]+(qs[3]-qs[1])*0.5,qs[3]+(qs[3]-qs[1])*1.0],'k--',lw=linewidth,dashes=(1.2,1.2)))
  else:
    a2 = ax.plot([x-w,x+w,x,x,x-w,x+w,x,x,x-w,x+w],[qs[0],qs[0],qs[0],qs[2],qs[2],qs[2],qs[2],qs[4],qs[4],qs[4]],'k-',lw=linewidth)
    print(str(qs[4]))
  return [p,a2]

tposts = [10,30*60]
stimSet = 'n1600_freq16.0'
    
cols = ['#000000','#3333FF','#FF3333']
dimcols = ['#CCCCCC','#BBBBFF','#FFBBBB']
dimgray = '#CCCCCC'
cols_compartments = ['#000000','#888888','#FF4444','#4444FF','#33FF33']

Nsum = 100

areas = ['ACC','PFC']

channels = ['CTRL', 'PConly']
channel_titles = ['CTRL','Combination']


FRThrCoeff = 2.0

filename_CTRL = '../synplast/nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,Gix0.5,0.5,1.5,1.5,0.984_onset24040000.0_n1600_freq16.0_dur3.0_flux150.0_Lflux5.0_Gluflux10.0_AChflux10.0.mat'
if exists(filename_CTRL):
  print('Loading '+filename_CTRL)
  conds_CTRL, times_CTRL = calcconds.calcconds_nrn(filename_CTRL)
  GluRconcs,GluRNs,times_GluR = calcconds.GluRSurfExpr(filename_CTRL)
  imin = argmin(conds_CTRL)
  imax = argmax(conds_CTRL)
  print('t_min = '+str(times_CTRL[imin]))
  print('t_max = '+str(times_CTRL[imax]))

for iax in range(0,len(toPlot)):
  axarr[iax].plot((times_CTRL-24040000)/1000/60,conds_CTRL/conds_CTRL[0],'-',lw=0.3,color=cols[0])
  for itoPlot in range(0,len(toPlot)):
    axnew[0].bar(0+5*itoPlot,conds_CTRL[0],facecolor=dimcols[0],width=0.6)
    axnew[1].bar(0+5*itoPlot,GluRconcs[2][0],facecolor=dimcols[0],width=0.6)
    
    axnew[2].bar(0+5*itoPlot,conds_CTRL[2],facecolor=dimcols[0],width=0.6)
    axnew[3].bar(0+5*itoPlot,conds_CTRL[1+30*6],facecolor=dimcols[0],width=0.6)

    axnew[4].bar(0+5*itoPlot,conds_CTRL[2]/conds_CTRL[0],facecolor=dimcols[0],width=0.6)
    axnew[5].bar(0+5*itoPlot,conds_CTRL[1+30*6]/conds_CTRL[0],facecolor=dimcols[0],width=0.6)

    drawdiscontinuity(axnew[0],33,0.1,0+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[0],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
    drawdiscontinuity(axnew[1],17,0.18,0+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[0],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
    drawdiscontinuity(axnew[2],27.4,0.06,0+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[0],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
    drawdiscontinuity(axnew[3],81.8,0.05,0+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[0],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
    drawdiscontinuity(axnew[4],0.795,0.0008,0+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[0],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
    drawdiscontinuity(axnew[5],2.15,0.0055,0+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[0],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
groupmeans_all = []
for iarea in range(0,len(areas)):
    blockeds = blockeds_PFC if iarea == 1 else blockeds_ACC
    printedLoading = False
    thisarea = areas[iarea]

    for iblocked in range(0,len(blockeds)):
      blocked = blockeds[iblocked][0]
      blockedCoeff = blockeds[iblocked][1]
      itoPlot = [i for i in range(0,len(toPlot)) if toPlot[i]==blocked]
      if len(itoPlot) == 0:
        itoPlot = len(blockeds)-1
        print(blocked+" not found in toPlot")
      else:
        itoPlot = itoPlot[0]
      filename = '../synplast/nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,'+blocked+'x0.5,0.5,1.5,1.5,'+blockedCoeff+'_onset24040000.0_n1600_freq16.0_dur3.0_flux150.0_Lflux5.0_Gluflux10.0_AChflux10.0.mat'
      if exists(filename):
        print('Loading '+filename)
        conds, times = calcconds.calcconds_nrn(filename)
        GluRconcs,GluRNs,times_GluR = calcconds.GluRSurfExpr(filename)
        axarr[itoPlot].plot((times-24040000)/1000/60,conds/conds[0],'-',lw=0.3,color=cols[1+iarea])
        axnew[0].bar(5*itoPlot+1+iarea,conds[0],facecolor=dimcols[1+iarea],width=0.6)
        print(areas[iarea]+' baseline '+str((conds[0]/conds_CTRL[0]-1)*100)+'% '+blocked+'x'+blockedCoeff+ '('+str(conds[0])+')')
        axnew[1].bar(5*itoPlot+1+iarea,GluRconcs[2][0],facecolor=dimcols[1+iarea],width=0.6)
        axnew[2].bar(5*itoPlot+1+iarea,conds[2],facecolor=dimcols[1+iarea],width=0.6)
        print(areas[iarea]+' abs-cond-10-sec '+str((conds[2]/conds_CTRL[2]-1)*100)+'% '+blocked+'x'+blockedCoeff+ '('+str(conds[2])+')')
        axnew[3].bar(5*itoPlot+1+iarea,conds[1+30*6],facecolor=dimcols[1+iarea],width=0.6)
        print(areas[iarea]+' abs-cond-30-min '+str((conds[1+30*6]/conds_CTRL[1+30*6]-1)*100)+'% '+blocked+'x'+blockedCoeff+ '('+str(conds[1+30*6])+')')
        axnew[4].bar(5*itoPlot+1+iarea,conds[2]/conds[0],facecolor=dimcols[1+iarea],width=0.6)
        print(areas[iarea]+' rel-cond-10-sec '+str(((conds[2]-conds[0])/(conds_CTRL[2]-conds_CTRL[0])-1)*100)+'% '+blocked+'x'+blockedCoeff+' ('+str((conds[2]/conds[0]-1)*100)+'% from '+str((conds_CTRL[2]/conds_CTRL[0]-1)*100)+'% in CTRL)')
        axnew[5].bar(5*itoPlot+1+iarea,conds[1+30*6]/conds[0],facecolor=dimcols[1+iarea],width=0.6)
        print(areas[iarea]+' rel-cond-30-min '+str(((conds[1+30*6]-conds[0])/(conds_CTRL[1+30*6]-conds_CTRL[0])-1)*100)+'% '+blocked+'x'+blockedCoeff+' ('+str((conds[1+30*6]/conds[0]-1)*100)+'% from '+str((conds_CTRL[1+30*6]/conds_CTRL[0]-1)*100)+'% in CTRL)')

        drawdiscontinuity(axnew[0],33-0.1*(1+iarea),0.1,1+iarea+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[1+iarea],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
        drawdiscontinuity(axnew[1],17-0.18*(1+iarea),0.18,1+iarea+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[1+iarea],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
        drawdiscontinuity(axnew[2],27.4-0.06*(1+iarea),0.06,1+iarea+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[1+iarea],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
        drawdiscontinuity(axnew[3],81.8-0.05*(1+iarea),0.05,1+iarea+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[1+iarea],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
        drawdiscontinuity(axnew[4],0.795-0.0008*(1+iarea),0.0008,1+iarea+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[1+iarea],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
        drawdiscontinuity(axnew[5],2.15-0.0055*(1+iarea),0.0055,1+iarea+5*itoPlot,0.32,1.0,0.8,blackcol=dimcols[1+iarea],whitecol='#EEEEEE' if itoPlot%2==1 else '#FFFFFF')
      else:
        print(filename+' does not exist')

axnew[0].set_ylabel('Baseline\ncond. (pS)',fontsize=6)
axnew[1].set_ylabel('Baseline\n[GluR1] at\nmemb. (nM)',fontsize=6)
axnew[2].set_ylabel('Cond. at\n10 s (pS)',fontsize=6)
axnew[3].set_ylabel('Cond. at\n30 min (pS)',fontsize=6)
axnew[4].set_ylabel('Rel. cond.\nat 10 s\n(fold change)',fontsize=6)
axnew[5].set_ylabel('Rel. cond.\nat 30 min\n(fold change)',fontsize=6)

        
axnew[0].set_ylim([32,39])
axnew[1].set_ylim([15,30])
axnew[2].set_ylim([27,31])
axnew[3].set_ylim([81.6,84.5])
axnew[4].set_ylim([0.79,0.83])
axnew[5].set_ylim([2.1,2.4])

axnew[0].set_yticks([34,36,38])
axnew[1].set_yticks([20,25,30])
axnew[2].set_yticks([28,29,30,31])
axnew[3].set_yticks([82,83,84])
axnew[4].set_yticks([0.8,0.81,0.82,0.83])
axnew[5].set_yticks([2.2,2.3,2.4])

drawdiscontinuity(axnew[0],33,0.1,-1.5,0.32,1.0,0.8)
drawdiscontinuity(axnew[1],17,0.18,-1.5,0.32,1.0,0.8,blackcol=dimcols[0])
drawdiscontinuity(axnew[2],27.4,0.06,-1.5,0.32,1.0,0.8,blackcol=dimcols[0])
drawdiscontinuity(axnew[3],81.8,0.05,-1.5,0.32,1.0,0.8,blackcol=dimcols[0])
drawdiscontinuity(axnew[4],0.795,0.0008,-1.5,0.32,1.0,0.8,blackcol=dimcols[0])
drawdiscontinuity(axnew[5],2.15,0.0055,-1.5,0.32,1.0,0.8,blackcol=dimcols[0])

for ipol in range(0,len(toPlot)):
  axnew[0].text(1+5*ipol,38,toPlotNames[ipol],ha='center',va='center',fontsize=5,rotation=45)
        
f.savefig("figS6.pdf")

        
axarr[0].set_ylabel('Rel. cond. (A.U.)',fontsize=6)
axarr[4].set_ylabel('Rel. cond. (A.U.)',fontsize=6)
axarr[8].set_ylabel('Rel. cond. (A.U.)',fontsize=6)
axarr[12].set_ylabel('Rel. cond. (A.U.)',fontsize=6)
for iax in range(0,13):
    axarr[iax].set_title(toPlotNames[iax],fontsize=6)
for iax in range(10,15):
    axarr[iax].set_xlabel('$t$ (min)',fontsize=6)
    
for iax in range(0,6):
    pos = axnew[iax].get_position()
    f.text(pos.x0 - 0.08, pos.y1 - 0.0, chr(ord('A')+iax), fontsize=9)
for iax in range(0,13):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.035, pos.y1 - 0.0, chr(ord('G')+iax), fontsize=9)

    
f.savefig("figS6.pdf")
