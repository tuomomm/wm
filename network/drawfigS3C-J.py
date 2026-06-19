from pylab import *
import scipy.io
import scipy.stats

#-rw-r--r--@ 1 phtuma  staff  20060 Mar  6 12:33 ACC_psc_data.mat
#-rw-r--r--@ 1 phtuma  staff  20061 Mar  6 12:33 rMFG1_psc_data.mat

files = ['ACC_psc_data.mat','rMFG1_psc_data.mat']
areas = ['ACC','rMFG']
iareas = [1,0]
cols = ['#FF3333','#3333FF']
dimcols = ['#FFBBBB','#BBBBFF']
dimgray = '#CCCCCC'

def boxoff(ax,whichxoff='top'):
    ax.spines[whichxoff].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

f,axs = subplots(2,4)
axarr = axs.reshape(prod(axs.shape),).tolist()

for iax in range(0,4):
  for iay in range(0,2):
    axs[iay,iax].set_position([0.09+0.23*iax,0.8-0.27*iay,0.21,0.15])
    axs[iay,iax].tick_params(axis='both', which='major', labelsize=4)
    for axis in ['top','bottom','left','right']:
      axs[iay,iax].spines[axis].set_linewidth(0.3)
      boxoff(axs[iay,iax])

for itarget in [0,1]: #0=target, 1=non-target?
 for ifile in range(0,len(files)):
  iarea = iareas[ifile]
  A = scipy.io.loadmat(files[iarea])
  for iload in [0,1]:

    means_HC = []
    means_SCZ = []
    stds_HC = []
    stds_SCZ = []
    
    for itime in range(0,7):
      means_HC.append(mean(A['psc_g1'][itime][itarget][iload])*100)
      stds_HC.append(std(A['psc_g1'][itime][itarget][iload])*100)
      means_SCZ.append(mean(A['psc_g2'][itime][itarget][iload])*100)
      stds_SCZ.append(std(A['psc_g2'][itime][itarget][iload])*100)

    axs[iload,iarea+itarget*2].plot([-2.7+i*2.5 if i>0 else -2.5 for i in range(0,7)],means_HC,'k-',lw=0.6,label='HC')
    axs[iload,iarea+itarget*2].plot([-2.3+i*2.5 if i>0 else -2.5 for i in range(0,7)],means_SCZ,'-',color=cols[ifile],lw=0.6,label='SCZ')
    for itime in range(0,7):
      axs[iload,iarea+itarget*2].plot([-2.7+itime*2.5]*2, [means_HC[itime]-stds_HC[itime],means_HC[itime]+stds_HC[itime]], 'k-',lw=0.6)
      axs[iload,iarea+itarget*2].plot([-2.3+itime*2.5]*2, [means_SCZ[itime]-stds_SCZ[itime],means_SCZ[itime]+stds_SCZ[itime]], '-',color=cols[ifile],lw=0.6)

      pval = scipy.stats.ranksums(A['psc_g1'][itime][itarget][iload],A['psc_g2'][itime][itarget][iload])[1]
      if pval < 0.05/6:
        axs[iload,iarea+itarget*2].plot([-2.7+2.5*itime,-2.7+2.5*itime,-2.3+2.5*itime,-2.3+2.5*itime],[bracketys[iarea],bracketys[iarea]+3,bracketys[iarea]+3,bracketys[iarea]],'k-',lw=0.3)
        axs[iload,iarea+itarget*2].text(mean([-2.7+2.5*itime,-2.3+2.5*itime]),bracketys[iarea]+5,'**',fontsize=6,ha='center',va='center')
        print('Significance, itarget='+str(itarget)+', itime='+str(itime)+', pval='+str(pval))
      elif pval < 0.05:
        axs[iload,iarea+itarget*2].plot([-2.7+2.5*itime,-2.7+2.5*itime,-2.3+2.5*itime,-2.3+2.5*itime],[bracketys[iarea],bracketys[iarea]+3,bracketys[iarea]+3,bracketys[iarea]],'k-',lw=0.3)
        axs[iload,iarea+itarget*2].text(mean([-2.7+2.5*itime,-2.3+2.5*itime]),bracketys[iarea]+5,'*',fontsize=6,ha='center',va='center')
        print('Significance, itarget='+str(itarget)+', itime='+str(itime)+', pval='+str(pval))
    axs[iload,iarea+itarget*2].set_xlabel('Time (s)',fontsize=6)
    if iarea == itarget == 0:
      axs[iload,iarea+itarget*2].set_ylabel('Percent signal\nchange (%)',fontsize=6)
    else:
      axs[iload,iarea+itarget*2].set_yticklabels([])
    axs[iload,iarea+itarget*2].set_title(areas[iarea]+', '+str(iload)+' back, '+('target' if itarget==0 else 'nontarget'),fontsize=6)
    axs[iload,iarea+itarget*2].set_ylim([-65,65])
    axs[iload,iarea+itarget*2].legend(fontsize=5,labelspacing=0.25,columnspacing=1.25,framealpha=1.0)

    pos = axs[iload,iarea+itarget*2].get_position()
    f.text(pos.x0 - 0.04, pos.y1 - 0.0, chr(ord('C')+iarea+2*itarget+4*iload), fontsize=9)

      
  f.savefig("figS3C-J.pdf")
