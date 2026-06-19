from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection
from os.path import exists
import scipy.stats
import calcconds

from scipy.stats import linregress

addition = ''

def sci_notation(num, decimal_digits=1):
    """Format a float into LaTeX scientific notation string."""
    exponent = int(np.floor(np.log10(abs(num))))
    coeff = round(num / 10**exponent, decimal_digits)
    return f"{coeff}$\\times10^{{{exponent}}}$"
  
def drawdiscontinuity(ax,y,yoffset,x=0,xoffset=0.1,lw=2.0,lw2=1.0):
  thisline = ax.plot([x-xoffset,x+xoffset],[y-yoffset,y],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y,y+yoffset],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-1.5*xoffset,x+1.5*xoffset],[y-0.75*yoffset,y+0.75*yoffset],'k-',color='#FFFFFF',zorder=100,linewidth=lw)
  thisline[0].set_clip_on(False)

f,axarr = subplots(1,6)
axnew = []
for iay in range(0,2):
  axarr[iay*3+0].set_position([0.08+0.27*0,0.75-0.32*iay,0.21,0.16])
  axarr[iay*3+1].set_position([0.08+0.27*1-0.04,0.75-0.32*iay,0.11,0.16])
  axarr[iay*3+2].set_position([0.08+0.27*1+0.10,0.75-0.32*iay,0.11,0.16])
axarr[4].set_position([0.08+0.27*1,0.75-0.32*1,0.21,0.16])
axarr[5].set_visible(False)
for iax in range(0,6):
  axarr[iax].tick_params(axis='both', which='major', labelsize=4)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)

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

stimSet = 'n1600_freq16.0'
    
cols = ['#3333FF','#FF3333']
dimcols = ['#BBBBFF','#FFBBBB']
dimgray = '#CCCCCC'


Nsum = 100
Ncellsperpop = [582,97,97,291,291]

areas = ['ACC','PFC']
channels = ['CTRL', 'PConly']
channel_titles = ['CTRL','Combination']

isubjs_HC_PFC = [0, 1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20, 22, 23, 25, 27, 28, 31, 32, 33, 34, 35, 36, 37, 39, 41, 42, 43, 47, 51, 52, 53, 54, 63, 64, 66, 69, 73, 80, 81, 82, 83, 84, 87, 88, 89, 90, 92, 93, 94, 96, 97, 100, 105, 107, 111, 112, 118, 119, 120, 125, 127, 128, 135, 139, 140, 142, 143, 145, 146, 152, 154, 156, 165, 168, 176, 181, 182, 187, 192, 194, 201, 203, 204, 205, 207, 210, 211, 217, 220, 222, 223, 225, 226, 232, 235, 236, 240, 241, 242, 243, 244, 247, 248, 249, 251, 252, 253, 254, 257, 258, 262, 263, 265, 266, 267, 268, 269, 270, 271, 272, 283, 290, 291, 293, 294, 295, 296, 297, 298, 299, 300, 301, 305, 306, 308, 313, 314, 315, 317, 319, 320, 321, 323, 324, 325, 326, 329, 330, 332, 333, 335, 336, 339, 340, 341, 342, 343, 344, 345, 347, 351, 352, 354, 355, 356, 359, 361, 362, 368, 371, 373, 374, 375, 376, 377, 380, 381, 382, 383, 384, 386, 387, 388, 389, 390, 393, 396, 397, 398, 399, 400, 402, 403, 405, 406, 407, 408, 415, 416, 417, 420, 421, 422, 423, 424, 425, 426, 427, 428]
isubjs_SCZ_PFC = [4, 6, 15, 17, 21, 24, 26, 29, 30, 38, 40, 44, 45, 46, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 62, 65, 67, 68, 70, 71, 72, 74, 75, 76, 77, 78, 79, 85, 86, 91, 95, 98, 99, 101, 102, 103, 104, 106, 108, 109, 110, 113, 114, 115, 116, 117, 121, 122, 123, 124, 126, 129, 130, 131, 132, 133, 134, 136, 137, 138, 141, 144, 147, 148, 149, 150, 151, 153, 155, 157, 158, 159, 160, 161, 162, 163, 164, 166, 167, 169, 170, 171, 172, 173, 174, 175, 177, 178, 179, 180, 183, 184, 185, 186, 188, 189, 190, 191, 193, 195, 196, 197, 198, 199, 200, 202, 206, 208, 209, 212, 213, 214, 215, 216, 218, 219, 221, 224, 227, 228, 229, 230, 231, 233, 234, 237, 238, 239, 245, 246, 250, 255, 256, 259, 260, 261, 264, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 285, 286, 287, 288, 289, 292, 302, 303, 304, 307, 309, 310, 311, 312, 316, 318, 322, 327, 328, 334, 337, 338, 346, 348, 349, 350, 353, 357, 358, 360, 363, 364, 365, 366, 367, 369, 370, 372, 378, 379, 385, 392, 394, 395, 401, 404, 409, 410, 411, 412, 413, 414, 418, 419]
isubjs_HC_ACC = [1, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 73, 75, 76, 77, 78, 79, 80, 83, 84, 86, 90, 91, 92, 95, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 124, 126, 127, 129, 130, 131, 132, 133, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152, 154, 159, 161, 163, 164, 167, 169, 170, 172, 174, 175, 176, 177, 182, 199, 208, 235, 240, 243, 249, 253, 256, 257, 260, 261, 262, 263, 265, 266, 267, 268, 271, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 329, 330, 331, 332, 343, 349, 350, 352, 353, 355, 360, 366, 369, 370, 371, 372, 374, 375, 376, 378, 381, 382, 383, 385, 388, 389, 392, 393, 395, 396, 397, 398, 401, 402, 409, 410, 414, 416, 418, 421, 422, 426, 427, 431, 432, 433, 434, 435, 437, 438, 442, 443, 444, 446, 447, 448, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480]
isubjs_SCZ_ACC = [0, 2, 3, 12, 17, 28, 31, 58, 69, 70, 71, 72, 74, 81, 82, 85, 87, 88, 89, 93, 94, 96, 97, 98, 99, 100, 116, 117, 118, 119, 120, 121, 122, 123, 125, 134, 151, 155, 156, 157, 158, 160, 162, 165, 166, 168, 171, 173, 178, 179, 180, 181, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 200, 201, 202, 203, 204, 205, 206, 207, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 236, 237, 238, 239, 241, 242, 244, 245, 246, 247, 248, 250, 251, 252, 254, 255, 258, 259, 264, 269, 270, 272, 273, 274, 275, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 344, 345, 346, 347, 348, 351, 354, 356, 357, 358, 359, 361, 362, 363, 364, 365, 367, 368, 373, 377, 379, 380, 384, 386, 387, 390, 391, 394, 399, 400, 403, 404, 405, 406, 407, 411, 412, 413, 415, 417, 419, 420, 423, 424, 425, 428, 429, 430, 436, 439, 440, 441, 445, 449, 450, 451, 452]

FRThrCoeff = 2.0

iarea = 0
if len(sys.argv) > 1:
  try:
    iarea = int(sys.argv[1])
  except:
    if sys.argv[1] in ['PFC','PFCmild']:
      iarea = 1

ipops = [3]
iipop = 0

tposts =[0, 0.1666, 30]
tposts_str = ['0','10 sec','30 min']
for itpost in range(0,3):
  tpost = tposts[itpost]
  ipop = ipops[iipop]

  groupmeans_all = []
  thisarea = areas[iarea]+'mild'
    
  #Subject-wise:
  isubjs_HC = isubjs_HC_PFC if iarea == 1 else isubjs_HC_ACC
  isubjs_SCZ = isubjs_SCZ_PFC if iarea == 1 else isubjs_SCZ_ACC
  groupmeans = []
  FRthisarea = []
  condsthisarea = []
  for igroup in [0,1]:
      isubjs = isubjs_HC if igroup == 0 else isubjs_SCZ
      FRvecs_all_means = []
      FRthisgroup = []
      condsthisgroup = []
      for iisubj in range(0,len(isubjs)):
        isOK = 1
        isubj = isubjs[iisubj]
        bins_all = []
        FRvecs = []
        myseed = 1
        filename = '../synplast/nrn_tstop27000000_tol1e-06_GluR1,GluR1_memb,GluR2,GluR2_memb,CMcomb_samps_nonIC_'+areas[iarea]+'x0.5,0.5,1.5,1.5,'+str(isubj)+'_onset24040000.0_'+stimSet+'_dur3.0_flux150.0_Lflux5.0_Gluflux10.0_AChflux10.0.mat'
        if exists(filename):
          print("Loading "+filename)
          A = scipy.io.loadmat(filename)
          conds, times = calcconds.calcconds_nrn(filename)
          minabs = min(abs(times-(24040000+tpost*60*1000)))
          its = [it for it in range(0,len(times)) if times[it]-(24040000+tpost*60*1000) == minabs][0]
        else:
          print(filename+' does not exist')
          isOK = 0
            
        filename = 'subjs/FRandMUA_sim_net_long_withSK3_'+thisarea+'subj'+str(isubj)+addition+'_altgain_PConly_seed'+str(myseed)+'.mat'
        if not exists(filename) and exists(filename.replace('subjs','subjs/oldCMcomb')):
          print('Using '+filename.replace('subjs','subjs/oldCMcomb')+' instead of '+filename)
          filename = filename.replace('subjs','subjs/oldCMcomb')
        elif not exists(filename) and exists(filename.replace('FRandMUA','FR')):
          print('Using '+filename.replace('FRandMUA','FR')+' instead of '+filename)
          filename = filename.replace('FRandMUA','FR')
        if exists(filename):
          print("Loading "+filename)
          A = scipy.io.loadmat(filename)
          FRs = zeros(int(len(A['FRvecs'][ipop])/Nsum)+1,)
          for ibin in range(0,len(FRs)):
            FRs[ibin] = sum(A['FRvecs'][ipop][Nsum*ibin:Nsum*(ibin+1)])/(Ncellsperpop[ipop]*Nsum/100) #Divided by the number of cells in the population, corrected by the bin size (if bin size is 1 s, i.e. Nsum=100, then we have the frequency in Hz without normalization, but if 10 s, i.e. Nsum=1000, then we have to also divide by 10 to get spikes/s)  
        elif iisubj == 0:
          print(filename+' does not exist')
          isOK = 0
        if isOK: #Add the data only if both data were available
          if itpost == 0:
            condsthisgroup.append(conds[0])
          else:
            condsthisgroup.append(conds[its]/conds[0])
          FRthisgroup.append(FRs[:])
        else:
          print('isubj '+str(isubj)+' not OK')
          qwe
      FRthisarea.append(FRthisgroup[:])
      groupmeans.append(mean(array(FRthisgroup),axis=0))
      condsthisarea.append(condsthisgroup[:])
    
  FRbeyondThrthisarea = []
  FRbaseline = mean(groupmeans[0][3:10])
  groupmeans_all.append(groupmeans[:])
  for igroup in [0,1]:
      isubjs = isubjs_HC if igroup == 0 else isubjs_SCZ
      FRthisgroup = FRthisarea[igroup]
      FRbeyondThrthisgroup = []
      for iisubj in range(0,len(FRthisgroup)):
        FRbeyondThrthissubj = []
        FRmean = FRthisgroup[iisubj]

        FRThr = FRThrCoeff*FRbaseline
        FRbeyondThr = sum([1 for i in range(0,len(FRmean)) if FRmean[i] >= FRThr])
        FRbeyondThrthissubj.append(FRbeyondThr)
        FRbeyondThrthisgroup.append(FRbeyondThrthissubj[:])
      FRbeyondThrthisarea.append(FRbeyondThrthisgroup[:])

  for iigroup in [0,1,2]:
      axarrind = iigroup if itpost == 0 else 2+itpost

      igroup = iigroup-1
      if igroup >= 0:
        axarr[axarrind].plot([x[0] for x in FRbeyondThrthisarea[igroup]],condsthisarea[igroup],'b.',ms=1.5,color=cols[iarea] if igroup else '#000000',label='R='+'{:.2f}'.format(r_value)+', p='+sci_notation(p_value))
        print('Plotting axarrind='+str(axarrind)+' conds = '+str(mean(condsthisarea[igroup]))+' igroup = '+str(igroup)+' itpost = '+str(itpost))
        if itpost == 0:
          axarr[0].plot([x[0] for x in FRbeyondThrthisarea[igroup]],condsthisarea[igroup],'b.',ms=1.5,color=cols[iarea] if igroup else '#000000',label='R='+'{:.2f}'.format(r_value)+', p='+sci_notation(p_value))
          print('Plotting axarrind='+str(axarrind)+' conds = '+str(mean(condsthisarea[igroup]))+' igroup = '+str(igroup)+' itpost = '+str(itpost)+' again in '+str(0))

      Xdata = [x[0] for x in FRbeyondThrthisarea[igroup]] if igroup >= 0 else [x[0] for x in FRbeyondThrthisarea[0]]+[x[0] for x in FRbeyondThrthisarea[1]]
      Ydata = condsthisarea[igroup] if igroup >= 0 else condsthisarea[0]+condsthisarea[1]
      
      slope, intercept, r_value, p_value, std_err = linregress(Xdata,Ydata)
      xlim = [0,100]

      if itpost:
        ylim = [1.7*(tpost>=5),2.7] if tpost >= 1 else [0.7,0.92]
      else:
        ylim = [29,50]

      axarr[axarrind].set_xlim([xlim[0],xlim[1]])
      print(areas[iarea]+' igroup='+str(igroup)+' P-val = '+str(p_value)+", R="+str(r_value))

      if iigroup == 0 or itpost == 0:
        axarr[axarrind].plot([xlim[0],xlim[1]],[intercept+slope*xlim[0],intercept+slope*xlim[1]],'-',color='#888888',lw=0.6) #'#000000' if iigroup==0 else cols[iarea],lw=0.5)
        print('Plotting regression axarrind='+str(axarrind))
        axarr[axarrind].set_title('R='+'{:.2f}'.format(r_value)+(', ' if iigroup==0 else '\n')+'p='+sci_notation(p_value),color=('#191980' if iarea==0 else '#801919') if igroup < 0 else ('#000000' if igroup == 0 else cols[iarea]),fontsize=6)

  print(areas[iarea]+' p-val FRdur HC vs SCZ = '+str(scipy.stats.ranksums([x[0] for x in FRbeyondThrthisarea[0]],[x[0] for x in FRbeyondThrthisarea[1]])[1])+' N='+str(len(FRbeyondThrthisarea[0]))+','+str(len(FRbeyondThrthisarea[1]))+', means='+str(mean([x[0] for x in FRbeyondThrthisarea[0]]))+','+str(mean([x[0] for x in FRbeyondThrthisarea[1]])))
  print(areas[iarea]+' p-val postHFScond HC vs SCZ = '+str(scipy.stats.ranksums(condsthisarea[0],condsthisarea[1])[1])+' N='+str(len(condsthisarea[0]))+','+str(len(condsthisarea[1]))+', means='+str(mean(condsthisarea[0]))+','+str(mean(condsthisarea[1])))

  for iax in range(0,3):
    axarr[iax].set_xlabel('High-activity\nduration (s)\n'+('Both' if iax%3==0 else ('SCZ' if iax%3 == 2 else 'HC')),fontsize=6)
    if iax%3==0:
      axarr[iax].set_ylabel(areas[iarea]+'\nBaseline cond',fontsize=6)
    else:
      axarr[iax].set_yticklabels([])
  for iax in range(0,5):
    pos = axarr[iax].get_position()
    f.text(pos.x0 - 0.02, pos.y1 + 0.015 + 0.02*(iax<3), chr(ord('A')+iax+2*(iarea==0)), fontsize=9)
  axarr[3].set_ylabel('Rel. post-HFS cond\n'+areas[iarea]+',$t$='+tposts_str[1],fontsize=6)
  axarr[4].set_ylabel(areas[iarea]+',$t$='+tposts_str[2],fontsize=6)
  for iax in [0,3,4]:
    axarr[iax].set_xlabel('High-activity duration (s)\nBoth HC and SCZ',fontsize=6)


if iarea == 0:
  f.savefig("fig4C-G.pdf")
else:
  f.savefig("figS7.pdf")
    

