from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection
from os.path import exists
import scipy.stats

f,axs = subplots(2,8)
axarr = axs.reshape(prod(axs.shape),).tolist()

for iay in range(0,2):
  for iax in range(0,5):
    axarr[iay*5+iax].set_position([0.1+0.18*iax,0.78-0.26*iay,0.12,0.14])
for iax in range(0,15):
    axarr[iax].tick_params(axis='both', which='major', labelsize=4)
    for axis in ['top','bottom','left','right']:
      axarr[iax].spines[axis].set_linewidth(0.3)

axarr[10].set_position([0.1,0.06,0.24,0.36])
axarr[11].set_position([0.4,0.06,0.24,0.36])
axarr[12].set_position([0.7,0.06,0.24,0.36])
axarr[13].set_position([0.4+0.19,0.06+0.2,0.04,0.15])
axarr[14].set_position([0.7+0.19,0.06+0.2,0.04,0.15])
axarr[15].set_visible(False)

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


cols = ['#FF3333','#3333FF']
dimcols = ['#FFBBBB','#BBBBFF']
dimgray = '#CCCCCC'
cols_compartments = ['#000000','#888888','#FF4444','#4444FF','#33FF33']

areas = ['PFC','ACC']
iareas = [1,0]

channels = ['CTRL', 'nax_PConly', 'iar_PConly', 'kap', 'cat', 'kdr_PConly', 'serca', 'cagk_ikc', 'km', 'cal_PConly', 'GABAB', 'PConly']
channel_titles = ['CTRL','Fast Na$^+$', 'HCN', 'A-type K$^+$', 'T-type Ca$^{2+}$', 'DR K$^+$', 'SERCA', 'BK', 'M-type K$^+$', 'L-type Ca$^{2+}$', 'GABAB', 'Combination']
channel_genes = ['','\n(SCN1B)','\n(HCN1)','\n(KCND3)','\n(CACNA1I)','\n(KCNB1)','\n(ATP2A2)','\n(KCNMA1)','\n(KCNQ3)','\n(CACNA1C,CACNA1D)','\n(KCNJ6)','']

isubjs_HC_PFC = [0, 1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20, 22, 23, 25, 27, 28, 31, 32, 33, 34, 35, 36, 37, 39, 41, 42, 43, 47, 51, 52, 53, 54, 63, 64, 66, 69, 73, 80, 81, 82, 83, 84, 87, 88, 89, 90, 92, 93, 94, 96, 97, 100, 105, 107, 111, 112, 118, 119, 120, 125, 127, 128, 135, 139, 140, 142, 143, 145, 146, 152, 154, 156, 165, 168, 176, 181, 182, 187, 192, 194, 201, 203, 204, 205, 207, 210, 211, 217, 220, 222, 223, 225, 226, 232, 235, 236, 240, 241, 242, 243, 244, 247, 248, 249, 251, 252, 253, 254, 257, 258, 262, 263, 265, 266, 267, 268, 269, 270, 271, 272, 283, 290, 291, 293, 294, 295, 296, 297, 298, 299, 300, 301, 305, 306, 308, 313, 314, 315, 317, 319, 320, 321, 323, 324, 325, 326, 329, 330, 332, 333, 335, 336, 339, 340, 341, 342, 343, 344, 345, 347, 351, 352, 354, 355, 356, 359, 361, 362, 368, 371, 373, 374, 375, 376, 377, 380, 381, 382, 383, 384, 386, 387, 388, 389, 390, 393, 396, 397, 398, 399, 400, 402, 403, 405, 406, 407, 408, 415, 416, 417, 420, 421, 422, 423, 424, 425, 426, 427, 428]
isubjs_SCZ_PFC = [4, 6, 15, 17, 21, 24, 26, 29, 30, 38, 40, 44, 45, 46, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 62, 65, 67, 68, 70, 71, 72, 74, 75, 76, 77, 78, 79, 85, 86, 91, 95, 98, 99, 101, 102, 103, 104, 106, 108, 109, 110, 113, 114, 115, 116, 117, 121, 122, 123, 124, 126, 129, 130, 131, 132, 133, 134, 136, 137, 138, 141, 144, 147, 148, 149, 150, 151, 153, 155, 157, 158, 159, 160, 161, 162, 163, 164, 166, 167, 169, 170, 171, 172, 173, 174, 175, 177, 178, 179, 180, 183, 184, 185, 186, 188, 189, 190, 191, 193, 195, 196, 197, 198, 199, 200, 202, 206, 208, 209, 212, 213, 214, 215, 216, 218, 219, 221, 224, 227, 228, 229, 230, 231, 233, 234, 237, 238, 239, 245, 246, 250, 255, 256, 259, 260, 261, 264, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 285, 286, 287, 288, 289, 292, 302, 303, 304, 307, 309, 310, 311, 312, 316, 318, 322, 327, 328, 334, 337, 338, 346, 348, 349, 350, 353, 357, 358, 360, 363, 364, 365, 366, 367, 369, 370, 372, 378, 379, 385, 392, 394, 395, 401, 404, 409, 410, 411, 412, 413, 414, 418, 419]
isubjs_HC_ACC = [1, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 73, 75, 76, 77, 78, 79, 80, 83, 84, 86, 90, 91, 92, 95, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 124, 126, 127, 129, 130, 131, 132, 133, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152, 154, 159, 161, 163, 164, 167, 169, 170, 172, 174, 175, 176, 177, 182, 199, 208, 235, 240, 243, 249, 253, 256, 257, 260, 261, 262, 263, 265, 266, 267, 268, 271, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 329, 330, 331, 332, 343, 349, 350, 352, 353, 355, 360, 366, 369, 370, 371, 372, 374, 375, 376, 378, 381, 382, 383, 385, 388, 389, 392, 393, 395, 396, 397, 398, 401, 402, 409, 410, 414, 416, 418, 421, 422, 426, 427, 431, 432, 433, 434, 435, 437, 438, 442, 443, 444, 446, 447, 448, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480]
isubjs_SCZ_ACC = [0, 2, 3, 12, 17, 28, 31, 58, 69, 70, 71, 72, 74, 81, 82, 85, 87, 88, 89, 93, 94, 96, 97, 98, 99, 100, 116, 117, 118, 119, 120, 121, 122, 123, 125, 134, 151, 155, 156, 157, 158, 160, 162, 165, 166, 168, 171, 173, 178, 179, 180, 181, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 200, 201, 202, 203, 204, 205, 206, 207, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 236, 237, 238, 239, 241, 242, 244, 245, 246, 247, 248, 250, 251, 252, 254, 255, 258, 259, 264, 269, 270, 272, 273, 274, 275, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 344, 345, 346, 347, 348, 351, 354, 356, 357, 358, 359, 361, 362, 363, 364, 365, 367, 368, 373, 377, 379, 380, 384, 386, 387, 390, 391, 394, 399, 400, 403, 404, 405, 406, 407, 411, 412, 413, 415, 417, 419, 420, 423, 424, 425, 428, 429, 430, 436, 439, 440, 441, 445, 449, 450, 451, 452]

FRThrCoeff = 2.0

for iiarea in range(0,len(iareas)):
  iarea = iareas[iiarea]
  thisarea_temp = areas[iarea]
  bins_all_all = []
  for ichannel in range(0,len(channels)):
    thischannel = channels[ichannel]
    thisarea = thisarea_temp
    if thischannel in ["PConly","kap","nax_PConly"]:
      thisarea = thisarea + "mild"
      print("channel = "+thischannel+", thisarea = "+thisarea)
    if 'PFC' in thisarea and thischannel in ['kap','cagk_ikc'] or 'ACC' in thisarea and thischannel in ['iar_PConly','nax_PConly']:
      continue
    bins_all = []
    for iseed in range(1,6):
      filename = 'sim_onepyr_long_withSK3_'+thisarea+'_'+thischannel+'_altgain_seed'+str(iseed)+'.mat'
      if thischannel == '':
        filename = 'sim_onepyr_long_withSK3_'+thisarea+'_altgain_seed'+str(iseed)+'.mat'
      elif thischannel == 'CTRL':
        filename = 'sim_onepyr_long_withSK3_CTRL_altgain_seed'+str(iseed)+'.mat'

      binwidth = 1000

      if not exists(filename):
        print(filename+' does not exist')
        continue
      print('Loading '+filename)
      A = scipy.io.loadmat(filename)

      bins = zeros(int(100000/binwidth)+1,)
      for ispike in range(0,len(A['spikes'][0])):
        bins[int(A['spikes'][0][ispike]/binwidth)] = bins[int(A['spikes'][0][ispike]/binwidth)]+1
      bins_all.append(bins[:])

    FRmean = mean(array(bins_all),axis=0)/(binwidth/1000)
    if thischannel == 'CTRL':
      FRbaseline = mean(FRmean[3:10])
    if ichannel > 0:
      FRstd = std(array(bins_all),axis=0)/(binwidth/1000)
      polygon = Polygon(array([(binwidth*array(list(range(0,int(100000/binwidth)+1))+list(range(int(100000/binwidth),-1,-1)))-10000)/1000,
                               r_[FRmean-FRstd,FRmean[::-1]+FRstd[::-1]]]).T)
      p = PatchCollection([polygon], cmap=matplotlib.cm.jet, zorder=1)
      p.set_facecolor(dimcols[iarea])
      p.set_edgecolor(None)
      axarr[ichannel-1].add_collection(p)
      axarr[ichannel-1].plot((binwidth*array(range(0,int(100000/binwidth)+1))-10000)/1000,FRmean,lw=0.1,color=cols[iarea],label=thisarea)

      axarr[ichannel-1].plot((binwidth*array(range(0,int(100000/binwidth)+1))-10000)/1000,mean(array(bins_all_all[0]),axis=0)/(binwidth/1000),lw=0.1,color='#000000',label='CTRL')
    bins_all_all.append(bins_all[:])

    if 'Comb' in channel_titles[ichannel]:
      axarr[ichannel-1].set_title(channel_titles[ichannel]+channel_genes[ichannel],fontsize=5)
    else:
      axarr[ichannel-1].set_title(channel_titles[ichannel]+channel_genes[ichannel],fontsize=5,pad=2.5)
    

  isubjs_HC = isubjs_HC_PFC if iarea == 0 else isubjs_HC_ACC
  isubjs_SCZ = isubjs_SCZ_PFC if iarea == 0 else isubjs_SCZ_ACC
  groupmeans = []
  
  bins_thisarea = []
  FRbeyondThrthisarea = []
  for igroup in [0,1]:
    isubjs = isubjs_HC if igroup == 0 else isubjs_SCZ
    bins_all_means = []
    bins_thisgroup = []
    FRvecs_all_means = []
    for iisubj in range(0,len(isubjs)):
      isubj = isubjs[iisubj]
      bins_all = []
      for myseed in [1]:
        filename = 'subjs/sim_onepyr_long_withSK3_'+thisarea+'subj'+str(isubj)+'_seed'+str(myseed)+'.mat'
        if not exists(filename):
          print(filename+"does not exist")
        if exists(filename):
          if myseed == 1:
            print('Loading '+filename)
          A = scipy.io.loadmat(filename)
          bins = zeros(int(100000/binwidth)+1,)
          for ispike in range(0,len(A['spikes'][0])):
            bins[int(A['spikes'][0][ispike]/binwidth)] = bins[int(A['spikes'][0][ispike]/binwidth)]+1
          bins_all.append(bins[:])
        else:
          print(filename+' does not exist')
      
      if len(bins_all) > 0:
        FRmean = mean(array(bins_all),axis=0)/(binwidth/1000)
        axarr[11+iiarea].plot((binwidth*array(range(0,int(100000/binwidth)+1))-10000)/1000,FRmean,lw=0.1,color=dimcols[iarea] if igroup==1 else dimgray,label=thisarea)
        FRvecs_all_means.append([mean([bins_all[isamp][it] for isamp in range(0,len(bins_all))],axis=0)/(binwidth/1000) for it in range(0,len(bins_all[0]))])
      else:
        FRvecs_all_means.append(nan+zeros(int(100000/binwidth)+1,))
        
      bins_thisgroup.append(bins_all[:])
    bins_thisarea.append(bins_thisgroup[:])
    groupmeans.append(mean(array(FRvecs_all_means),axis=0))

  FRbaseline = mean(groupmeans[0][3:10])
  FRThr = FRThrCoeff*FRbaseline
  print("FRThr = "+str(FRThr))

  for igroup in [0,1]:
    isubjs = isubjs_HC if igroup == 0 else isubjs_SCZ
    bins_thisgroup = bins_thisarea[igroup]
    FRbeyondThrthisgroup = []
    for iisubj in range(0,len(isubjs)):
      isubj = isubjs[iisubj]
      bins_all = bins_thisgroup[iisubj]
      if len(bins_all) > 0:
        FRmean = mean(array(bins_all),axis=0)/(binwidth/1000)
        FRbeyondThr = sum([1 for i in range(0,len(FRmean)) if FRmean[i] >= FRThr])
        FRbeyondThrthisgroup.append(FRbeyondThr)
      else:
        FRbeyondThrthisgroup.append(nan)

        
      bins_all_means.append(FRmean[:])
    print(thisarea+' '+str(igroup)+' done')
    FRbeyondThrthisarea.append(FRbeyondThrthisgroup[:])
  for igroup in [0,1]:
    axarr[11+iiarea].plot((binwidth*array(range(0,int(100000/binwidth)+1))-10000)/1000,groupmeans[igroup],lw=0.75,color=cols[iarea] if igroup else '#000000',label=thisarea)
    mybar(axarr[13+iiarea],igroup,FRbeyondThrthisarea[igroup],facecolor=cols[iarea] if igroup else '#888888')
    
  pval = scipy.stats.ranksums(FRbeyondThrthisarea[0],FRbeyondThrthisarea[1])[1]
  print(thisarea+' : '+str(median(FRbeyondThrthisarea[0]))+' +- '+str(std(FRbeyondThrthisarea[0]))+' vs '+str(median(FRbeyondThrthisarea[1]))+' +- '+str(std(FRbeyondThrthisarea[1]))+', pval = '+str(pval))
  if pval < 0.05:
    axarr[13+iiarea].plot([0,0,1,1],[103,108,108,103],'k-',lw=0.1)
    axarr[13+iiarea].plot(0.5,110,'k*',mew=0.9,lw=0.9,ms=3)    
  print("p-val = "+str(pval))
  axarr[11+iiarea].set_title('Subject-wise, '+areas[iarea],fontsize=5)
  axarr[13+iiarea].set_ylabel('time\nabove '+'{:.2f}'.format(FRThr)+'\nHz (s)',fontsize=4.5)
  
  
for iax in range(13,15):
  axarr[iax].set_xticks([])
  axarr[iax].set_xlim([-0.5,1.5])
  axarr[iax].set_ylim([-2,117])

for iax in range(0,13):
  axarr[iax].set_xlim([-9.000,89.000])
  axarr[iax].set_xlabel('time (s)',fontsize=5)
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.04, pos.y1 - 0.0, chr(ord('A')+iax), fontsize=9)

for iax in [0,5,10,11,12]:  
  axarr[iax].set_ylabel('FR (spikes/s)',fontsize=5)

ml = mytools.mylegend(f,[0.2,0.975,0.6,0.0248],['k-','k-','k-'],['Control','ACC-like variant','PFC-like variant'],3,2,0.5,0.35,['#000000']+cols[::-1],linewidths=[0.5,0.5,0.5],myfontsize=5)
ml.spines['top'].set_visible(False)
ml.spines['right'].set_visible(False)
ml.spines['bottom'].set_visible(False)
ml.spines['left'].set_visible(False)


  
f.savefig("fig1.pdf")
