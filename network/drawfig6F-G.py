from pylab import *
import scipy.io
import mytools
from matplotlib.collections import PatchCollection
from os.path import exists
import scipy.stats

addition = ''

def drawdiscontinuity(ax,y,yoffset,x=0,xoffset=0.1,lw=2.0,lw2=1.0):
  thisline = ax.plot([x-xoffset,x+xoffset],[y-yoffset,y],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-xoffset,x+xoffset],[y,y+yoffset],'k-',linewidth=lw2)
  thisline[0].set_clip_on(False)
  thisline = ax.plot([x-1.5*xoffset,x+1.5*xoffset],[y-0.75*yoffset,y+0.75*yoffset],'k-',color='#FFFFFF',zorder=100,linewidth=lw)
  thisline[0].set_clip_on(False)

f,axarr = subplots(1,2)
axnew = []
for iax in range(0,2):
    axarr[iax].set_position([0.08+0.31*iax,0.78,0.25,0.18])
    axnew.append(f.add_axes([0.08+0.31*iax+0.21,0.78+0.065,0.03,0.105]))
for iax in range(0,2):
  axarr[iax].tick_params(axis='both', which='major', labelsize=4)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.3)
for iax in range(0,len(axnew)):
  axnew[iax].tick_params(axis='both', which='major', labelsize=4)
  for axis in ['top','bottom','left','right']:
    axnew[iax].spines[axis].set_linewidth(0.3)

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

ipop = 3 #Stimulated pyramidal cells

Nsum = 100
Ncellsperpop = [582,97,97,291,291]

areas = ['PFC','ACC']
iareas = [1,0] # Plot first ACC, then PFC
channels = ['CTRL', 'PConly']
channel_titles = ['CTRL','Combination']

isubjs_HC_PFC = [0, 1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20, 22, 23, 25, 27, 28, 31, 32, 33, 34, 35, 36, 37, 39, 41, 42, 43, 47, 51, 52, 53, 54, 63, 64, 66, 69, 73, 80, 81, 82, 83, 84, 87, 88, 89, 90, 92, 93, 94, 96, 97, 100, 105, 107, 111, 112, 118, 119, 120, 125, 127, 128, 135, 139, 140, 142, 143, 145, 146, 152, 154, 156, 165, 168, 176, 181, 182, 187, 192, 194, 201, 203, 204, 205, 207, 210, 211, 217, 220, 222, 223, 225, 226, 232, 235, 236, 240, 241, 242, 243, 244, 247, 248, 249, 251, 252, 253, 254, 257, 258, 262, 263, 265, 266, 267, 268, 269, 270, 271, 272, 283, 290, 291, 293, 294, 295, 296, 297, 298, 299, 300, 301, 305, 306, 308, 313, 314, 315, 317, 319, 320, 321, 323, 324, 325, 326, 329, 330, 332, 333, 335, 336, 339, 340, 341, 342, 343, 344, 345, 347, 351, 352, 354, 355, 356, 359, 361, 362, 368, 371, 373, 374, 375, 376, 377, 380, 381, 382, 383, 384, 386, 387, 388, 389, 390, 393, 396, 397, 398, 399, 400, 402, 403, 405, 406, 407, 408, 415, 416, 417, 420, 421, 422, 423, 424, 425, 426, 427, 428]
isubjs_SCZ_PFC = [4, 6, 15, 17, 21, 24, 26, 29, 30, 38, 40, 44, 45, 46, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 62, 65, 67, 68, 70, 71, 72, 74, 75, 76, 77, 78, 79, 85, 86, 91, 95, 98, 99, 101, 102, 103, 104, 106, 108, 109, 110, 113, 114, 115, 116, 117, 121, 122, 123, 124, 126, 129, 130, 131, 132, 133, 134, 136, 137, 138, 141, 144, 147, 148, 149, 150, 151, 153, 155, 157, 158, 159, 160, 161, 162, 163, 164, 166, 167, 169, 170, 171, 172, 173, 174, 175, 177, 178, 179, 180, 183, 184, 185, 186, 188, 189, 190, 191, 193, 195, 196, 197, 198, 199, 200, 202, 206, 208, 209, 212, 213, 214, 215, 216, 218, 219, 221, 224, 227, 228, 229, 230, 231, 233, 234, 237, 238, 239, 245, 246, 250, 255, 256, 259, 260, 261, 264, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 285, 286, 287, 288, 289, 292, 302, 303, 304, 307, 309, 310, 311, 312, 316, 318, 322, 327, 328, 334, 337, 338, 346, 348, 349, 350, 353, 357, 358, 360, 363, 364, 365, 366, 367, 369, 370, 372, 378, 379, 385, 392, 394, 395, 401, 404, 409, 410, 411, 412, 413, 414, 418, 419]
isubjs_HC_ACC = [1, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 73, 75, 76, 77, 78, 79, 80, 83, 84, 86, 90, 91, 92, 95, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 124, 126, 127, 129, 130, 131, 132, 133, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 152, 154, 159, 161, 163, 164, 167, 169, 170, 172, 174, 175, 176, 177, 182, 199, 208, 235, 240, 243, 249, 253, 256, 257, 260, 261, 262, 263, 265, 266, 267, 268, 271, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 329, 330, 331, 332, 343, 349, 350, 352, 353, 355, 360, 366, 369, 370, 371, 372, 374, 375, 376, 378, 381, 382, 383, 385, 388, 389, 392, 393, 395, 396, 397, 398, 401, 402, 409, 410, 414, 416, 418, 421, 422, 426, 427, 431, 432, 433, 434, 435, 437, 438, 442, 443, 444, 446, 447, 448, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480]
isubjs_SCZ_ACC = [0, 2, 3, 12, 17, 28, 31, 58, 69, 70, 71, 72, 74, 81, 82, 85, 87, 88, 89, 93, 94, 96, 97, 98, 99, 100, 116, 117, 118, 119, 120, 121, 122, 123, 125, 134, 151, 155, 156, 157, 158, 160, 162, 165, 166, 168, 171, 173, 178, 179, 180, 181, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 200, 201, 202, 203, 204, 205, 206, 207, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 236, 237, 238, 239, 241, 242, 244, 245, 246, 247, 248, 250, 251, 252, 254, 255, 258, 259, 264, 269, 270, 272, 273, 274, 275, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 344, 345, 346, 347, 348, 351, 354, 356, 357, 358, 359, 361, 362, 363, 364, 365, 367, 368, 373, 377, 379, 380, 384, 386, 387, 390, 391, 394, 399, 400, 403, 404, 405, 406, 407, 411, 412, 413, 415, 417, 419, 420, 423, 424, 425, 428, 429, 430, 436, 439, 440, 441, 445, 449, 450, 451, 452]


FRThrCoeff = 2.0

ipops = [3]
for iipop in [0]:
  ipop = ipops[iipop]
  
  groupmeans_all = []
  for iiarea in range(0,len(iareas)):
    iarea = iareas[iiarea]
    thisarea = areas[iarea]+'mild'
    
    #Subject-wise:
    isubjs_HC = isubjs_HC_PFC if iarea == 0 else isubjs_HC_ACC
    isubjs_SCZ = isubjs_SCZ_PFC if iarea == 0 else isubjs_SCZ_ACC
    groupmeans = []
    FRthisarea = []
    for igroup in [0,1]:
      isubjs = isubjs_HC if igroup == 0 else isubjs_SCZ
      FRvecs_all_means = []
      FRthisgroup = []
      for iisubj in range(0,len(isubjs)):
        isubj = isubjs[iisubj]
        bins_all = []
        FRvecs = []
        for myseed in [1]:
          filename = 'subjs/FRandMUA_sim_net_long_withSK3_'+thisarea+'subj'+str(isubj)+addition+'_altgain_PConly_CACNA1Ionly_seed'+str(myseed)+'.mat'
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
            FRvecs.append(FRs[:])
          elif iisubj == 0:
            print(filename+' does not exist')
        FRthisgroup.append(FRvecs[:])

        if iisubj%8 == 0:
          axarr[iiarea].plot(10*Nsum*array(range(0,len(FRs)))/1000-10, [mean([FRvecs[isamp][it] for isamp in range(0,len(FRvecs))]) for it in range(0,len(FRs))],color=dimcols[iarea] if igroup==1 else dimgray,lw=0.1)
        FRvecs_all_means.append([mean([FRvecs[isamp][it] for isamp in range(0,len(FRvecs))]) for it in range(0,len(FRs))])
      FRthisarea.append(FRthisgroup[:])
      groupmeans.append(mean(array(FRvecs_all_means),axis=0))
    
    FRbeyondThrthisarea = []
    FRbaseline = mean(groupmeans[0][3:10])
    groupmeans_all.append(groupmeans[:])
    for igroup in [0,1]:
      isubjs = isubjs_HC if igroup == 0 else isubjs_SCZ
      FRthisgroup = FRthisarea[igroup]
      FRbeyondThrthisgroup = []
      for iisubj in range(0,len(isubjs)):
        FRbeyondThrthissubj = []
        FRmean = [mean([FRthisgroup[iisubj][isamp][it] for isamp in range(0,len(FRthisgroup[iisubj]))]) for it in range(0,len(FRs))]

        FRThr = FRThrCoeff*FRbaseline
        FRbeyondThr = sum([1 for i in range(0,len(FRmean)) if FRmean[i] >= FRThr])
        FRbeyondThrthissubj.append(FRbeyondThr)
        FRbeyondThrthisgroup.append(FRbeyondThrthissubj[:])
      FRbeyondThrthisarea.append(FRbeyondThrthisgroup[:])

    for igroup in [0,1]:
      axarr[iiarea].plot(10*Nsum*array(range(0,len(FRs)))/1000-10,groupmeans[igroup],lw=0.75,color=cols[iarea] if igroup else '#000000',label=thisarea)
    
      mybar(axnew[iiarea],igroup,FRbeyondThrthisarea[igroup],facecolor=cols[iarea] if igroup else '#888888')
      if max(array(FRbeyondThrthisarea[igroup])) > 90:
        drawdiscontinuity(axnew[iiarea],89,5,x=igroup,xoffset=0.15,lw=0.6,lw2=0.3)
      print('Median time above 2*baseline, iipop='+str(iipop)+', iiarea='+str(iiarea)+', igroup='+str(igroup)+': '+str(median(FRbeyondThrthisarea[igroup]))+', mean : '+str(mean(FRbeyondThrthisarea[igroup])))

    pval = scipy.stats.ranksums(FRbeyondThrthisarea[0],FRbeyondThrthisarea[1])[1]
    if pval < 0.05:
      axnew[iiarea].plot([0,0,1,1],[60,65,65,60],'k-',lw=0.1)
      axnew[iiarea].plot(0.5,70,'k*',mew=0.9,lw=0.9,ms=3)
    else:
      axnew[iiarea].plot([0,0,1,1],[60,65,65,60],'k-',lw=0.1,color='#FFFFFF')
      axnew[iiarea].plot(0.5,70,'k*',mew=0.9,lw=0.9,ms=3,color='#FFFFFF')
    
    print("iipop="+str(iipop)+" iarea="+str(iarea)+" "+areas[iarea]+" p-val = "+str(pval))

    axnew[iiarea].set_ylabel('time\nabove '+'{:.2f}'.format(FRThr)+'\nHz (s)',fontsize=4.8)
    axnew[iiarea].set_xticks([])
    axnew[iiarea].set_xlim([-0.5,1.5])
    axnew[iiarea].set_ylim([-5,78])
    axnew[iiarea].set_yticks([0,30,60])
  for iax in range(0,2):
    axarr[iax].set_xlabel('time (s)',fontsize=6)
    axarr[iax].set_xlim([-9,89])
    pos = axarr[iax].get_position()
    f.text(pos.x0 - 0.05, pos.y1 - 0.0, chr(ord('F')+iax), fontsize=9)

for iax in [0,1]:
  axarr[iax].set_ylim([0,10])

axarr[0].set_title('Subject-wise CACNA1I-ACC models', fontsize=6)
axarr[1].set_title('Subject-wise CACNA1I-PFC models', fontsize=6)

axarr[0].set_ylabel('Stimulated PC\nFR (Hz)',fontsize=6)
  
f.savefig("fig6F-G.pdf")

