# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 15:27:59 2016

@author: porio
"""

#from __future__ import division

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import neo
import os
import Wavelets
import sys

def datosSurr(serie,N=10):
    L=len(serie)
    if L%2!=0:
        serie=np.append(serie,serie[-1])
    fftserie=np.fft.fft(serie)
#    ang=np.angle(fftserie)
    amp=np.abs(fftserie)
    
    angSurr=np.random.uniform(-np.pi,np.pi,size=(N,L))
    angSurr[:,L//2:]= - angSurr[:,L//2:0:-1]
    angSurr[:,L//2]=0
    
    fftdatosSurr=np.cos(angSurr)*amp + 1j*np.sin(angSurr)*amp
    
    datosS=np.real(np.fft.ifft(fftdatosSurr,axis=-1))
    
    return datosS

def datosSurr2(serie,N=10):
    L=serie.shape[-1]
    if L%2!=0:
        serie=np.append(serie,serie[:,-1][:,None],-1)
        L+=1
    fftserie=np.fft.fft(serie,axis=-1)
#    ang=np.angle(fftserie)
    amp=np.abs(fftserie)
    
    angSurr=np.random.uniform(-np.pi,np.pi,size=(N,)+serie.shape)
    angSurr[:,:,L//2:]= - angSurr[:,:,L//2:0:-1]
    angSurr[:,:,L//2]=0
    
    fftdatosSurr=np.cos(angSurr)*amp + 1j*np.sin(angSurr)*amp
    
    datosS=np.real(np.fft.ifft(fftdatosSurr,axis=-1))
    
    return datosS

def SurrXC(sig1,sig2,N=10):
    if len(sig1)!=len(sig2):
        return
    sig1ss=datosSurr(sig1,N)
    sig2ss=datosSurr(sig2,N)
    
    XCsurr=np.array([np.correlate(S1,S2,mode='same')/np.sqrt(np.dot(S1,S1)*np.dot(S2,S2))
            for S1 in sig1ss for S2 in sig2ss])
    
    Mxc = np.mean(XCsurr,0)
    SDxc = np.std(XCsurr,0)
    
    return (Mxc,SDxc)

def GaussSmooth(datawt,scales,fact=1):
    smoothWT=[]
#    fact=5
    for D,SC in zip(datawt,scales):
        wind=scipy.signal.gaussian(np.int(np.ceil(SC*fact))*2,SC*fact)
        smoothWT.append(scipy.signal.convolve(D,wind,mode='same'))
    return np.array(smoothWT)

def addcolumn(filename,newdata,header=None,formt=',%8g'):
    with open(filename,'r') as df:
        data=df.readlines()
    data=[d.strip() for d in data]
    
    newcol=[formt%di for di in newdata]
    if header is not None:
        newcol.insert(0,','+header)
    
    if len(newcol)!=len(data):
        return
    
    newdata=[da+db+'\n' for da,db in zip(data,newcol)]
    
    with open(filename,'w') as df:
        df.writelines(newdata)
        
def mAVG(data,window=10):
    center=window//2
    N=data.shape[-1] - center
    data2=np.rollaxis(data,-1)
    MAdata=np.array([np.mean(data2[i:i+window+1],axis=0) for i in range(0,N,window)])
    MAdata=np.rollaxis(MAdata,0,len(MAdata.shape))
    return MAdata

#%%

folder="TRUCHA_x/"
#filenames = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and '.wcp' in f]
# listfile="blacklist2b.txt"
listfile="blacklist2bTW.txt"

bltext=[]
with open(folder + listfile,'r') as fl:
    bllines=fl.read().splitlines()

for ln in bllines:
    bltext.append(ln.split(','))

print('the number of files in folder is =', len(bltext))
filenames=[bl[0] for bl in bltext]

names = ['EOG','Vv','Dp','OB']  #chequear el orden de los canales!!

pairs=((0,1),(1,2),(0,2))
labels=["%s-%s"%(names[pair[0]+1],names[pair[1]+1]) for pair in pairs]

#filelist=range(len(filenames))    #todos los archivos
#filelist=(1,5,8,15,22) #varios archivos
#filelist=(60,) #un solo archivo
#filelist=range(20,40) #varios
if len(sys.argv)>1:
    filelist=[int(a) for a in sys.argv[1:]]
else:
    filelist=(131,130,120,118,) #un solo archivo
#%%
print(filelist)
    
#%%


for i in filelist:    
    
    r = neo.io.WinWcpIO(filename='%s%s'%(folder,filenames[i]))
    bl = r.read_block()
    seg_size = np.shape(bl.segments)[0]
    print(filenames[i], 'The number of segments is =', seg_size)
    ts = np.transpose(np.array(bl.segments[0].analogsignals))
    
        
    dt=bl.segments[0].analogsignals[0].sampling_period
    data_series = np.array([b.analogsignals for b in bl.segments]).squeeze()
        
    t_start=bl.segments[0].analogsignals[0].t_start
    t_stop=bl.segments[0].analogsignals[0].t_stop
    dt=bl.segments[0].analogsignals[0].sampling_period
    
    time=np.arange(t_start-t_start,t_stop-t_start,dt)
    
    outdirname=folder + filenames[i].replace('.wcp','')
    if not os.path.isdir(outdirname):
        os.mkdir(outdirname)
    else:
        print("skipping"+outdirname)
        break
#        os.rename(outdirname,outdirname+'-all')
#        os.mkdir(outdirname)
    

    Fc=50 #[Hz] Frecuencia de corte
    decimate=5
    
    b,a=scipy.signal.butter(3,Fc*dt,btype='lowpass')
    filteredSignal=scipy.signal.filtfilt(b,a,data_series,axis=-1)
    sigMean=np.mean(filteredSignal,-1)
    filteredSignal=filteredSignal-sigMean[:,:,None]
    
    desiredFreqs=np.arange(0.5,40,0.5)  #Frecuencias a analizar con CWT - lineal
    #desiredFreqs=np.logspace(np.log10(0.5),np.log10(30),num=40,endpoint=True) #Log scale
    
    desiredPer=1 / (decimate * desiredFreqs * dt)
    dScales = np.array(desiredPer / Wavelets.Morlet.fourierwl)
    
    t1zoom=0
    t2zoom=25
    BLlimit=10
    SpecZoom=[12,15] #Time limit for spectrum calc
    
    zoom=(time>=t1zoom)*(time<t2zoom)
    zoom2=(time>=BLlimit)*(time<t2zoom)
    
    zoom=zoom.nonzero()[0][::decimate]  #Decimate to ~200 Hz
#    if len(zoom)%2!=0:
#        zoom=zoom[1:]
    zoom2=zoom2.nonzero()[0][::decimate]  #Decimate to ~200 Hz
#    if len(zoom2)%2!=0:
#        zoom2=zoom2[1:]

    
    WTsignal=[]
    PWsignal=[]
    swcount=0
    for sweep in filteredSignal: # Loop over epochs
        WTsweep=[]
        PWsweep=[]
        for chann in sweep[1:]: # Loop over channel
            cwt=Wavelets.Morlet(chann[zoom],scales=dScales) # chann[zoom]=(1,time); morlet=(time,freqs)
            WTsweep.append(cwt.getdata())
            PWsweep.append(cwt.getnormpower()) # Colecta morlets (time,freqs,chan)
        WTsignal.append(WTsweep)
        PWsignal.append(PWsweep) # "Mega array" (time,freqs,chan,sweep)
        print("listo CWT sweep %g"%swcount)
        swcount +=1
    
    WTsignal=np.array(WTsignal)
#    PWsignal=np.sqrt(np.array(PWsignal))
    PWsignal=np.array(PWsignal)
    
    BLindex=(time[::decimate]<BLlimit).nonzero()[0]
    zoom2index=((time[::decimate]>BLlimit) * (time[::decimate]<t2zoom)).nonzero()[0]
    SpecZoomIndex=((time[::decimate]>SpecZoom[0]) * (time[::decimate]<SpecZoom[1])).nonzero()[0]
    PWBLmean=np.mean(PWsignal[:,:,:,BLindex],-1)
    PWBLSD=np.std(PWsignal[:,:,:,BLindex],-1)
    
    PWsignal=(PWsignal - PWBLmean[:,:,:,None])/PWBLSD[:,:,:,None]
    #normWTsignal=WTsignal / dScales[None,None,:,None]
    #PWsignal=(normWTsignal * np.conjugate(normWTsignal)).real
    
    #plt.figure(0)
    #plt.clf()
    #for j in range(4):
    #    plt.subplot(4,1,j+1)
    #    plt.plot(time,data_series[:,j,:].T,color=colors[j])
    #    plt.ylabel(names[j])
    #        plt.yticks([])
    
    
    
    #%%
    
    
    XCmeas=["XCorrMax","XCorrMaxZ","XCorrMaxT","XCorrMin","XCorrMinZ","XCorrMinT"]
    XCheader=["%s %s"%(l1,l2)  for l2 in labels for l1 in XCmeas]
    
    XCorrResultFile=outdirname + "/XCorrTable.csv"
    if not os.path.isfile(XCorrResultFile):
        with open(XCorrResultFile,'w') as df:
            df.write("sweep,")
            for head in XCheader:
                df.write("%s,"%head)
            df.write("\n")
            
    CoherMaxFFile = [outdirname + "/CoherMaxF%s.csv"%lab for lab in labels]
    CoherFreqFile= [outdirname + "/CoherFreq%s.csv"%lab for lab in labels]
    PhDiffFile = [outdirname + "/PhaseDiff%s.csv"%lab for lab in labels]
    PLVFile= [outdirname + "/PLV%s.csv"%lab for lab in labels]
    SignalSpecFile=[outdirname + "/Spectrum%s.csv"%lab for lab in names[1:]]
    SignalMaxFTFile=[outdirname + "/MaxF%s.csv"%lab for lab in names[1:]]
    
    
    XcorrWindow=[12,15]  #Ventana para calcular cross-correlación
    maxLag = 0.5        #Tiempo en segundos del máximo delay a explorar

    
    CoherspecZoom=np.logical_and(time[zoom]>=XcorrWindow[0],time[zoom]<=XcorrWindow[1]).nonzero()[0]
    
    Mtime=mAVG(time[zoom2])
    MPWsignal=mAVG(PWsignal[:,:,:,zoom2index])
    
    for outf in CoherMaxFFile + SignalMaxFTFile:
        if not os.path.isfile(outf):
            with open(outf,'w') as df:
                df.write("time\n")
                for tt in Mtime:
                    df.write("%8g\n"%tt)
                    
    for outf in CoherFreqFile + SignalSpecFile + PLVFile:
        if not os.path.isfile(outf):
            with open(outf,'w') as df:
                df.write("freq\n")
                for tt in desiredFreqs:
                    df.write("%8g\n"%tt)

        
    cmap=plt.get_cmap('hsv')
    cmap.set_under('k',1)
    
    blacklist = [int(bsw) for bsw in bltext[i][1:]] #sweeps a ignorar
    whitelist = [x for x in range(len(filteredSignal)) if x not in blacklist]
                 
#    whitelist = (4,8,10,12,14,16,17,18)  #solo hacer algunos sweeps      

    TotalCoher=[[],[],[]]    
    TotalCoherSurr=[[],[],[]]     
    TotalXCorr=[[],[],[]]
    
    for sw in whitelist:
        with open(XCorrResultFile,'a') as df:
                df.write("%g,"%sw)

        print("File %g - %s sweep %g"%(i,filenames[i],sw))
        plt.figure(1,figsize=(15,10))
        plt.clf()
        
        for s in range(4):
            ax=plt.subplot(4,3,3*s+1)
            plt.plot(time,filteredSignal[sw,s,:],lw=0.5)
            ax.text(0.05,0.95,names[s],va='top',ha='left',
                    transform=ax.transAxes,size='large')
        
        for s in range(4):
            plt.subplot(4,3,3*s+2)
            plt.plot(time[zoom2],filteredSignal[sw,s,zoom2],lw=0.5)
            
        #PWmax=np.max(PWsignal[sw])    
        for s in range(1,4):
            plt.subplot(4,3,3*s+3)
            PWmax=np.max(PWsignal[sw,s-1])
            plt.imshow((PWsignal[sw,s-1][:,zoom2index]),interpolation='none',aspect='auto',origin='lower',
                       extent=(BLlimit,t2zoom,min(desiredFreqs),max(desiredFreqs)),
                       vmin=-PWmax,vmax=PWmax,cmap='RdBu')
            plt.colorbar(fraction=0.05,pad=0.02)
            
            PwSpec=np.mean(PWsignal[sw,s-1][:,SpecZoomIndex],-1)
            addcolumn(SignalSpecFile[s-1],PwSpec,header="sweep%g"%sw)
            
            MaxFTi=np.argmax(MPWsignal[sw,s-1],axis=0)
            addcolumn(SignalMaxFTFile[s-1],desiredFreqs[MaxFTi],header="MFreq%s"%sw)
            addcolumn(SignalMaxFTFile[s-1],MPWsignal[sw,s-1,MaxFTi,range(MPWsignal.shape[-1])],header="MFVal%s"%sw)
        
            
            
        plt.figtext(0.5, 0.98, "%s sweep %g"%(filenames[i],sw),ha='center',size='x-large')
    
        plt.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.05)
        
        plt.savefig(outdirname+'/'+'sweep%02g'%sw,dpi=300)
        
        plt.close()


        XcZoom=(time>=XcorrWindow[0])*(time<=XcorrWindow[1])
        XcL=len(time[XcZoom])
        maxLagI=np.int32(maxLag/dt)
        lagZoomI=range(XcL//2-maxLagI,XcL//2+maxLagI)
        lagZoomT=(np.array(lagZoomI)- lagZoomI[0])*np.float(dt)-maxLag
        
        print("Surrogate data for sweep %g"%sw)
        S1surr=datosSurr2(filteredSignal[sw,1:,zoom2].T,N=12)
        
        WTsurr=[]
        
        print("Surrogate Wavelet transforms for sweep %g"%sw)
        
        for sweep in S1surr:
            WTsweep=[]
            for chann in sweep:
                cwt=Wavelets.Morlet(chann,scales=dScales)
                WTsweep.append(cwt.getdata())
            WTsurr.append(WTsweep)
        
        WTsurr=np.array(WTsurr)
        
        plt.figure(2,figsize=(18,10))
        plt.clf()
        
        for s in range(3):
            S1=filteredSignal[sw,pairs[s][0]+1,XcZoom]
            S2=filteredSignal[sw,pairs[s][1]+1,XcZoom]
            
            XCorr=np.correlate(S1,S2,mode='same')/np.sqrt(np.dot(S1,S1)*np.dot(S2,S2))
        
            maxXC=lagZoomT[np.argmax(XCorr[lagZoomI])]
            minXC=lagZoomT[np.argmin(XCorr[lagZoomI])]
            
            meanXCsurr,SDXCsurr=SurrXC(S1,S2,N=12)
            meanXCsurr=meanXCsurr[lagZoomI]
            SDXCsurr=SDXCsurr[lagZoomI]    
            
            TotalXCorr[s].append(XCorr[lagZoomI])
            
            ax1=plt.subplot(3,4,4*s+1)
            ax1.plot(lagZoomT,XCorr[lagZoomI])
            if s==0:
                plt.title("cross-correlation")
            if s==2:
                plt.xlabel("Time delay")
            ax1.plot(lagZoomT,meanXCsurr,'b-',lw=2)
            ax1.plot(lagZoomT,meanXCsurr+2*SDXCsurr,'b--')
            ax1.plot(lagZoomT,meanXCsurr-2*SDXCsurr,'b--')
                
            plt.text(0.01,0.99,"%s\nmin: %.3g s\nmax:%.3g s"%(labels[s],minXC,maxXC),
                     va='top',transform=ax1.transAxes)
            
            W1=WTsignal[sw,pairs[s-1][0]][:,zoom2index]
            W2=WTsignal[sw,pairs[s-1][1]][:,zoom2index]
            Wx1x2 = W1 * np.conj(W2)
            Wx1x1 = np.abs(W1 * np.conj(W1))
            Wx2x2 = np.abs(W2 * np.conj(W2))
            Wphcoher = np.angle(Wx1x2)
            
            print("sweep %g, channels %g and %g Smoothing spectrograms "%(sw,pairs[s][0],pairs[s][1]))
#            smWx1x2=GaussSmooth(np.real(Wx1x2),dScales) + 1j*GaussSmooth(np.imag(Wx1x2),dScales)
            fact=5
            
            sScales=20*np.ones_like(dScales) #scale(s) for smoothing
            
            smWx1x2=GaussSmooth(Wx1x2,sScales,fact)
            smWx1x1=GaussSmooth(Wx1x1,sScales,fact)
            smWx2x2=GaussSmooth(Wx2x2,sScales,fact)

            Wcoher = (np.abs(smWx1x2))**2 / (smWx1x1*smWx2x2)
            
            TotalCoher[s].append(Wcoher)
#            Wcoher = (np.abs(Wx1x2))**2 / (Wx1x1*Wx2x2)
            
#            Wcoher = 2*(np.abs(smWx1x2))**2 / (smWx1x1**4 + smWx2x2**4)
#            Wcoher = 2*(np.abs(Wx1x2))**2 / (Wx1x1**4 + Wx2x2**4)

            #Calculation of significance level
            
            print("sweep %g, channels %g and %g Surrogate coherence "%(sw,pairs[s][0],pairs[s][1]))
            W1surr=WTsurr[:,pairs[s-1][0]]
            W2surr=WTsurr[:,pairs[s-1][1]]
            Wx1x2sur = W1surr[None,:,:,:] * np.conj(W2surr[:,None,:,:])
            Wx1x1sur = np.abs(W1surr * np.conj(W1surr))
            Wx2x2sur = np.abs(W2surr * np.conj(W2surr))
            
            print("sweep %g, channels %g and %g Smoothing surrogate cross-spectrogram "%(sw,pairs[s][0],pairs[s][1]))
            smWx1x2s=np.array([[GaussSmooth(Ws,sScales,fact) for Ws in Wx1x2surN] for Wx1x2surN in Wx1x2sur])
            print("sweep %g, channel %g Smoothing surrogate spectrogram "%(sw,pairs[s][0]))
            smWx1x1s=np.array([GaussSmooth(Ws,sScales,fact) for Ws in Wx1x1sur])
            print("sweep %g, channel %g Smoothing surrogate spectrogram "%(sw,pairs[s][1]))
            smWx2x2s=np.array([GaussSmooth(Ws,sScales,fact) for Ws in Wx2x2sur])
            
            Wcohersurr = (np.abs(smWx1x2s))**2 / (smWx1x1s[None,:,:,:]*smWx2x2s[:,None,:,:])

            TotalCoherSurr[s].append(Wcohersurr[0][:120//len(whitelist)])
#            MWcohersurr = np.mean(Wcohersurr,(0,1))
#            SDcohersurr = np.std(Wcohersurr,(0,1))
            
#           CoherThresh=np.max(WcoherMean)*0.05
#            CoherThresh=np.mean(Wcoher)*5
#            CoherThresh=MWcohersurr + 2*SDcohersurr
            CoherThresh=np.percentile(Wcohersurr,95,axis=(0,1))
            
            CoherMask= (Wcoher>0.5)*(Wcoher>CoherThresh)
            
            Ang1 = np.angle(WTsignal[sw,pairs[s][0]][:,zoom2index])
            Ang2 = np.angle(WTsignal[sw,pairs[s][1]][:,zoom2index])
            AngDiff = Ang1-Ang2
            AngDiff = ( AngDiff + np.pi) % (2 * np.pi ) - np.pi
    #        AngDiff[AngDiff>np.pi] = AngDiff[AngDiff>np.pi]%np.pi
    #        AngDiff[AngDiff<np.pi] = AngDiff[AngDiff<np.pi]%np.pi
            AngDiff[np.logical_not(CoherMask)]= np.nan#-10
                        
            Wphcoher[np.logical_not(CoherMask)]=np.nan #-10
            
            plt.subplot(3,4,4*s+2)
            plt.imshow(Wcoher,interpolation='none',aspect='auto',origin='lower',
                       extent=(BLlimit,t2zoom,min(desiredFreqs),max(desiredFreqs)),cmap='jet')#,vmin=-np.pi,vmax=np.pi)
            plt.colorbar(fraction=0.05,pad=0.02)
            plt.contour(1*CoherMask,levels=[0.5],colors='k',antialiased=True,
                        extent=(BLlimit,t2zoom,min(desiredFreqs),max(desiredFreqs)))
            
            if s==0:
                plt.title("Coherence")
            if s==2:
                plt.xlabel("Time (s)")            
            plt.ylabel("Frequency (Hz)")
                
            plt.subplot(3,4,4*s+3)
            plt.imshow(AngDiff,interpolation='none',aspect='auto',origin='lower',
                       extent=(BLlimit,t2zoom,min(desiredFreqs),max(desiredFreqs)),
                       cmap=cmap,vmin=-np.pi,vmax=np.pi)
            if s==0:
                plt.title("Phase difference")
            if s==2:
                plt.xlabel("Time (s)")            
            
            plt.subplot(3,4,4*s+4)
            plt.imshow(Wphcoher,interpolation='none',aspect='auto',origin='lower',
                       extent=(BLlimit,t2zoom,min(desiredFreqs),max(desiredFreqs)),
                       vmin=-np.pi,vmax=np.pi,cmap=cmap)
            plt.text(11,35,labels[s])
            plt.colorbar(fraction=0.05,pad=0.02)
            if s==0:
                plt.title("Coherence phase")
            if s==2:
                plt.xlabel("Time (s)")
                
            
            XCorrP=np.diff(XCorr[lagZoomI])
            XCorr_idn=(np.diff(1*(XCorrP>0))==1).nonzero()[0]
            XCorr_idp=(np.diff(1*(XCorrP>0))==-1).nonzero()[0]
            if len(XCorr_idp)>1:
                closestMax_i=np.argmin(np.abs((XCorr_idp+1)-len(lagZoomI)/2.))
                closestMaxT=lagZoomT[XCorr_idp+1][closestMax_i]
                closestMaxV=XCorr[lagZoomI][XCorr_idp+1][closestMax_i]
                closestMaxConf=(closestMaxV - meanXCsurr[XCorr_idp+1][closestMax_i])/SDXCsurr[XCorr_idp+1][closestMax_i]
            else:
                closestMax_i=np.nan
                closestMaxT=np.nan
                closestMaxV=np.nan
                closestMaxConf=np.nan
                
            if len (XCorr_idn)>1:
                closestMin_i=np.argmin(np.abs((XCorr_idn+1)-len(lagZoomI)/2.))
                closestMinT=lagZoomT[XCorr_idn+1][closestMin_i]
                closestMinV=XCorr[lagZoomI][XCorr_idn+1][closestMin_i]
                closestMinConf=(closestMinV - meanXCsurr[XCorr_idn+1][closestMin_i])/SDXCsurr[XCorr_idn+1][closestMin_i]
            else:
                closestMin_i=np.nan
                closestMinT=np.nan
                closestMinV=np.nan
                closestMinConf=np.nan
            
            with open(XCorrResultFile,'a') as df:
                df.write(("%g,"*6)%(closestMaxV,closestMaxConf,closestMaxT,
                                    closestMinV,closestMinConf,closestMinT))
                
            mWcoher=mAVG(Wcoher*CoherMask)
            
            mWphcoher=[np.arctan2(np.sum(np.sin(Wphcoher[:,ii:ii+11]),-1),np.sum(np.cos(Wphcoher[:,ii:ii+11]),-1))
                        for ii in range(0,Wphcoher.shape[-1]-5,10)]
            mWphcoher=np.array(mWphcoher).T
        
#            center=window//2
#            N=data.shape[-1] - center
#            data2=np.rollaxis(data,-1)
#            MAdata=np.array([np.mean(data2[i:i+window+1],axis=0) for i in xrange(0,N,window)])
#            MAdata=np.rollaxis(MAdata,0,len(MAdata.shape))
                    
            mCoheri=np.argmax(mWcoher,0)
            addcolumn(CoherMaxFFile[s],desiredFreqs[mCoheri],header="maxF %g"%sw)
            addcolumn(CoherMaxFFile[s],mWcoher[mCoheri,range(mWphcoher.shape[-1])],header="maxFval %g"%sw)
#            addcolumn(CoherMaxFTAngFile[s],mWphcoher[mCoheri,range(mWphcoher.shape[-1])],header="PhCoher %g"%sw)
            
            CoherSpec=np.mean((Wcoher*CoherMask)[:,CoherspecZoom],-1)
            maxFreqi=np.argmax(CoherSpec)
            
            addcolumn(CoherFreqFile[s],CoherSpec,header="sweep %s"%sw)
#            addcolumn(CoherMaxFAngFile[s],mWphcoher[maxFreqi,:],header="sweep %g"%sw)
                
            PLVfreq=np.abs(np.nanmean(np.exp(1j*AngDiff),axis=-1))
            addcolumn(PLVFile[s],PLVfreq,header="sweep %g"%sw)
            
        plt.figtext(0.5, 0.98, "%s sweep %g plot 2"%(filenames[i],sw),ha='center',size='x-large')
        
    
    
        plt.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.05)
        
        
        with open(XCorrResultFile,'a') as df:
                df.write("\n")
        
        
        plt.savefig(outdirname+'/'+'sweep%02g plot2'%sw,dpi=300)

        plt.close()     
#%%        
    PWsignalAVG=np.mean(PWsignal[whitelist,:,:,:],axis=0)
    MPWsignalAVG=mAVG(PWsignalAVG[:,:,zoom2index])
    
    TotalCoher=np.array([np.mean(xx,axis=0) for xx in TotalCoher])
    
    TotalCoherSurr=np.array(TotalCoherSurr)
    dims=TotalCoherSurr.shape
    TotalCoherSurr=[np.reshape(TCS,(-1,dims[3],dims[4])) for TCS in TotalCoherSurr]
    CoherThresh=np.array([np.percentile(TotalCoherSurr[s],95,axis=(0)) for s in range(len(TotalCoherSurr))])
    
    np.save(outdirname+'/'+'Specgrams',PWsignalAVG[:,:,zoom2index])
    np.save(outdirname+'/'+'Cohergrams',TotalCoher)
    np.save(outdirname+'/'+'Coher95Thres',CoherThresh)
    
    TotalXCorr=np.array([np.mean(xx,axis=0) for xx in TotalXCorr])
    with open(XCorrResultFile,'a') as df:
        df.write("AVG,")
#%%    
    plt.figure(1,figsize=(15,10))
    plt.clf()
        
    for s in range(3):
        ax=plt.subplot(3,3,3*s+1)
        PWmax=np.max(PWsignalAVG[s])
        plt.imshow((PWsignalAVG[s][:,zoom2index]),interpolation='none',aspect='auto',origin='lower',
                   extent=(BLlimit,t2zoom,min(desiredFreqs),max(desiredFreqs)),
                   vmin=-PWmax,vmax=PWmax,cmap='RdBu')
        plt.colorbar(fraction=0.05,pad=0.02)
        ax.text(0.05,0.95,names[s+1],va='top',ha='left',
                    transform=ax.transAxes,size='large')
        
        PwAVGSpec=np.mean(PWsignalAVG[s][:,SpecZoomIndex],-1)
        addcolumn(SignalSpecFile[s],PwAVGSpec,header="AVG")
        
                
        MaxFTi=np.argmax(MPWsignalAVG[s],axis=0)
        addcolumn(SignalMaxFTFile[s],desiredFreqs[MaxFTi],header="MFreqAVG")
        addcolumn(SignalMaxFTFile[s],MPWsignalAVG[s,MaxFTi,range(MPWsignalAVG.shape[-1])],header="MFValAVG")
    
#        CoherThresh=np.percentile(TotalCoherSurr[s],95,axis=(0))
            
        CoherMask= TotalCoher[s]>CoherThresh[s]
        
        ax=plt.subplot(3,3,3*s+2)
        plt.imshow((TotalCoher[s]),interpolation='none',aspect='auto',origin='lower',
                   extent=(BLlimit,t2zoom,min(desiredFreqs),max(desiredFreqs)),
                   vmin=0,vmax=1,cmap='jet')
        plt.colorbar(fraction=0.05,pad=0.02)
        plt.text(11,35,labels[s],color='white',size='large')
        plt.contour(1*CoherMask,levels=[0.5],colors='k',antialiased=True,
                        extent=(BLlimit,t2zoom,min(desiredFreqs),max(desiredFreqs)))
        
        CoherSpec=np.mean(TotalCoher[s,:,CoherspecZoom],-1)
        addcolumn(CoherFreqFile[s],CoherSpec,header="AVG")

        ax1=plt.subplot(3,3,3*s+3)
        plt.plot(lagZoomT,TotalXCorr[s])
        plt.ylim((-1,1))
        if s==0:
            plt.title("cross-correlation")
        if s==2:
            plt.xlabel("Time delay")
            
        maxXC=lagZoomT[np.argmax(TotalXCorr[s])]
        minXC=lagZoomT[np.argmin(TotalXCorr[s])]
        
        plt.text(0.01,0.99,"%s\nmin: %.3g s\nmax:%.3g s"%(labels[s],minXC,maxXC),
                 va='top',transform=ax1.transAxes)
        
        
        XCorrP=np.diff(TotalXCorr[s])
        XCorr_idn=(np.diff(1*(XCorrP>0))==1).nonzero()[0]
        XCorr_idp=(np.diff(1*(XCorrP>0))==-1).nonzero()[0]
        if len(XCorr_idp)>1:
            closestMax_i=np.argmin(np.abs((XCorr_idp+1)-len(lagZoomI)/2.))
            closestMaxT=lagZoomT[XCorr_idp+1][closestMax_i]
            closestMaxV=XCorr[lagZoomI][XCorr_idp+1][closestMax_i]
            closestMaxConf=(closestMaxV - meanXCsurr[XCorr_idp+1][closestMax_i])/SDXCsurr[XCorr_idp+1][closestMax_i]
        else:
            closestMax_i=np.nan
            closestMaxT=np.nan
            closestMaxV=np.nan
            closestMaxConf=np.nan
            
        if len (XCorr_idn)>1:
            closestMin_i=np.argmin(np.abs((XCorr_idn+1)-len(lagZoomI)/2.))
            closestMinT=lagZoomT[XCorr_idn+1][closestMin_i]
            closestMinV=XCorr[lagZoomI][XCorr_idn+1][closestMin_i]
            closestMinConf=(closestMinV - meanXCsurr[XCorr_idn+1][closestMin_i])/SDXCsurr[XCorr_idn+1][closestMin_i]
        else:
            closestMin_i=np.nan
            closestMinT=np.nan
            closestMinV=np.nan
            closestMinConf=np.nan
        
        with open(XCorrResultFile,'a') as df:
            df.write(("%g,"*6)%(closestMaxV,closestMaxConf,closestMaxT,
                                closestMinV,closestMinConf,closestMinT))
            
    with open(XCorrResultFile,'a') as df:
        df.write("\n")
        
    plt.tight_layout()
    
    plt.savefig(outdirname+'/'+'Averages plot',dpi=300)
        
    #    ax.set_yscale('log')
    #    ax.set_yticks((1,3,5,10,30))
    #    ax.set_yticklabels(('1','3','5','10','30'))
