#!/usr/bin/python3 
# math_util namespace 

import os 
import sys 
import math
import numpy as np
import scipy 
 
# from utildf import logger 
import algorithm

#______________________________________________________________________________
def factorize(n):
   # factorize the number and return a list of values 
   # this is Pollard's rho integer factorization algorithm
   # uses the helper function get_factor (below) 
   # reference: https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm 
   factors = []
   while n>1:
       next = get_factor(n)
       factors.append(next)
       n //= next
   return factors
#___________________________________________________________________________
def get_factor(n):
    x_fixed = 2
    cycle_size = 2
    x = 2
    factor = 1
    while factor == 1:
        for count in range(cycle_size):
            if factor > 1: break
            x = (x * x + 1) % n
            factor = math.gcd(x - x_fixed, n)
        cycle_size *= 2
        x_fixed = x
    return factor

#_______________________________________________________________________________
def getSpectrogram(y,Fs,mode='psd',scale='dB',noverlap=-1,nfft=-1):
   # compute a spectrogram (time,freq,pwr) using matplotlib 
   # input 
   # - y        = time series of data 
   # - Fs       = sampling frequency [Hz] 
   # - mode     = type of power computation (e.g., psd) 
   # - scale    = dB or linear 
   # - noverlap = number of segments overlapping in calculation 
   # - nfft     = number of FFT bins 
   # output 
   # - time     = array of times [sec]    
   # - freq     = array of frequencies [Hz]   
   # - pxx      = array of power in dB [dB/sqrt(Hz)] or linear [V^2/Hz] 

   if nfft==-1: 
      # default to Fs 
      NFFT = int(Fs) 
   else: 
      NFFT = nfft 

   if noverlap==-1:
      nol = int(nfft/2) 
   else: 
      nol = noverlap

   # use scipy since it doesn't try to plot the data 
   freq,time,pxx = scipy.signal.spectrogram(y,fs=Fs,nfft=NFFT,nperseg=int(Fs),noverlap=nol,mode=mode)

   if scale=='dB':
      pxx = 10.*np.log10(pxx)

   return time,freq,pxx
#_______________________________________________________________________________
def getFFT(y,Fs):
   # get FFT of signal 
   # input: 
   # - y  = sampled signal in time domain [V] 
   # - Fs = sampling frequency [Hz] 
   # output: 
   # freq = frequency bins [Hz] 
   # pwr  = power density [V^2/Hz]  

   # fft 
   ffty  = np.fft.fft(y)
   # power spetrum after real2complex transform
   NY = len(y)
   sf = 2.0/(NY*NY)

   # get squared amplitude of ffty  
   pwr[i] = np.zeros( len(ffty) )
   for i in range(0,NY):
      arg    = ffty.real**2 + ffty.imag**2
      pwr[i] = sf*arg

   # get frequencies
   Ts   = 1./Fs  # sampling period 
   freq = np.fft.fftfreq(NY,Ts)

   return freq,pwr
#_______________________________________________________________________________
def getIFFT(y:np.array):
   # get the inverse FFT of signal y(f)
   return np.fft.ifft(y)  
#_______________________________________________________________________________
def getPSD_2(y:np.array,Fs:int,window='hamming',nperseg:int=-1,onesided:bool=False,
        scaling:str='density',dB:bool=False,debug:bool=False):
   # compute the power spectral density 
   # input 
   # - y        = numpy array of data samples 
   # - Fs       = sampling frequency [Hz]
   # - window   = windowing method (e.g., boxcar, hamming)   
   # - nperseg  = number of samples per segment in welch calculation 
   # - scaling  = density (output in V^2/Hz) or spectrum (output in V^2)  
   # - onesided = True  => return one-sided spectrum for real data [0,...,Fs/2]  
   #              False => return two-sided spectrum for complex data [0,...,Fs/2,-Fs/2,...,-1]  
   # output: 
   # - FREQ = array of frequencies [Hz] 
   # - PWR  = array of power densities based on scaling input  

   if nperseg==-1:
      # default to sampling frequency => 1 Hz bins 
      nps = Fs
   else: 
      nps = nperseg  

   freq,pwr = scipy.signal.welch(y,fs=Fs,window=window,nperseg=nps,
                                 scaling=scaling,axis=-1,average='mean',
                                 return_onesided=onesided)

   Fs2  = Fs/2.
   N    = len(freq)
   FREQ = np.zeros(N)
   P    = np.zeros(N)
   ii=0

   if onesided:
      # one-sided, save the data as-is 
      FREQ = freq
      P    = pwr
   else:
      # not one-sided! re-arrange to be freq = [-Fs/2,...,0,...,Fs/2]  
      for i in range(0,N):
         diff = freq[i]-Fs2
         if abs(diff)<2:
            ii = i
            # if debug:
            #    msg = "Found Fs/2! freq[{0}] = {1}".format(i,freq[i])
            #    logger.debug(funcName,msg)
            #    msg = "-Fs/2: freq[{0}] = {1}".format(i+1,freq[i+1])
            #    logger.debug(funcName,msg)
      # extract negative frequencies and corresponding psd
      f2 = freq[ii+1:]
      p2 = pwr[ii+1:]
      # extract positive frequencies and corresponding psd
      f1 = freq[:ii+1]
      p1 = pwr[:ii+1]
      if debug:
         print("First half of freq:")
         print(f1)
         print("Second half of freq:")
         print(f2)
      # build single arrays for each 
      FREQ = np.append(f2,f1)
      P    = np.append(p2,p1)

   if dB==True:
      PWR = 10.*np.log10(P)
   else:
      PWR = P

   return FREQ,PWR
#_______________________________________________________________________________
def getPSD(y,Fs,window='hamming',nperseg=-1,onesided=False,scale='linear',debug=False):
   # compute the power spectral density 
   # input 
   # - y        = numpy array of data samples 
   # - Fs       = sampling frequency [Hz]
   # - window   = windowing method (e.g., boxcar, hamming)   
   # - nperseg  = number of samples per segment in welch calculation 
   # - scale    = linear or dB 
   # - onesided = True  => return one-sided spectrum for real data [0,...,Fs/2]  
   #              False => return two-sided spectrum for complex data [0,...,Fs/2,-Fs/2,...,-1]  
   # output: 
   # - FREQ = array of frequencies [Hz] 
   # - PWR  = array of power densities [V^2/Hz] or [dB/Hz] based on scale input  

   funcName = "util.getPSD"
   
   if nperseg==-1:
      # default to sampling frequency => 1 Hz bins 
      nps = Fs
   else: 
      nps = nperseg  

   freq,pwr = scipy.signal.welch(y,fs=Fs,window=window,nperseg=nps,
                                 scaling='density',axis=-1,average='mean',
                                 return_onesided=onesided)

   Fs2  = Fs/2.
   N    = len(freq)
   FREQ = np.zeros(N)
   P    = np.zeros(N)
   ii=0

   if onesided:
      # one-sided, save the data as-is 
      FREQ = freq
      P    = pwr
   else:
      # not one-sided! re-arrange to be freq = [-Fs/2,...,0,...,Fs/2]  
      for i in range(0,N):
         diff = freq[i]-Fs2
         if abs(diff)<2:
            ii = i
            # if debug:
            #    msg = "Found Fs/2! freq[{0}] = {1}".format(i,freq[i])
            #    logger.debug(funcName,msg)
            #    msg = "-Fs/2: freq[{0}] = {1}".format(i+1,freq[i+1])
            #    logger.debug(funcName,msg)
      # extract negative frequencies and corresponding psd
      f2 = freq[ii+1:]
      p2 = pwr[ii+1:]
      # extract positive frequencies and corresponding psd
      f1 = freq[:ii+1]
      p1 = pwr[:ii+1]
      if debug:
         print("First half of freq:")
         print(f1)
         print("Second half of freq:")
         print(f2)
      # build single arrays for each 
      FREQ = np.append(f2,f1)
      P    = np.append(p2,p1)

   if scale=='dB':
      PWR = 10.*np.log10(P)
   elif scale=='linear':
      PWR = P

   return FREQ,PWR
#_______________________________________________________________________________
def linearInterp_timeInt(v0,V,F,fs,dt,axis='x',debug=False,thr=1E-3):
   # linear interpolation to estimate F(T,x0,Y) and then integrate/average over 
   # the time axis to obtain F(Y).  Or can obtain F(X).  
   # input 
   # - v0 = desired coordinate 
   # - V  = array X or Y; X = horizontal position array, Y = vertical position array 
   # - F  = 3D array F(T,X,Y); T = time index
   # - fs = sample frequency [Hz]  
   # - dt = full duration of time of the signal
   # - axis = axis of interpolation (x or y) 
   # output 
   # - G  = array interpolated to F(T,v0,Y), then integrated over time 
   #        to produce G(Y) (or similarly G(X))  

   funcName = "math_util.linearInterp_timeInt"

   # first do the linear interpolation 
   FF = linearInterp(v0=v0,V=V,F=F,axis=axis,debug=debug) # FF is function of time and one axis 

   # now integrate over time
   NP = len(FF[0,:])     # number of points in position axis 
   NT = len(FF[:,0])     # length of time axis for first position index. time axis size is constant 

   ts = 1./fs # size of time step
   SUM=0 

   # if debug:
   #    msg = "time step = {0:.3E} sec, dt = {1:.3E} sec".format(ts,dt) 
   #    logger.debug(funcName,msg) 
     
   # output array  
   G = np.zeros(NP) 

   # integrate over time 
   for i in range(0,NP): 
      for j in range(0,NT):
         SUM = SUM + ts*FF[j,i]
      # save result and clear   
      G[i] = (1./dt)*SUM 
      # clear the sum
      SUM = 0 

   return G
#_______________________________________________________________________________
def linearInterp_timeRMS(v0,V,F,axis='x',debug=False,thr=1E-3):
   # linear interpolation to estimate F(T,x0,Y) and then take RMS of  
   # the time axis to obtain F(Y).  Or can obtain F(X).  
   # input 
   # - v0 = desired coordinate 
   # - V  = array X or Y; X = horizontal position array, Y = vertical position array 
   # - F  = 3D array F(T,X,Y); T = time index
   # - axis = axis of interpolation (x or y) 
   # output 
   # - G  = array interpolated to F(T,v0,Y), then integrated over time 
   #        to produce G(Y) (or similarly G(X))  

   funcName = "math_util.linearInterp_timeRMS"

   # first do the linear interpolation 
   FF = linearInterp(v0=v0,V=V,F=F,axis=axis,debug=debug) # FF is function of time and one axis 

   # now integrate over time
   NP = len(FF[0,:])     # number of points in position axis 
   NT = len(FF[:,0])     # length of time axis for first position index. time axis size is constant 
     
   # output array  
   G = np.zeros(NP) 

   # compute RMS over time axis  
   rms = 0  
   for i in range(0,NP): 
      g    = FF[:,i] 
      G[i] = np.sqrt( g.dot(g)/g.size )
      # clean up 
      del g 

   return G
#_______________________________________________________________________________
def linearInterp(v0,V,F,axis='x',debug=False,thr=1E-3): 
   # linear interpolation to estimate F(T,x0,Y) 
   # input 
   # - v0 = desired coordinate 
   # - V  = array X or Y; X = horizontal position array, Y = vertical position array 
   # - F  = 3D array F(T,X,Y); T = time index 
   # - axis = axis of interpolation (x or y) 
   # output 
   # - the 2D array F(T,X) or F(T,Y) evaluated at v0  

   funcName = "math_util.linearInterp" 

   ilo,ihi = algorithm.binarySearch(V,v0) 

   # get bounding values 
   vl = V[ilo] 
   vh = V[ihi] 
   
   # get arrays evaluated at the identified boundaries  
   if axis=='x':
      FL = F[:,ilo,:]
      FH = F[:,ihi,:]
   elif axis=='y': 
      FL = F[:,:,ilo]
      FH = F[:,:,ihi]
   else: 
      msg = "Invalid interpolation axis = {0}".format(axis) 
      print(msg)
      # logger.error(funcName,msg)
      sys.exit(1)  

   # put everything together 
   b  = (v0-vl)/(vh-vl)
   F0 = FL + b*(FH-FL)  

   # if debug: 
   #     msg = "V[{0:d}] = {1}, v0 = {2}, V[{3:d}] = {4}".format(ilo,vl,v0,ihi,vh) 
   #     logger.debug(funcName,msg) 
   #     msg = "FL = {0}".format(FL) 
   #     logger.debug(funcName,msg)  
   #     msg = "FH = {0}".format(FH) 
   #     logger.debug(funcName,msg)  
   #     msg = "F0 = {0}".format(F0) 
   #     logger.debug(funcName,msg)  

   return F0
#_______________________________________________________________________________
def linearInterp(x0:float,X:np.array,Y:np.array):
   # linear interpolation to find y(x0) given the arrays X and Y  
   # find bounding indices for x0 in X  
   lo,hi = algorithm.binarySearch(X,x0)
   # bounding values  
   x_lo  = X[lo]
   x_hi  = X[hi] 
   y_lo  = Y[lo] 
   y_hi  = Y[hi] 
   # construct result
   b     = (x0-x_lo)/(x_hi-x_lo) 
   y     = y_lo + b*(y_hi-y_lo)
   return y
#_______________________________________________________________________________
def bilinearInterp(x0,y0,X,Y,F,debug=False,thr=1E-3):
   # bilinear interpolation to estimate F(x0,y0)  
   # input 
   # - desired coordinate x0, y0 
   # - arrays X, Y, F representing the grid F(T,X,Y) with T = time array 
   # optional input
   # - debug: print debug info (default = false)  
   # - thr: threshold for determining corners of F(T,x0,y0) given input x0,y0  
   # output 
   # - the 1D array F(T) evaluated at (x0,y0) 

   funcName = "math_util.bilinearInterp"
 
   # find bounding values for x0 and y0  
   ixlo,ixhi = algorithm.binarySearch(X,x0)
   iylo,iyhi = algorithm.binarySearch(Y,y0)
   
   # get bounding (x,y) values
   xl = X[ixlo];
   xh = X[ixhi];
   yl = Y[iylo];
   yh = Y[iyhi];
 
   # weight factors for interpolation 
   xwl = (x0-xl)/(xh-xl);
   xwh = (xh-x0)/(xh-xl);
   ywl = (y0-yl)/(yh-yl);
   ywh = (yh-y0)/(yh-yl);

   # construct arrays F(T,xlo,ylo), F(T,xlo,yhi), F(T,xhi,ylo), F(T,xhi,yhi)
   Fll = F[:,ixlo,iylo]
   Flh = F[:,ixlo,iyhi]
   Fhl = F[:,ixhi,iylo] 
   Fhh = F[:,ixhi,iyhi] 

   # put everything together 
   F00 = ywh*(xwh*Fll + xwl*Fhl) + ywl*(xwh*Flh + xwl*Fhh) 
  
   # if debug: 
   #    msg = "x[{0:d}] = {1}, x0 = {2}, x[{3:d}] = {4}".format(ixlo,xl,x0,ixhi,xh) 
   #    logger.debug(funcName,msg) 
   #    msg = "y[{0:d}] = {1}, y0 = {2}, y[{3:d}] = {4}".format(iylo,yl,y0,iyhi,yh) 
   #    logger.debug(funcName,msg)
   #    msg = "Fll = {0}".format(Fll) 
   #    logger.debug(funcName,msg)  
   #    msg = "Flh = {0}".format(Flh) 
   #    logger.debug(funcName,msg)  
   #    msg = "Fhl = {0}".format(Fhl) 
   #    logger.debug(funcName,msg)  
   #    msg = "Fhh = {0}".format(Fhh) 
   #    logger.debug(funcName,msg)  
   #    msg = "F00 = {0}".format(F00) 
   #    logger.debug(funcName,msg)  
 
   return F00 
