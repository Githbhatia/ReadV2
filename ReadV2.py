import numpy as np
import matplotlib.pyplot as plt
from itertools import islice
from tkinter import *
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter import messagebox
from RS_function import RS_function
import urllib.request as ur
import json as js
import webbrowser
import os
import math

def openMap():
    webbrowser.open('http://www.google.com/maps/place/'+ str(latitude) +','+str(longitude)+'/@'+ str(latitude) +','+str(longitude)+',12z', new=2)

def openStation():
    station = os.path.basename(filenames[0])
    station = station[:station.find(".v2")].strip()
    if station[:2] == "ce" or station[:2] == "CE" :
        webbrowser.open('https://www.strongmotioncenter.org/cgi-bin/CESMD/stationhtml.pl?stationID='+station+'&network=CGS')
    else:
        messagebox.showinfo('ReadV2', 'Non CSMIP freefield station, no info available')

def absmaxND(a):
    amax = np.max(a)
    amin = np.min(a)
    return np.where(-amin > amax, -amin, amax)

def argmedian(x):
  return np.argpartition(x, len(x) // 2)[len(x) // 2]

def chunkstring(string, length):
    return (float(string[0+i:length+i]) for i in range(0, len(string), length))

def readchunk(f, numofLines):
    x=[]
    for line in islice(f, 0,  numofLines):
        x = x + list(chunkstring(line[0:len(line)-1], 10))
    #print(x)
    return x

def lines(points):
    if points % 8 == 0:
        nLines = int(points/8) 
    else:
        nLines = int(points/8)+1
    return nLines

def scaleValue(units):
    if units =="cm/sec2":
        return 1/980.665
    else:
        return 1.0

def maxaccel(x, t):
    ymax = max(x)
    xpos = x.index(ymax); xmax = t[xpos]
    return [xmax, ymax]

def maxSA(x, t):
    ymax = np.amax(x)
    xpos = np.argmax(x); xmax = t[xpos]
    return [xmax, ymax]

def minaccel(x, t):
    ymin = min(x)
    xpos = x.index(ymin); xmin = t[xpos]
    return [xmin, ymin]

def accelim(x,y,z):
    xmax = max([abs(i) for i in x])
    ymax = max([abs(i) for i in y])
    zmax = max([abs(i) for i in z])
    return max([xmax,ymax,zmax])



def on_click():
    #%% Parameters of the response spectra

    if EOF == 1:
        noSubplotsRows = 1 + canvas.plotVel.get() + canvas.plotDisp.get() + canvas.plotFFT.get()
        noSubplotsCols = 1
        subplotCounter = 1
        plt.figure(1, figsize=(14,10))
        plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
        plt.grid()
        plt.title(nameCh1)
        plt.xlabel('Time (secs)')
        plt.ylabel('Accel(g)')
        plt.plot(T1,scaledAccel1, label="Channel1", color= 'Red', linewidth=1.0)
        amax=maxaccel(scaledAccel1, T1); plt.annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(scaledAccel1, T1); plt.annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
        plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])

        if str(canvas.plotVel.get()) =="1":
            subplotCounter+=1
            plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Velocity '+ unitsVel1)
            plt.plot(T1,vel1, label="Channel1", color= 'Blue', linewidth=1.0)
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])

        if str(canvas.plotDisp.get()) =="1":
            subplotCounter+=1
            plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel("Disp "+ unitsDispl1)
            plt.plot(T1,displ1, label="Channel1", color= 'Green', linewidth=1.0)
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])

        if str(canvas.plotFFT.get()) =="1":
            fourierTransform = np.fft.fft(scaledAccel1)/len(scaledAccel1)   # Normalize amplitude
            fourierTransform = fourierTransform[range(int(len(scaledAccel1)/2))] # Exclude sampling frequency
            tpCount     = len(scaledAccel1)
            values      = np.arange(int(tpCount/2))
            timePeriod  = tpCount*dtAccel1
            frequencies = values/timePeriod

            subplotCounter+=1
            plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            #plt.title('FFT of Acceleration')
            plt.grid()
            plt.plot(frequencies, abs(fourierTransform),color= 'Red', linewidth=1.0)
            plt.xlabel('Frequency')
            plt.ylabel('Normalized Amplitude (FFT of Accel)')
            plt.xlim([0, 6])
            amax=[frequencies[np.argmax(abs(fourierTransform))], max(abs(fourierTransform))]; plt.annotate(str(round(amax[0],3)) +"Hz", xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            

    else: 
        pb.start()
        pb.configure(maximum=1.05+ canvas.createRS.get()+ canvas.plotVel.get() + canvas.plotDisp.get() + canvas.createRS2.get() +canvas.arias.get())
        win.update_idletasks()
        tT = np.concatenate( (np.arange(0.05, 0.1, 0.005) , np.arange (0.1, 0.5, 0.01) , np.arange (0.5, 1, 0.02) , np.arange (1, float(canvas.entry_endPeriod.get()), 0.05) ) ) # Time vector for the spectral response
        freq = 1/tT # Frequenxy vector
        xi = float(canvas.entry_Damping.get()) # Damping factor
        df = 1.0/dtAccel1
        Sfin=[]
    
        noSubplotsRows = 1 + canvas.createRS.get()+ canvas.plotVel.get() + canvas.plotDisp.get() + canvas.createRS2.get() +canvas.arias.get();noSubplotsCols = 3;subplotCounter = 1
        yaxislimit = round(accelim(scaledAccel1, scaledAccel2, scaledAccel3)*1.1,2)
        nyaxislimit = 0.0 - yaxislimit
        plt.figure(figsize=(18,8 + 2*noSubplotsRows))

        plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
        plt.grid()
        plt.title(nameCh1)
        plt.xlabel('Time (secs)')
        plt.ylabel('Accel(g)')
        plt.plot(T1,scaledAccel1, label="Channel1", color= 'Red', linewidth=1.0)
        amax=maxaccel(scaledAccel1, T1); plt.annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(scaledAccel1, T1); plt.annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
        plt.ylim([nyaxislimit, yaxislimit])
        plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])

        subplotCounter+=1
        plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
        plt.grid()
        plt.title(nameCh2)
        plt.xlabel('Time (secs)')
        plt.ylabel('Accel(g)')
        plt.plot(T1,scaledAccel2, label="Channel2", color= 'Red', linewidth=1.0)
        amax=maxaccel(scaledAccel2, T2); plt.annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(scaledAccel2, T2); plt.annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
        plt.ylim([nyaxislimit, yaxislimit])
        plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])

        subplotCounter+=1
        plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
        plt.grid()
        plt.title(nameCh3)
        plt.xlabel('Time (secs)')
        plt.ylabel('Accel(g)')
        plt.plot(T1,scaledAccel3, label="Channel3", color= 'Red', linewidth=1.0)
        amax=maxaccel(scaledAccel3, T3); plt.annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        amin=minaccel(scaledAccel3, T3); plt.annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
        plt.ylim([nyaxislimit, yaxislimit])
        plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])
        
        if str(canvas.createRS2.get()) =="1":
            pb.step(1)
            win.update_idletasks()
            subplotCounter+=1
            ax =plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            Sfin= RS_function(accel1, df, tT, xi, Resp_type = 'SA')
            S=Sfin[0,:]*scaleValue(unitsAccel1)
            plt.xlabel('Period (secs)')
            plt.ylabel('SA (g)')
            ax.plot(tT,S,color= 'Red', linewidth=1.0)
            amax=[tT[np.argmax(abs(S))], max(abs(S))]; plt.annotate(str(round(amax[0],3)) +"sec, "+str(round(amax[1],2)) +"g" , xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            ax.text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)
            

            subplotCounter+=1
            ax=plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            Sfin= RS_function(accel2, df, tT, xi, Resp_type = 'SA')
            S=Sfin[0,:]*scaleValue(unitsAccel2)
            plt.xlabel('Period (secs)')
            plt.ylabel('SA (g)')
            ax.plot(tT,S,color= 'Red', linewidth=1.0)
            amax=[tT[np.argmax(abs(S))], max(abs(S))]; plt.annotate(str(round(amax[0],3)) +"sec, "+str(round(amax[1],2)) +"g" , xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            ax.text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)
            

            subplotCounter+=1
            ax=plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            Sfin= RS_function(accel3, df, tT, xi, Resp_type = 'SA')
            S=Sfin[0,:]*scaleValue(unitsAccel3)
            plt.xlabel('Period (secs)')
            plt.ylabel('SA (g)')
            ax.plot(tT,S,color= 'Red', linewidth=1.0)
            amax=[tT[np.argmax(abs(S))], max(abs(S))]; plt.annotate(str(round(amax[0],3)) +"sec, "+str(round(amax[1],2)) +"g" , xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            ax.text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)
            
        
        if str(canvas.createRS.get()) =="1":
            pb.step(1)
            win.update_idletasks()
            #canvas.Labelplot["text"]="This may take some time"
            subplotCounter+=1
            ax =plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            Sfin= RS_function(accel1, df, tT, xi, Resp_type = 'PSASD')
            S=Sfin[0,:]*scaleValue(unitsAccel1)
            area= round(np.trapezoid(Sfin[0,:],Sfin[1,:])/10000,2)
            #print(area)
            plt.xlabel('Peak D (cm)')
            plt.ylabel('Peak PSA (g)')
            ax.plot(Sfin[1,:],S,color= 'Red', linewidth=1.0)
            SfinClosed = np.append(np.insert(Sfin[1,:],0,0.0),0.0)
            SClosed = np.append(np.insert(S,0,0.0),0.0)
            plt.fill(SfinClosed,SClosed, "r", alpha=0.5)
            x_left, x_right = ax.get_xlim()
            y_low, y_high = ax.get_ylim()
            ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))
            #amax=maxSA(Sfin, Sfin2); plt.annotate(str(round(amax[0],3)) + " secs", xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            radialPeriods(1/scaleValue(unitsAccel1),plt, ax)
            ax.text(x_right/3, y_high/3, str(area) + r"$(m/s)^2$", horizontalalignment='center', fontsize=10, color ='Blue')
            ax.text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)

            subplotCounter+=1
            ax=plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            Sfin= RS_function(accel2, df, tT, xi, Resp_type = 'PSASD')
            S=Sfin[0,:]*scaleValue(unitsAccel2)
            area= round(np.trapezoid(Sfin[0,:],Sfin[1,:])/10000,2)
            #print(area)
            plt.xlabel('Peak D (cm)')
            plt.ylabel('Peak PSA (g)')
            ax.plot(Sfin[1,:],S,color= 'Red', linewidth=1.0)
            SfinClosed = np.append(np.insert(Sfin[1,:],0,0.0),0.0)
            SClosed = np.append(np.insert(S,0,0.0),0.0)
            plt.fill(SfinClosed,SClosed, "r", alpha=0.5)
            x_left, x_right = ax.get_xlim()
            y_low, y_high = ax.get_ylim()
            ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))
            #amax=maxSA(Sfin, Sfin2); plt.annotate(str(round(amax[0],3)) + " secs", xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            radialPeriods(1/scaleValue(unitsAccel2),plt, ax)
            ax.text(x_right/3, y_high/3, str(area) + r"$(m/s)^2$", horizontalalignment='center', fontsize=10, color ='Blue')
            ax.text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)
            
            subplotCounter+=1
            ax=plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            Sfin= RS_function(accel3, df, tT, xi, Resp_type = 'PSASD')
            S=Sfin[0,:]*scaleValue(unitsAccel3)
            area= round(np.trapezoid(Sfin[0,:],Sfin[1,:])/10000,2)
            #print(area)
            plt.xlabel('Peak D (cm)')
            plt.ylabel('Peak PSA (g)')
            SfinClosed = np.append(np.insert(Sfin[1,:],0,0.0),0.0)
            SClosed = np.append(np.insert(S,0,0.0),0.0)
            plt.fill(SfinClosed,SClosed, "r", alpha=0.5)
            ax.plot(Sfin[1,:],S,color= 'Red', linewidth=1.0)
            x_left, x_right = ax.get_xlim()
            y_low, y_high = ax.get_ylim()
            ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))
            #amax=maxSA(Sfin, Sfin2); plt.annotate(str(round(amax[0],3)) + " secs", xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            radialPeriods(1/scaleValue(unitsAccel3),plt, ax)
            ax.text(x_right/3, y_high/3, str(area) + r"$(m/s)^2$", horizontalalignment='center', fontsize=10, color ='Blue')
            ax.text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)

        if str(canvas.arias.get()) =="1":
            arias1 = np.cumsum(np.square(accel1)*dtAccel1*np.pi/2/980.665/100)
            arias2 = np.cumsum(np.square(accel2)*dtAccel2*np.pi/2/980.665/100)
            arias3 = np.cumsum(np.square(accel3)*dtAccel3*np.pi/2/980.665/100)
            normarias1 = arias1/np.max(arias1)
            normarias2 = arias2/np.max(arias2)
            normarias3 = arias3/np.max(arias3)
            


            ariasmax= max(np.max(arias1), np.max(arias2), np.max(arias3))

            subplotCounter+=1
            ax = plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Arias Intensity (m/s)')
            plt.plot(T1,arias1, label="Channel1", color= 'Red', linewidth=1.0)
            plt.ylim([0, ariasmax])
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])
            ax.text(0.97, 0.97, 'D5-75 = ' + str(round(T1[np.argmax(normarias1 > 0.75)] - T1[np.argmax(normarias1 > 0.05)],3)), horizontalalignment='right', verticalalignment='top', fontsize=9, color ='Black',transform=ax.transAxes)
            ax.text(0.97, 0.92, 'D5-95 = ' + str(round(T1[np.argmax(normarias1 > 0.95)] - T1[np.argmax(normarias1 > 0.05)],3)), horizontalalignment='right', verticalalignment='top', fontsize=9, color ='Black',transform=ax.transAxes)

            subplotCounter+=1
            ax = plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Arias Intensity (m/s)')
            plt.plot(T1,arias2, label="Channel1", color= 'Red', linewidth=1.0)
            plt.ylim([0, ariasmax])
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])
            ax.text(0.97, 0.97, 'D5-75 = ' + str(round(T1[np.argmax(normarias2 > 0.75)] - T1[np.argmax(normarias2 > 0.05)],3)), horizontalalignment='right', verticalalignment='top', fontsize=9, color ='Black',transform=ax.transAxes)
            ax.text(0.97, 0.92, 'D5-95 = ' + str(round(T1[np.argmax(normarias2 > 0.95)] - T1[np.argmax(normarias2 > 0.05)],3)), horizontalalignment='right', verticalalignment='top', fontsize=9, color ='Black',transform=ax.transAxes)

            subplotCounter+=1
            ax = plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Arias Intensity (m/s)')
            plt.plot(T1,arias3, label="Channel1", color= 'Red', linewidth=1.0)
            plt.ylim([0, ariasmax])
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])
            ax.text(0.97, 0.97, 'D5-75 = ' + str(round(T1[np.argmax(normarias3 > 0.75)] - T1[np.argmax(normarias3 > 0.05)],3)), horizontalalignment='right', verticalalignment='top', fontsize=9, color ='Black',transform=ax.transAxes)
            ax.text(0.97, 0.92, 'D5-95 = ' + str(round(T1[np.argmax(normarias3 > 0.95)] - T1[np.argmax(normarias3 > 0.05)],3)), horizontalalignment='right', verticalalignment='top', fontsize=9, color ='Black',transform=ax.transAxes)
    
        if str(canvas.plotVel.get()) =="1":
            pb.step(1)
            win.update_idletasks()
            subplotCounter+=1
            plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Velocity '+unitsVel1)
            plt.plot(T1,vel1, label="Channel1", color= 'Blue', linewidth=1.0)
            amax=T1[np.argmax(vel1)],np.max(vel1); plt.annotate(str(round(amax[1],1)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            amin=T1[np.argmin(vel1)],np.min(vel1); plt.annotate(str(round(amin[1],1)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])

            subplotCounter+=1
            plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Velocity '+unitsVel2)
            plt.plot(T1,vel2, label="Channel2", color= 'Blue', linewidth=1.0)
            amax=T1[np.argmax(vel2)],np.max(vel2); plt.annotate(str(round(amax[1],1)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            amin=T1[np.argmin(vel2)],np.min(vel2); plt.annotate(str(round(amin[1],1)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])

            subplotCounter+=1
            plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Velocity '+unitsVel3)
            plt.plot(T1,vel3, label="Channel3", color= 'Blue', linewidth=1.0)
            amax=T1[np.argmax(vel3)],np.max(vel3); plt.annotate(str(round(amax[1],1)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
            amin=T1[np.argmin(vel3)],np.min(vel3); plt.annotate(str(round(amin[1],1)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])
        
        if str(canvas.plotDisp.get()) =="1":
            pb.step(1)
            win.update_idletasks()
            subplotCounter+=1
            plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Disp '+unitsDispl1)
            plt.plot(T1,displ1, label="Channel1", color= 'Green', linewidth=1.0)
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])

            subplotCounter+=1
            plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Disp '+unitsDispl2)
            plt.plot(T1,displ2, label="Channel2", color= 'Green', linewidth=1.0)
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])

            subplotCounter+=1
            plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            plt.grid()
            plt.xlabel('Time (secs)')
            plt.ylabel('Disp '+unitsDispl3)
            plt.plot(T1,displ3, label="Channel3", color= 'Green', linewidth=1.0)
            plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])


       

    if EOF == 0 and str(canvas.plotOrbit.get()) =="1":    
        plt.figure(2,figsize=(14,10))
        noSubplotsRows = 1;noSubplotsCols = 1;subplotCounter = 1
        locanvasdex1=int(float(canvas.entry_Lowxlim.get())/dtDispl1); highIndex1=int(float(canvas.entry_Highxlim.get())/dtDispl1); 
        locanvasdex2=int(float(canvas.entry_Lowxlim.get())/dtDispl2); highIndex2=int(float(canvas.entry_Highxlim.get())/dtDispl2); 
        locanvasdex3=int(float(canvas.entry_Lowxlim.get())/dtDispl3); highIndex3=int(float(canvas.entry_Highxlim.get())/dtDispl3); 
        if str(canvas.plotOrbit.get()) =="1":
            ax = plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter, projection='3d')
            if "Up" in nameCh1:
                ax.set_xlabel(nameCh2 + " displacement", fontsize=6)
                ax.set_ylabel(nameCh3 + " displacement", fontsize=6)
                ax.set_zlabel(nameCh1 + " displacement", fontsize=6)
                ax.plot(displ2[locanvasdex2:highIndex2],displ3[locanvasdex3:highIndex3],displ1[locanvasdex1:highIndex1], label="Orbitplot", color= 'Black', linewidth=0.5)
            elif "Up" in nameCh2:
                ax.set_xlabel(nameCh1 + " displacement", fontsize=6)
                ax.set_ylabel(nameCh3 + " displacement", fontsize=6)
                ax.set_zlabel(nameCh2 + " displacement", fontsize=6)
                ax.plot(displ1[locanvasdex1:highIndex1],displ3[locanvasdex3:highIndex3],displ2[locanvasdex2:highIndex2], label="Orbitplot", color= 'Black', linewidth=0.5)
            elif "Up" in nameCh3:
                ax.set_xlabel(nameCh1 + " displacement", fontsize=6)
                ax.set_ylabel(nameCh2 + " displacement", fontsize=6)
                ax.set_zlabel(nameCh3 + " displacement", fontsize=6)
                ax.plot(displ1[locanvasdex1:highIndex1],displ2[locanvasdex2:highIndex2],displ3[locanvasdex3:highIndex3], label="Orbitplot", color= 'Black', linewidth=0.5)
    pb.stop()
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05)    
    plt.show()

def on_clickRot():
    #%% Parameters of the response spectra

    horRec=np.zeros((2,len(scaledAccel1)))
       
    if "Up" in nameCh1 or "HNZ" in nameCh1:
        if "360" in nameCh2 or "180" in nameCh2:
            horRec[0,:] = scaledAccel3
            horRec[1,:] = scaledAccel2
            horRec1=nameCh3;horRec2=nameCh2
        else:
            horRec[0,:] = scaledAccel2
            horRec[1,:] = scaledAccel3
            horRec1=nameCh2;horRec2=nameCh3

    elif "Up" in nameCh2 or "HNZ" in nameCh2:
        if "360" in nameCh1 or "180" in nameCh1:
            horRec[0,:] = scaledAccel3
            horRec[1,:] = scaledAccel1
            horRec1=nameCh3;horRec2=nameCh1
        else:
            horRec[0,:] = scaledAccel1
            horRec[1,:] = scaledAccel3
            horRec1=nameCh1;horRec2=nameCh3

    elif "Up" in nameCh3 or "HNZ" in nameCh3:
        if "360" in nameCh1 or "180" in nameCh1:
            horRec[0,:] = scaledAccel2
            horRec[1,:] = scaledAccel1
            horRec1=nameCh2;horRec2=nameCh1
        else:
            horRec[0,:] = scaledAccel1
            horRec[1,:] = scaledAccel2
            horRec1=nameCh1;horRec2=nameCh2
    
    if "360" in horRec2:
        horRec2 = horRec2.replace("360 Deg", "NS")
    elif "180" in horRec2:
        horRec[1,:]=[x*-1 for x in horRec[1,:]]
        horRec2 = horRec2.replace("180 Deg", "NS")
    
    if "90" in horRec1:
        horRec1 = horRec1.replace("90 Deg", "EW")
    elif "270" in horRec1:
        horRec[0,:]=[x*-1 for x in horRec[0,:]]
        horRec1 = horRec1.replace("270 Deg", "EW")

    plt.close(3)
    pb.start()
    pb.configure(maximum=1.05+(canvas.createRS.get()+canvas.createRS2.get()+canvas.createTrip.get()) )
    win.update_idletasks()
    plt.figure(3,figsize=(14,15))
    noSubplotsRows = 1 + math.ceil((canvas.createRS.get()+canvas.createRS2.get()+canvas.createTrip.get())/2);noSubplotsCols = 2;subplotCounter = 1

    ax = plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
    plt.title("Orbit plot for acceleration")
    plt.plot(horRec[0,:], horRec[1,:])
    rotmaxLoc = np.argmax(np.sqrt(np.square(horRec[0,:])+np.square(horRec[1,:])))
    resAccelmax = np.sqrt(np.square(horRec[0,rotmaxLoc])+np.square(horRec[1,rotmaxLoc]))
    resAngle = np.arctan2(horRec[1,rotmaxLoc],horRec[0,rotmaxLoc])
    #print(resAccelmax, horRec[0,rotmaxLoc]*np.cos(resAngle)+horRec[1,rotmaxLoc]*np.sin(resAngle) )
    plt.plot([0,horRec[0,rotmaxLoc]], [0, horRec[1,rotmaxLoc]], color='red',linewidth=2.0 )
    plt.annotate(round(resAccelmax,3), xy=(horRec[0,rotmaxLoc], horRec[1,rotmaxLoc]), xytext=(horRec[0,rotmaxLoc], horRec[1,rotmaxLoc]), fontsize=10, color= 'Blue')
    plt.xlabel(horRec1); plt.ylabel(horRec2)
    maxLimit = max(np.max(horRec), np.abs(np.min(horRec)))/0.95
    plt.xlim(-maxLimit, maxLimit)
    plt.ylim(-maxLimit, maxLimit)
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))

    resAccelmax = (horRec[0,:]*np.cos(resAngle)+horRec[1,:]*np.sin(resAngle))
    rotType="Resultant Acceleration in direction of maximum acceleration"
    rotatedplots(plt, ax, T1, resAccelmax, noSubplotsRows,noSubplotsCols, subplotCounter,rotType)

def rotatedplots(plt, ax, T1, resAccelmax, noSubplotsRows,noSubplotsCols, subplotCounter,rotType):
    if str(canvas.includeASCE.get())=="1":
        asceSpect = ASCE722Spectra()
    subplotCounter+=1
    plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
    plt.grid()
    plt.title(rotType)
    plt.xlabel('Time (secs)')
    plt.ylabel('Accel(g)')
    plt.plot(T1,resAccelmax , label="Resultant", color= 'Red', linewidth=1.0)
    amax=[T1[np.argmax(resAccelmax)], np.max(resAccelmax)]; plt.annotate(str(round(amax[1],3)), xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
    amin=[T1[np.argmin(resAccelmax)], np.min(resAccelmax)]; plt.annotate(str(round(amin[1],3)), xy=(amin[0], amin[1]), xytext=(amin[0], amin[1]), verticalalignment='top')
    plt.xlim([float(canvas.entry_Lowxlim.get()), float(canvas.entry_Highxlim.get())])
    
    resAccelmax = resAccelmax/scaleValue(unitsAccel1)
    tT = np.concatenate( (np.arange(0.05, 0.1, 0.005) , np.arange (0.1, 0.5, 0.01) , np.arange (0.5, 1, 0.02) , np.arange (1, float(canvas.entry_endPeriod.get()), 0.05) ) ) # Time vector for the spectral response
    freq = 1/tT # Frequenxy vector
    xi = float(canvas.entry_Damping.get()) # Damping factor
    df = 1.0/dtAccel1
    Sfin=[]
    if str(canvas.createRS2.get()) =="1":
        pb.step(1)
        win.update_idletasks()
        subplotCounter+=1
        ax =plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
        Sfin= RS_function(resAccelmax, df, tT, xi, Resp_type = 'SA')
        S=Sfin[0,:]*scaleValue(unitsAccel1)
        plt.xlabel('Period (secs)')
        plt.ylabel('SA (g)')
        ax.plot(tT,S,color= 'Red', linewidth=1.0, label = "Resultant Response Spectrum")
        amax=[tT[np.argmax(abs(S))], max(abs(S))]; plt.annotate(str(round(amax[0],3)) +"sec, "+str(round(amax[1],2)) +"g" , xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        if str(canvas.includeASCE.get())=="1":
            ax.plot(asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["periods"],\
                asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"],color= 'blue',linewidth=0.5,label="ASCE7-22 Multiperiod Spectrum")
            ax.plot(asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["periods"],\
                asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["ordinates"],color= 'green',linewidth=0.5,label="ASCE7-22 2-period Spectrum")
            ax.plot(asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["periods"],\
                asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"],color= 'blue', linestyle="--",linewidth=0.5,label="MCE Multiperiod Spectrum")
            ax.plot(asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["periods"],\
                asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["ordinates"],color= 'green',linestyle="--", linewidth=0.5,label="MCE 2-period Spectrum")
            ax.set_xlim(0,float(canvas.entry_endPeriod.get()))
        plt.legend(loc="center right",fontsize = 'x-small')
        ax.text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)

    
    if str(canvas.createRS.get()) =="1":
        #canvas.Labelplot["text"]="This may take some time"
        pb.step(1)
        win.update_idletasks()
        subplotCounter+=1
        ax =plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
        Sfin= RS_function(resAccelmax, df, tT, xi, Resp_type = 'PSASD')
        S=Sfin[0,:]*scaleValue(unitsAccel1)
        area= round(np.trapezoid(Sfin[0,:],Sfin[1,:])/10000,2)
        #print(area)
        plt.xlabel('Peak D (cm)')
        plt.ylabel('Peak PSA (g)')
        ax.plot(Sfin[1,:],S,color= 'Red', linewidth=1.0, label = "Resultant Response Spectrum")
        SfinClosed = np.append(np.insert(Sfin[1,:],0,0.0),0.0)
        SClosed = np.append(np.insert(S,0,0.0),0.0)
        plt.fill(SfinClosed,SClosed, "r", alpha=0.5)
        if str(canvas.includeASCE.get())=="1":
            mPS, tPS, MmPS, MtPS = convertADRS(asceSpect)
            ax.plot(mPS[0,:],mPS[1,:], color='blue', linewidth =0.5,label="ASCE7-22 Multiperiod Spectrum")
            ax.plot(tPS[0,:],tPS[1,:], color='green', linewidth =0.5,label="ASCE7-22 2-period Spectrum")
            ax.plot(MmPS[0,:],MmPS[1,:], color='blue',linestyle="--", linewidth =0.5,label="MCE Multiperiod Spectrum")
            ax.plot(MtPS[0,:],MtPS[1,:], color='green',linestyle="--", linewidth =0.5,label="MCE 2-period Spectrum")
            #plt.legend(loc="center right",fontsize = 'xx-small')
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()
        ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))
        #amax=maxSA(Sfin, Sfin2); plt.annotate(str(round(amax[0],3)) + " secs", xy=(amax[0], amax[1]), xytext=(amax[0], amax[1]))
        radialPeriods(1/scaleValue(unitsAccel1),plt, ax)
        ax.text(x_right/4, y_high/4, str(area) + r"$(m/s)^2$", horizontalalignment='center', fontsize=10, color ='Blue')
        ax.text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)
    
    if str(canvas.createTrip.get()) =="1":
        pb.step(1)
        win.update_idletasks()
        subplotCounter+=1
        ax =plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
        Sfin= RS_function(resAccelmax, df, tT, xi, Resp_type = 'SA')
        S=Sfin[0,:]/(2*np.pi/tT)
        plt.xlabel('Period (secs)')
        plt.ylabel('Psuedo Velocity '+ unitsVel1)
        plt.grid()
        ax.plot(tT,S,color= 'Red', linewidth=1.0, label = "Resultant Response Spectrum")

        x_left = np.min(tT); x_right =max(tT)
        if str(canvas.includeASCE.get())=="1":
            tasce =np.delete(np.array(asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["periods"]),0)
            aasce = (1/scaleValue(unitsAccel1)) *np.delete(np.array(asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"]),0)/(2*np.pi/tasce)
            ax.plot(tasce,aasce,color= 'blue',linewidth=1.0,label="ASCE7-22 Multiperiod Spectrum")
            tasce =np.delete(np.array(asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["periods"]),0)
            aasce = (1/scaleValue(unitsAccel1)) *np.delete(np.array(asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["ordinates"]),0)/(2*np.pi/tasce)
            ax.plot(tasce,aasce,color= 'green',linewidth=1.0,label="ASCE7-22 2-period Spectrum")
            tasce =np.delete(np.array(asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["periods"]),0)
            aasce = (1/scaleValue(unitsAccel1)) *np.delete(np.array(asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]),0)/(2*np.pi/tasce)
            ax.plot(tasce,aasce,color= 'blue', linestyle="--",linewidth=1.0,label="MCE Multiperiod Spectrum")
            tasce =np.delete(np.array(asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["periods"]),0)
            aasce = (1/scaleValue(unitsAccel1)) *np.delete(np.array(asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["ordinates"]),0)/(2*np.pi/tasce)
            ax.plot(tasce,aasce,color= 'green', linestyle="--",linewidth=1.0,label="MCE 2-period Spectrum")

        
        plt.xscale("log")
        plt.yscale("log")
        #plt.legend(loc="center",fontsize = 'x-small')
        ax.text(0.03, 0.90, 'Damping=' + str(round(xi,3)), horizontalalignment='left', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)
        tripartitegrids(1/scaleValue(unitsAccel1),plt,ax,x_left,x_right)
        ax.set_xlim(x_left,x_right)

    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05)    
    pb.stop()
    plt.show()

def on_clickRotAngle():
    #%% Parameters of the response spectra
    

    horRec=np.zeros((2,len(scaledAccel1)))

    if "Up" in nameCh1 or "HNZ" in nameCh1:
        if "360" in nameCh2 or "180" in nameCh2:
            horRec[0,:] = scaledAccel3
            horRec[1,:] = scaledAccel2
            horRec1=nameCh3;horRec2=nameCh2
        else:
            horRec[0,:] = scaledAccel2
            horRec[1,:] = scaledAccel3
            horRec1=nameCh2;horRec2=nameCh3

    elif "Up" in nameCh2 or "HNZ" in nameCh2:
        if "360" in nameCh1 or "180" in nameCh1:
            horRec[0,:] = scaledAccel3
            horRec[1,:] = scaledAccel1
            horRec1=nameCh3;horRec2=nameCh1
        else:
            horRec[0,:] = scaledAccel1
            horRec[1,:] = scaledAccel3
            horRec1=nameCh1;horRec2=nameCh3

    elif "Up" in nameCh3 or "HNZ" in nameCh3:
        if "360" in nameCh1 or "180" in nameCh1:
            horRec[0,:] = scaledAccel2
            horRec[1,:] = scaledAccel1
            horRec1=nameCh2;horRec2=nameCh1
        else:
            horRec[0,:] = scaledAccel1
            horRec[1,:] = scaledAccel2
            horRec1=nameCh1;horRec2=nameCh2
    
    if "360" in horRec2:
        horRec2 = horRec2.replace("360 Deg", "NS")
    elif "180" in horRec2:
        horRec[1,:]=[x*-1 for x in horRec[1,:]]
        horRec2 = horRec2.replace("180 Deg", "NS")
    
    if "90" in horRec1:
        horRec1 = horRec1.replace("90 Deg", "EW")
    elif "270" in horRec1:
        horRec[0,:]=[x*-1 for x in horRec[0,:]]
        horRec1 = horRec1.replace("270 Deg", "EW")

    plt.close(7)
    pb.start()
    pb.configure(maximum=1.05+(canvas.createRS.get()+canvas.createRS2.get()+canvas.createTrip.get()) )
    win.update_idletasks()
    plt.figure(7,figsize=(14,15))
    noSubplotsRows = 1 + math.ceil((canvas.createRS.get()+canvas.createRS2.get()+canvas.createTrip.get())/2);noSubplotsCols = 2;subplotCounter = 1
    ax = plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
    plt.title("Orbit plot for acceleration")
    plt.plot(horRec[0,:], horRec[1,:])
    
    #resAccelmax = np.sqrt(np.square(horRec[0,rotmaxLoc])+np.square(horRec[1,rotmaxLoc]))
    resAngle = float(canvas.entry_Angle.get())/180.0 * np.pi
    resAccelmax = (horRec[0,:]*np.cos(resAngle)+horRec[1,:]*np.sin(resAngle))
    rotmax = np.max(resAccelmax)
    #print(resAccelmax, horRec[0,rotmaxLoc]*np.cos(resAngle)+horRec[1,rotmaxLoc]*np.sin(resAngle) )
    plt.plot([0,rotmax*np.cos(resAngle)], [0, rotmax*np.sin(resAngle)], color='red',linewidth=2.0 )
    #plt.annotate(round(resAccelmax,3), xy=(horRec[0,rotmaxLoc], horRec[1,rotmaxLoc]), xytext=(horRec[0,rotmaxLoc], horRec[1,rotmaxLoc]), fontsize=10, color= 'Blue')
    plt.xlabel(horRec1); plt.ylabel(horRec2)
    maxLimit = max(np.max(horRec), np.abs(np.min(horRec)))/0.95
    plt.xlim(-maxLimit, maxLimit)
    plt.ylim(-maxLimit, maxLimit)
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))

    
    
    rotType="Resultant Acceleration in direction of specified angle"
    rotatedplots(plt, ax, T1, resAccelmax, noSubplotsRows,noSubplotsCols, subplotCounter,rotType)


def on_clickRotD50():
    #%% Parameters of the response spectra

    horRec=np.zeros((2,len(scaledAccel1)))
    res = messagebox.askyesno("RotD50","This is a computational intensive process \nand will take some time (3 to 10 mins)\nContinue?")
    if res:
        if "Up" in nameCh1 or "HNZ" in nameCh1:
            if "360" in nameCh2 or "180" in nameCh2:
                horRec[0,:] = scaledAccel3
                horRec[1,:] = scaledAccel2
                horRec1=nameCh3;horRec2=nameCh2
            else:
                horRec[0,:] = scaledAccel2
                horRec[1,:] = scaledAccel3
                horRec1=nameCh2;horRec2=nameCh3

        elif "Up" in nameCh2 or "HNZ" in nameCh2:
            if "360" in nameCh1 or "180" in nameCh1:
                horRec[0,:] = scaledAccel3
                horRec[1,:] = scaledAccel1
                horRec1=nameCh3;horRec2=nameCh1
            else:
                horRec[0,:] = scaledAccel1
                horRec[1,:] = scaledAccel3
                horRec1=nameCh1;horRec2=nameCh3

        elif "Up" in nameCh3 or "HNZ" in nameCh3:
            if "360" in nameCh1 or "180" in nameCh1:
                horRec[0,:] = scaledAccel2
                horRec[1,:] = scaledAccel1
                horRec1=nameCh2;horRec2=nameCh1
            else:
                horRec[0,:] = scaledAccel1
                horRec[1,:] = scaledAccel2
                horRec1=nameCh1;horRec2=nameCh2
        
        if "360" in horRec2:
            horRec2 = horRec2.replace("360 Deg", "NS")
        elif "180" in horRec2:
            horRec[1,:]=[x*-1 for x in horRec[1,:]]
            horRec2 = horRec2.replace("180 Deg", "NS")
        
        if "90" in horRec1:
            horRec1 = horRec1.replace("90 Deg", "EW")
        elif "270" in horRec1:
            horRec[0,:]=[x*-1 for x in horRec[0,:]]
            horRec1 = horRec1.replace("270 Deg", "EW")

        plt.close(5)
        pb.start()
        pb.configure(maximum=182)
        if str(canvas.includeASCE.get())=="1":
            asceSpect = ASCE722Spectra()
        plt.figure(6,figsize=(14,15))
        noSubplotsRows = 2 ;noSubplotsCols = 1 + canvas.createTrip.get();subplotCounter = 1
        ax = plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)


        tT = np.concatenate( (np.arange(0.05, 0.1, 0.01) , np.arange (0.1, 0.5, 0.05) , np.arange (0.5, 1, 0.05) , np.arange (1, float(canvas.entry_endPeriod.get()), 0.1) ) ) # Time vector for the spectral response
        freq = 1/tT # Frequenxy vector
        xi = float(canvas.entry_Damping.get()) # Damping factor
        df = 1.0/dtAccel1
        Sfin=[]
        rotmax = np.zeros((180,len(tT)))
        rotmaxlimit = np.zeros((360,2))
        for i in range(0,180,1):
            resAngle = i/180.0 * np.pi
            resAccel = np.zeros((len(horRec[0,:])))
            resAccel = (horRec[0,:]*np.cos(resAngle)+horRec[1,:]*np.sin(resAngle))
            rotmaxlimit[i,:] = absmaxND(resAccel)*np.cos(resAngle), absmaxND(resAccel)*np.sin(resAngle)
            rotmaxlimit[i+180,:] = -absmaxND(resAccel)*np.cos(resAngle), -absmaxND(resAccel)*np.sin(resAngle)
            Sfin= RS_function(resAccel, df, tT, xi, Resp_type = 'SA')
            rotmax[i,:]= Sfin[0,:]
            pb.step(1)
            win.update_idletasks()
        
        rotD50Spec = np.zeros((len(tT)));rotD100Spec = np.zeros((len(tT)));rotD00Spec = np.zeros((len(tT)))
        for i in range(0,len(tT),1):
            rotD50Spec[i] = np.median(rotmax[:,i])
            rotD100Spec[i] = np.max(rotmax[:,i])
            rotD00Spec[i] = np.min(rotmax[:,i])


        plt.title("Orbit plot for acceleration")
        plt.plot(horRec[0,:], horRec[1,:])
        plt.plot(np.append(rotmaxlimit[:,0],rotmaxlimit[-1,0]), np.append(rotmaxlimit[:,1],rotmaxlimit[-1,1]))
   
        plt.xlabel(horRec1); plt.ylabel(horRec2)
        maxLimit = max(np.max(horRec), np.abs(np.min(horRec)))/0.95
        plt.xlim(-maxLimit, maxLimit)
        plt.ylim(-maxLimit, maxLimit)
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()
        ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))


        Sfin1= RS_function(horRec[0,:], df, tT, xi, Resp_type = 'SA')
        Sfin2= RS_function(horRec[1,:], df, tT, xi, Resp_type = 'SA')
        geomeanSpectra = np.sqrt(np.array(Sfin1[0:])*np.array(Sfin2[0,:]))

        subplotCounter+=1
        ax =plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
        plt.xlabel('Period (secs)')
        plt.ylabel('SA (g)')
        ax.plot(tT,rotD50Spec,color= 'Red', linewidth=1.0, label = "RotD50 Response Spectrum")
        ax.plot(tT,rotD100Spec,color= 'Red', linestyle="--", linewidth=1.0, label = "RotD100 Response Spectrum")
        ax.plot(tT,rotD00Spec,color= 'Red', linestyle='-.', linewidth=1.0, label = "RotD00 Response Spectrum")
        ax.plot(tT,geomeanSpectra[0,:],color= 'k', linewidth=1.0, label = "Geomean Spectra")
        if str(canvas.includeASCE.get())=="1":
            ax.plot(asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["periods"],\
                asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"],color= 'blue',linewidth=0.5,label="ASCE7-22 Multiperiod Spectrum")
            ax.plot(asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["periods"],\
                asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["ordinates"],color= 'green',linewidth=0.5,label="ASCE7-22 2-period Spectrum")
            ax.plot(asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["periods"],\
                asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"],color= 'blue', linestyle="--",linewidth=0.5,label="MCE Multiperiod Spectrum")
            ax.plot(asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["periods"],\
                asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["ordinates"],color= 'green',linestyle="--", linewidth=0.5,label="MCE 2-period Spectrum")
            ax.set_xlim(0,float(canvas.entry_endPeriod.get()))
        plt.legend(loc="center right",fontsize = 'x-small')
        ax.text(0.97, 0.97, 'Damping=' + str(round(xi,3)), horizontalalignment='right', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)

        if str(canvas.createTrip.get()) =="1":
            subplotCounter+=1
            ax =plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
            S=rotD50Spec[:]*(1/scaleValue(unitsAccel1))/(2*np.pi/tT)
            plt.xlabel('Period (secs)')
            plt.ylabel('Psuedo Velocity cm/sec')
            plt.grid()
            ax.plot(tT,S,color= 'Red', linewidth=1.0, label = "RotD50 Response Spectrum")

            x_left = np.min(tT); x_right =max(tT)
            if str(canvas.includeASCE.get())=="1":
                tasce =np.delete(np.array(asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["periods"]),0)
                aasce = (1/scaleValue(unitsAccel1)) *np.delete(np.array(asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"]),0)/(2*np.pi/tasce)
                ax.plot(tasce,aasce,color= 'blue',linewidth=1.0,label="ASCE7-22 Multiperiod Spectrum")
                tasce =np.delete(np.array(asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["periods"]),0)
                aasce = (1/scaleValue(unitsAccel1)) *np.delete(np.array(asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["ordinates"]),0)/(2*np.pi/tasce)
                ax.plot(tasce,aasce,color= 'green',linewidth=1.0,label="ASCE7-22 2-period Spectrum")
                tasce =np.delete(np.array(asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["periods"]),0)
                aasce = (1/scaleValue(unitsAccel1)) *np.delete(np.array(asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]),0)/(2*np.pi/tasce)
                ax.plot(tasce,aasce,color= 'blue', linestyle="--",linewidth=1.0,label="MCE Multiperiod Spectrum")
                tasce =np.delete(np.array(asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["periods"]),0)
                aasce = (1/scaleValue(unitsAccel1)) *np.delete(np.array(asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["ordinates"]),0)/(2*np.pi/tasce)
                ax.plot(tasce,aasce,color= 'green', linestyle="--",linewidth=1.0,label="MCE 2-period Spectrum")

        
            plt.xscale("log")
            plt.yscale("log")
            #plt.legend(loc="center",fontsize = 'x-small')
            ax.text(0.03, 0.90, 'Damping=' + str(round(xi,3)), horizontalalignment='left', verticalalignment='top', fontsize=6, color ='Black',transform=ax.transAxes)
            tripartitegrids(1/scaleValue(unitsAccel1),plt,ax,x_left,x_right)
            ax.set_xlim(x_left,x_right)
        plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05)    
        pb.stop()
        
        if str(canvas.createAz3D.get()) =="1":
            plt.figure(7)
            surf = np.zeros((3, 180*len(tT)))
            k=0
            for i in range(0,180,1):
                for j in range(0, len(tT)):
                    surf[2,k] = rotmax[i,j]
                    surf[0,k] = i
                    surf[1,k] = tT[j] 
                    k+=1
            zmax = np.max(surf[2,:])
            ax = plt.subplot(1,1,1, projection='3d')
            my_cmap = plt.get_cmap('hot')
            ax.set_box_aspect(aspect = (1,2,1))
            surfplot=ax.plot_trisurf(surf[0,:],surf[1,:],surf[2,:], cmap = my_cmap, edgecolor ='none')
            plt.figure(7).colorbar(surfplot, shrink = 0.5, aspect = 5)
            ax.tricontour(surf[0,:],surf[1,:],surf[2,:], zdir='z', offset = -1, cmap = my_cmap)
            ax.set(zlim=(-1,zmax))
            ax.set_xlabel("Azimuth Angle", fontsize=8)
            ax.set_ylabel("Period (secs)", fontsize=8)
            ax.set_zlabel("Spectra Acceleration (g)", fontsize=8)
        plt.show()
 
def on_clickRotDisp():
    #%% Parameters of the response spectra

    horRec=np.zeros((2,len(scaledAccel1)))
    horDisp=np.zeros((2,len(displ1)))
  
    if "Up" in nameCh1 or "HNZ" in nameCh1:
        if "360" in nameCh2 or "180" in nameCh2:
            horRec[0,:] = scaledAccel3; horDisp[0,:] = displ3
            horRec[1,:] = scaledAccel2; horDisp[1,:] = displ2
            horRec1=nameCh3;horRec2=nameCh2
        else:
            horRec[0,:] = scaledAccel2; horDisp[0,:] = displ2
            horRec[1,:] = scaledAccel3; horDisp[1,:] = displ3
            horRec1=nameCh2;horRec2=nameCh3

    elif "Up" in nameCh2 or "HNZ" in nameCh2:
        if "360" in nameCh1 or "180" in nameCh1:
            horRec[0,:] = scaledAccel3; horDisp[0,:] = displ3
            horRec[1,:] = scaledAccel1; horDisp[1,:] = displ1
            horRec1=nameCh3;horRec2=nameCh1
        else:
            horRec[0,:] = scaledAccel1; horDisp[0,:] = displ1
            horRec[1,:] = scaledAccel3; horDisp[1,:] = displ3
            horRec1=nameCh1;horRec2=nameCh3

    elif "Up" in nameCh3 or "HNZ" in nameCh3:
        if "360" in nameCh1 or "180" in nameCh1:
            horRec[0,:] = scaledAccel2; horDisp[0,:] = displ2
            horRec[1,:] = scaledAccel1; horDisp[1,:] = displ1
            horRec1=nameCh2;horRec2=nameCh1
        else:
            horRec[0,:] = scaledAccel1; horDisp[0,:] = displ1
            horRec[1,:] = scaledAccel2; horDisp[1,:] = displ2
            horRec1=nameCh1;horRec2=nameCh2
    
    if "360" in horRec2:
        horRec2 = horRec2.replace("360 Deg", "NS")
    elif "180" in horRec2:
        horRec[1,:]=[x*-1 for x in horRec[1,:]];horDisp[1,:]=[x*-1 for x in horDisp[1,:]]
        horRec2 = horRec2.replace("180 Deg", "NS")
    
    if "90" in horRec1:
        horRec1 = horRec1.replace("90 Deg", "EW")
    elif "270" in horRec1:
        horRec[0,:]=[x*-1 for x in horRec[0,:]];horDisp[0,:]=[x*-1 for x in horDisp[0,:]]
        horRec1 = horRec1.replace("270 Deg", "EW")

    plt.close(4)
    pb.start()
    pb.configure(maximum=1.05+(canvas.createRS.get()+canvas.createRS2.get()+canvas.createTrip.get()) )
    win.update_idletasks()
    plt.figure(4,figsize=(14,15))
    noSubplotsRows = 1 + math.ceil((canvas.createRS.get()+canvas.createRS2.get()+canvas.createTrip.get())/2);noSubplotsCols = 2;subplotCounter = 1
    ax = plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
    plt.title("Orbit plot for displacement")
    plt.plot(horDisp[0,:], horDisp[1,:])
    rotmaxLoc = np.argmax(np.sqrt(np.square(horDisp[0,:])+np.square(horDisp[1,:])))
    resDispmax = np.sqrt(np.square(horDisp[0,rotmaxLoc])+np.square(horDisp[1,rotmaxLoc]))
    resAngle = np.arctan2(horDisp[1,rotmaxLoc],horDisp[0,rotmaxLoc])
    #print(resAccelmax, horRec[0,rotmaxLoc]*np.cos(resAngle)+horRec[1,rotmaxLoc]*np.sin(resAngle) )
    plt.plot([0,horDisp[0,rotmaxLoc]], [0, horDisp[1,rotmaxLoc]], color='red',linewidth=2.0 )
    plt.annotate(round(resDispmax,3), xy=(horDisp[0,rotmaxLoc], horDisp[1,rotmaxLoc]), xytext=(horDisp[0,rotmaxLoc], horDisp[1,rotmaxLoc]), fontsize=10, color= 'Blue')
    plt.xlabel(horRec1); plt.ylabel(horRec2)
    maxLimit = max(np.max(horDisp), np.abs(np.min(horDisp)))/0.95
    plt.xlim(-maxLimit, maxLimit)
    plt.ylim(-maxLimit, maxLimit)
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))

    resAccelmax = (horRec[0,:]*np.cos(resAngle)+horRec[1,:]*np.sin(resAngle))
    
    rotType="Resultant Acceleration in direction of maximum displacement"
    rotatedplots(plt, ax, T1, resAccelmax, noSubplotsRows,noSubplotsCols, subplotCounter,rotType)

def on_clickRotVel():
    #%% Parameters of the response spectra

    horRec=np.zeros((2,len(scaledAccel1)))
    horVel=np.zeros((2,len(displ1)))
  
    if "Up" in nameCh1 or "HNZ" in nameCh1:
        if "360" in nameCh2 or "180" in nameCh2:
            horRec[0,:] = scaledAccel3; horVel[0,:] = vel3
            horRec[1,:] = scaledAccel2; horVel[1,:] = vel2
            horRec1=nameCh3;horRec2=nameCh2
        else:
            horRec[0,:] = scaledAccel2; horVel[0,:] = vel2
            horRec[1,:] = scaledAccel3; horVel[1,:] = vel3
            horRec1=nameCh2;horRec2=nameCh3

    elif "Up" in nameCh2 or "HNZ" in nameCh2:
        if "360" in nameCh1 or "180" in nameCh1:
            horRec[0,:] = scaledAccel3; horVel[0,:] = vel3
            horRec[1,:] = scaledAccel1; horVel[1,:] = vel1
            horRec1=nameCh3;horRec2=nameCh1
        else:
            horRec[0,:] = scaledAccel1; horVel[0,:] = vel1
            horRec[1,:] = scaledAccel3; horVel[1,:] = vel3
            horRec1=nameCh1;horRec2=nameCh3

    elif "Up" in nameCh3 or "HNZ" in nameCh3:
        if "360" in nameCh1 or "180" in nameCh1:
            horRec[0,:] = scaledAccel2; horVel[0,:] = vel2
            horRec[1,:] = scaledAccel1; horVel[1,:] = vel1
            horRec1=nameCh2;horRec2=nameCh1
        else:
            horRec[0,:] = scaledAccel1; horVel[0,:] = vel1
            horRec[1,:] = scaledAccel2; horVel[1,:] = vel2
            horRec1=nameCh1;horRec2=nameCh2
    
    if "360" in horRec2:
        horRec2 = horRec2.replace("360 Deg", "NS")
    elif "180" in horRec2:
        horRec[1,:]=[x*-1 for x in horRec[1,:]];horVel[1,:]=[x*-1 for x in horVel[1,:]]
        horRec2 = horRec2.replace("180 Deg", "NS")
    
    if "90" in horRec1:
        horRec1 = horRec1.replace("90 Deg", "EW")
    elif "270" in horRec1:
        horRec[0,:]=[x*-1 for x in horRec[0,:]];horVel[0,:]=[x*-1 for x in horVel[0,:]]
        horRec1 = horRec1.replace("270 Deg", "EW")

    plt.close(5)
    pb.start()
    pb.configure(maximum=1.05+(canvas.createRS.get()+canvas.createRS2.get()+canvas.createTrip.get()) )
    win.update_idletasks()

    plt.figure(5,figsize=(14,15))
    noSubplotsRows = 1 + math.ceil((canvas.createRS.get()+canvas.createRS2.get()+canvas.createTrip.get())/2);noSubplotsCols = 2;subplotCounter = 1
    ax = plt.subplot(noSubplotsRows,noSubplotsCols,subplotCounter)
    plt.title("Orbit plot for Velocity")
    plt.plot(horVel[0,:], horVel[1,:])
    rotmaxLoc = np.argmax(np.sqrt(np.square(horVel[0,:])+np.square(horVel[1,:])))
    resVelmax = np.sqrt(np.square(horVel[0,rotmaxLoc])+np.square(horVel[1,rotmaxLoc]))
    resAngle = np.arctan2(horVel[1,rotmaxLoc],horVel[0,rotmaxLoc])
    #print(resAccelmax, horRec[0,rotmaxLoc]*np.cos(resAngle)+horRec[1,rotmaxLoc]*np.sin(resAngle) )
    plt.plot([0,horVel[0,rotmaxLoc]], [0, horVel[1,rotmaxLoc]], color='red',linewidth=2.0 )
    plt.annotate(round(resVelmax,3), xy=(horVel[0,rotmaxLoc], horVel[1,rotmaxLoc]), xytext=(horVel[0,rotmaxLoc], horVel[1,rotmaxLoc]), fontsize=10, color= 'Blue')
    plt.xlabel(horRec1); plt.ylabel(horRec2)
    maxLimit = max(np.max(horVel), np.abs(np.min(horVel)))/0.95
    plt.xlim(-maxLimit, maxLimit)
    plt.ylim(-maxLimit, maxLimit)
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high)))

    resAccelmax = (horRec[0,:]*np.cos(resAngle)+horRec[1,:]*np.sin(resAngle))
    rotType="Resultant Acceleration in direction of maximum Velocity"
    rotatedplots(plt, ax, T1, resAccelmax, noSubplotsRows,noSubplotsCols, subplotCounter,rotType)

def onclick2():
    plt.close('all')
    win.destroy()

def saveFile():
    messagebox.showinfo('Save Folder', 'Select folder to save in :\n Format:Time<space>Acceleration')
    dir_name = fd.askdirectory()
    j=0
    saveCh1value = str(canvas.SaveCh1.get())
    saveCh2value = str(canvas.SaveCh2.get())
    saveCh3value = str(canvas.SaveCh3.get())
    #print(saveCh1value)
    if saveCh1value == "1":
        index = len(T1)
        #print(index)
        with open(dir_name +'/Ch1.txt', 'w') as f:
            while j < index:
                f.write(str(round(T1[j],3))+ " " + str(scaledAccel1[j])+"\n")
                j+= 1
        #print(j)
    j=0
    if saveCh2value == "1":
        index = len(T2)
        with open(dir_name +'/Ch2.txt', 'w') as f:
            while j < index:
                f.write(str(round(T2[j],3))+ " " + str(scaledAccel2[j])+"\n")
                j+= 1
    j=0
    if saveCh3value == "1":
        index = len(T3)
        with open(dir_name +'/Ch3.txt', 'w') as f:
            while j < index:
                f.write(str(round(T3[j],))+ " " + str(scaledAccel3[j])+"\n")
                j+= 1

def radialPeriods(scale, plt, ax):
    periodSeries = np.concatenate(( np.arange(0.1,1.0,0.1) , np.arange(1.0,2.0,0.5), np.arange(2.0,5.0,1) ))
    #print(periodSeries)

    dispLimit, AccelLimit = ax.transData.inverted().transform(ax.transAxes.transform((0.95,0.95)))

    w = 2*np.pi/periodSeries
    w2 = np.square(w)
    a1 = dispLimit*w2/scale
    d1 = a1*scale/w2
    d2 = AccelLimit*scale/w2
    a2 = d2*w2/scale

    for i, items in enumerate(d1):
        if a1[i] < AccelLimit:
            plt.plot([0, d1[i]],[0, a1[i]], linestyle="--", color= 'Blue',linewidth=0.4)
            plt.annotate(round(periodSeries[i],1), xy=(d1[i], a1[i]), xytext=(d1[i], a1[i]), fontsize=5, color= 'Blue')
        else:
            plt.plot([0, d2[i]],[0, a2[i]], linestyle="--", color= 'Blue',linewidth=0.4)
            plt.annotate(round(periodSeries[i],1), xy=(d2[i], a2[i]), xytext=(d2[i], a2[i]), fontsize=5, color = 'Blue')

def tripartitegrids(scale,plt,ax,xl,xr):
    aSeries = np.array([0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100])
    dSeries = np.array([0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100])
    bset= np.array([0.01,0.1,1,10,100])
    bset2 = np.array([0.06,0.08,0.09,0.6,0.8,0.9,6,7,8,9,60,80,90])
    periodLimit, velLimit =  ax.transData.inverted().transform(ax.transAxes.transform((0.95,0.95)))
    periodLimit0, velLimit0 =  ax.transData.inverted().transform(ax.transAxes.transform((0.0,0.0)))
   

    for i, items in enumerate(aSeries):
        t0 =2*np.pi*velLimit0/(aSeries[i]*scale)
        t1= 2*np.pi*velLimit/(aSeries[i]*scale)
        
        t=t1;v=velLimit
        m = str(aSeries[i])+ "g"
        if aSeries[i] in bset:
            plt.plot([t0,t],[velLimit0,velLimit], linestyle="--", color= 'k',linewidth=0.3)
            if t < 0.95*xr and t > 1.05*xl:
                plt.annotate(m, xy=(t,v), xytext=(t,v), fontsize=5, color= 'k')
        elif aSeries[i] in bset2:
            plt.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'c',linewidth=0.3)
        else:
            plt.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'c',linewidth=0.3)
            if t < 0.95*xr and t > 1.05*xl:
                plt.annotate(m, xy=(t,v), xytext=(t,v), fontsize=5, color= 'c')

    for i, items in enumerate(dSeries):
        t0 =2*np.pi*dSeries[i]/velLimit0
        t1= 2*np.pi*dSeries[i]/velLimit
        t=t0; v=velLimit0
        m= str(dSeries[i])+"cm"
        if dSeries[i] in bset:        
            plt.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'k',linewidth=0.3)
            if t < 0.95*xr and t > 1.05*xl:
                plt.annotate(m, xy=(t,v), xytext=(t,v), ha='left', va="top",fontsize=5, color= 'k')
        elif dSeries[i] in bset2:
            plt.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'm',linewidth=0.3)
        else:
            plt.plot([t0,t1],[velLimit0,velLimit], linestyle="--", color= 'm',linewidth=0.3)
            if t < 0.95*xr and t > 1.05*xl:
                plt.annotate(m, xy=(t,v), xytext=(t,v), ha='left', va="top",fontsize=5, color= 'm')

def myabout():
        messagebox.showinfo('ReadV2', 'Read COSMOS V2 files, plot and analyze instrument records of earthquakes \nWritten by HXB')

def startlimAccel():
    a1 = next(i for i, x in enumerate(accel1) if abs(x) >1)
    startTime  = max(a1*dtAccel1- 2, 0)
    if EOF==0:
        a2 = next(i for i, x in enumerate(accel2) if abs(x) >1)
        a3 = next(i for i, x in enumerate(accel3) if abs(x) >1)
        startTime  = max(min(a1*dtAccel1,a2*dtAccel2,a3*dtAccel3) - 2, 0)
    return round(startTime,2)

def endlimAccel():
    a1 = next(i for i, x in reversed(list(enumerate(accel1))) if abs(x) >1)
    endTime  = max(a1*dtAccel1+ 2, 0)
    if EOF==0:
        a2 = next(i for i, x in reversed(list(enumerate(accel2))) if abs(x) >1)
        a3 = next(i for i, x in reversed(list(enumerate(accel3))) if abs(x) >1)
        endTime  = max(min(a1*dtAccel1,a2*dtAccel2,a3*dtAccel3) +2, 0)
    return round(endTime,2)

def ASCE722Spectra():
    lat = str(canvas.entry_Lat.get())
    longt= str(canvas.entry_Long.get())
    riskct = str(canvas.SelectedRiskCategory.get())
    sitecl = str(canvas.SelectedSiteClass.get())
    url = 'https://earthquake.usgs.gov/ws/designmaps/asce7-22.json?latitude='+ lat + '&longitude=' + longt\
 +'&riskCategory='+ riskct +'&siteClass=' + sitecl + '&title=Example'
    #print(url)
    try:
        response = ur.urlopen(url)
    except ur.error.HTTPError as e:
        if hasattr(e, 'reason'):
            print('We failed to reach a server.')
            print('Reason: ', e.reason)
            return()
        elif hasattr(e, 'code'):
            print('The server couldn\'t fulfill the request.')
            print('Error code: ', e.code)
            return()

    rdata = js.loads(response.read())


    #output = 'Output for Latitude = ' + str(latitude) + ' Longitude = ' + str(longitude

    return rdata
    
def convertADRS(asceSpect):
    multiPeriodSpectra = np.zeros((2,len(asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["periods"])))
    multiPeriodSpectra[0,:] = asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["periods"]
    multiPeriodSpectra[1,:] = asceSpect["response"]["data"]["multiPeriodDesignSpectrum"]["ordinates"]
    if multiPeriodSpectra[0,-1] > float(canvas.entry_endPeriod.get()):
        a1 =next(i for i,x in enumerate(multiPeriodSpectra[0,:]) if x >= float(canvas.entry_endPeriod.get()))+1
    else:
        a1 = len(multiPeriodSpectra[0,:])
    mPS=np.zeros((2,a1-1))
    mPS[0,:] = (1/scaleValue(unitsAccel1))*multiPeriodSpectra[1,1:a1]/np.square(2*np.pi/multiPeriodSpectra[0,1:a1])
    mPS[1,:] = multiPeriodSpectra[1,1:a1]

    MCEmultiPeriodSpectra = np.zeros((2,len(asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["periods"])))
    MCEmultiPeriodSpectra[0,:] = asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["periods"]
    MCEmultiPeriodSpectra[1,:] = asceSpect["response"]["data"]["multiPeriodMCErSpectrum"]["ordinates"]
    if MCEmultiPeriodSpectra[0,-1] > float(canvas.entry_endPeriod.get()):
        a1 =next(i for i,x in enumerate(MCEmultiPeriodSpectra[0,:]) if x >= float(canvas.entry_endPeriod.get()))+1
    else:
        a1 = len(MCEmultiPeriodSpectra[0,:])
    MmPS=np.zeros((2,a1-1))
    MmPS[0,:] = (1/scaleValue(unitsAccel1))*MCEmultiPeriodSpectra[1,1:a1]/np.square(2*np.pi/multiPeriodSpectra[0,1:a1])
    MmPS[1,:] = MCEmultiPeriodSpectra[1,1:a1]

    twopSpectra = np.zeros((2,len(asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["periods"])))
    twopSpectra[0,:] = asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["periods"]
    twopSpectra[1,:] = asceSpect["response"]["data"]["twoPeriodDesignSpectrum"]["ordinates"]
    if twopSpectra[0,-1] > float(canvas.entry_endPeriod.get()):
        a2 =next(i for i,x in enumerate(twopSpectra[0,:]) if x >= float(canvas.entry_endPeriod.get()))+1
    else:
        a2 = len(twopSpectra[0,:])
    tPS=np.zeros((2,a2-1))
    tPS[0,:] = (1/scaleValue(unitsAccel1))*twopSpectra[1,1:a2]/np.square(2*np.pi/twopSpectra[0,1:a2])
    tPS[1,:] =twopSpectra[1,1:a2]
    
    MCEtwopSpectra = np.zeros((2,len(asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["periods"])))
    MCEtwopSpectra[0,:] = asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["periods"]
    MCEtwopSpectra[1,:] = asceSpect["response"]["data"]["twoPeriodMCErSpectrum"]["ordinates"]
    if MCEtwopSpectra[0,-1] > float(canvas.entry_endPeriod.get()):
        a2 =next(i for i,x in enumerate(MCEtwopSpectra[0,:]) if x >= float(canvas.entry_endPeriod.get()))+1
    else:
        a2 = len(MCEtwopSpectra[0,:])
    MtPS=np.zeros((2,a2-1))
    MtPS[0,:] = (1/scaleValue(unitsAccel1))*MCEtwopSpectra[1,1:a2]/np.square(2*np.pi/MCEtwopSpectra[0,1:a2])
    MtPS[1,:] =MCEtwopSpectra[1,1:a2]

    return mPS, tPS, MmPS, MtPS

def readFile():
    global filenames
    global latitude, longitude
    global nameCh1,nameCh2,nameCh3,dtAccel1,dtAccel2,dtAccel3,dtDispl1,dtDispl2,dtDispl3,dtAccel1,dtVel2,dtVel3
    global numofPointsAccel1, numofPointsAccel2, numofPointsAccel3
    global T1,T2,T3
    global scaledAccel1, scaledAccel2, scaledAccel3, accel1, accel2, accel3, vel1, vel2, vel3,displ1, displ2, displ3
    global EOF
    global unitsAccel1, unitsAccel2, unitsAccel3, unitsVel1, unitsVel2, unitsVel3, unitsDispl1, unitsDispl2, unitsDispl3

    messagebox.showinfo('ReadV2', 'Select COSMOS V2 file (freefield or single channel record)')
    filetypes = (
            ('text files', '*.v2'),('text files', '*.V2'),
            ('All files', '*.*')
        )
    filenames = fd.askopenfilenames(
            title='Open files',
            initialdir='./',
            filetypes=filetypes)

    if len(filenames)==0:
        messagebox.showinfo('Error', 'File not selected, exiting')
        exit()
        
    #print(filenames)
    f=open(filenames[0])
    #f=open("./ReadV2/ce89781.v2")
        # islice (f,start, none)

    for line in islice(f, 5, 6):    
        latlong= line[17:40].strip()
        latitude =float(latlong[:latlong.find(",")-1])
        if latlong[len(latlong)-1:len(latlong)]=="W":
            longitude =float("-"+latlong[latlong.find(",")+1: len(latlong)-1].strip())
        else:
            longitude =float(latlong[latlong.find(",")+1: len(latlong)-1].strip())
        #print(latitude, longitude)
    for line in islice(f, 17, 18):
        nameCh1=line[26:].strip()
    for line in islice(f, 0, 1):
        nameCh1=nameCh1 + line.strip()
        #print(nameCh1)

    for line in islice(f, 20, 21):
        #print(line)
        numofPointsAccel1 = int(line[0: line.find("points")].strip())
        dtAccel1 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
        unitsAccel1 = line[line.find(", in") + 4: line.find(". (")].strip()
    numofLines = lines(numofPointsAccel1)
    accel1 = readchunk(f,numofLines)

    for line in islice(f, 0,1):
        #print(line)
        numofPointsVel1 = int(line[0: line.find("points")].strip())
        dtVel1 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
        unitsVel1 = line[line.find(", in") + 4: line.find(".  (")].strip()
    numofLines = lines(numofPointsVel1)
    vel1 = readchunk(f,numofLines)

    for line in islice(f, 0,1):
        #print(line)
        numofPointsDispl1 = int(line[0: line.find("points")].strip())
        dtDispl1 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
        unitsDispl1 = line[line.find(", in") + 4: line.find(".   ")].strip()
    numofLines = lines(numofPointsDispl1)
    displ1 = readchunk(f,numofLines)

    T1 = np.arange(0.0,numofPointsAccel1*dtAccel1, dtAccel1)
    scale = scaleValue(unitsAccel1) 
    scaledAccel1 = [value*scale for value in accel1]

    EOF=1
    for line in islice(f,1,24):
        if line != "":
            EOF = 0

    for line in islice(f, 0, 1):
        nameCh2=line[26:].strip()
    for line in islice(f, 0, 1):
        nameCh2=nameCh2 + line.strip()
        #print(nameCh2)

    for line in islice(f,1,20):
        i=0

    if EOF ==0:
        for line in islice(f, 0, 1):
            #print(line)
            numofPointsAccel2 = int(line[0: line.find("points")].strip())
            dtAccel2 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
            unitsAccel2 = line[line.find(", in") + 4: line.find(". (")].strip()
        numofLines = lines(numofPointsAccel2)
        accel2 = readchunk(f,numofLines)

        for line in islice(f, 0,1):
            #print(line)
            numofPointsVel2 = int(line[0: line.find("points")].strip())
            dtVel2 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
            unitsVel2 = line[line.find(", in") + 4: line.find(".  (")].strip()
        numofLines = lines(numofPointsVel2)
        vel2 = readchunk(f,numofLines)

        for line in islice(f, 0,1):
            #print(line)
            numofPointsDispl2 = int(line[0: line.find("points")].strip())
            dtDispl2 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
            unitsDispl2 = line[line.find(", in") + 4: line.find(".   ")].strip()
        numofLines = lines(numofPointsDispl2)
        displ2 = readchunk(f,numofLines)  

        for line in islice(f,1,24):
            i=0
        for line in islice(f, 0, 1):
            nameCh3=line[26:].strip()
        for line in islice(f, 0, 1):
            nameCh3=nameCh3 + line.strip()
            #print(nameCh3)
        for line in islice(f,1,20):
            i=0    

        for line in islice(f, 0, 1):
            #print(line)
            numofPointsAccel3 = int(line[0: line.find("points")].strip())
            dtAccel3 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
            unitsAccel3 = line[line.find(", in") + 4: line.find(". (")].strip()
        numofLines = lines(numofPointsAccel3)
        accel3 = readchunk(f,numofLines)

        for line in islice(f, 0,1):
            #print(line)
            numofPointsVel3 = int(line[0: line.find("points")].strip())
            dtVel3 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
            unitsVel3 = line[line.find(", in") + 4: line.find(".  (")].strip()
        numofLines = lines(numofPointsVel3)
        vel3 = readchunk(f,numofLines)

        for line in islice(f, 0,1):
            #print(line)
            numofPointsDispl3 = int(line[0: line.find("points")].strip())
            dtDispl3 = float(line[line.find("at ") + 3: line.find(" sec")].strip())
            unitsDispl3 = line[line.find(", in") + 4: line.find(".   ")].strip()
        numofLines = lines(numofPointsDispl3)
        displ3 = readchunk(f,numofLines)  

        T2 = np.arange(0.0,numofPointsAccel2*dtAccel2, dtAccel2)
        T3 = np.arange(0.0,numofPointsAccel1*dtAccel3, dtAccel3)
        scale = scaleValue(unitsAccel2) 
        scaledAccel2 = [value*scale for value in accel2]
        scale = scaleValue(unitsAccel3)  
        scaledAccel3 = [value*scale for value in accel3]

readFile()
win= Tk()
style = ttk.Style()
style.configure('W.TButton', font =
               ('calibri', 10, 'bold'))

style.map('TButton', foreground = [('active', '!disabled', 'red')],
                     background = [('active', 'black')])

rr=0
if EOF==1:
    win.geometry("480x430")
else:
    win.geometry("480x1050")
win.title("Read Cosmos V2 Files")

win.menubar = Menu()
win.menubar.add_command(label="Quit", command=lambda:onclick2())
win.menubar.add_command(label="About", command=lambda:myabout())  
win.config(menu = win.menubar) 

container = ttk.Frame(win)
container.pack()
canvas = Canvas(container)

canvas.grid(row=0,column=0,sticky='news')
canvas.grid_columnconfigure(1, minsize=100)
pb = ttk.Progressbar(canvas,orient='horizontal', mode='determinate', length=400)
pb.value = IntVar()


canvas.Label= Label(canvas, text="Read in memory :\n"+ filenames[0]).grid(row=rr,column=0, columnspan=2, sticky="w"); rr+=1
if EOF==1:
    canvas.Label= Label(canvas,text=nameCh1).grid(row=rr,column=0, columnspan=2, sticky="w"); rr+=1
else:
    canvas.Label= Label(canvas,text="Chan read="+nameCh1).grid(row=rr,column=0, columnspan=2, sticky="w"); rr+=1
    canvas.Label= Label(canvas,text="Chan read="+nameCh2).grid(row=rr,column=0, columnspan=2, sticky="w"); rr+=1
    canvas.Label= Label(canvas,text="Chan read="+nameCh3).grid(row=rr,column=0, columnspan=2, sticky="w"); rr+=1
canvas.Label= Label(canvas,text="Latitude = "+str(round(latitude,3))+", Longitude = "+str(round(longitude,3))).grid(row=rr,column=0, columnspan=2, sticky="w"); rr+=1
ttk.Button(canvas, text="Show location on a map", command=openMap).grid(row=rr,column=0,padx=7,columnspan = 1)
ttk.Button(canvas, text="Show Station", command=openStation).grid(row=rr,column=1,padx=7,columnspan = 1); rr+=1

canvas.LabelLx= Label(canvas, text="Start Time (secs) = ", justify="right").grid(row=rr,column=0,sticky="e")
canvas.entry_Lowxlim  = Entry(canvas)
canvas.entry_Lowxlim.insert(0,str(startlimAccel()))
canvas.entry_Lowxlim.grid(row=rr,column=1,sticky="ew"); rr+=1
canvas.LabelLx= Label(canvas, text="End Time (secs) = ", justify="right").grid(row=rr,column=0,sticky="e")
canvas.entry_Highxlim  = Entry(canvas)
canvas.entry_Highxlim.grid(row=rr,column=1,pady=7,sticky="ew"); rr+=1
#canvas.entry_Highxlim.insert(0,str(numofPointsAccel1*dtAccel1))
canvas.entry_Highxlim.insert(0,str(endlimAccel()))

ttk.Separator(canvas, orient='horizontal').grid(row=rr, column=0, columnspan=2, pady =10, sticky="ew"); rr+=1
canvas.SaveCh1 = IntVar(); canvas.SaveCh2 = IntVar(); canvas.SaveCh3 = IntVar(); canvas.createRS =IntVar(); canvas.createRS2 =IntVar();canvas.plotVel = IntVar();canvas.plotDisp = IntVar();canvas.plotOrbit = IntVar();canvas.plotFFT = IntVar()
canvas.includeASCE = IntVar();canvas.createTrip=IntVar();canvas.createAz3D=IntVar(); canvas.arias=IntVar()
if EOF==1:
    ttk.Checkbutton(canvas, text="Plot Velocity?", variable=canvas.plotVel).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="Plot Displacement?", variable=canvas.plotDisp).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="Plot FFT?", variable=canvas.plotFFT).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Separator(canvas, orient='horizontal').grid(row=rr, column=0, columnspan=2, pady =10, sticky="ew"); rr+=1
    ttk.Button(canvas, text="Plot", command=on_click).grid(row=rr,column=1,padx=7,columnspan = 1); rr+=1
else:
    canvas.label=Label(canvas, text="Parameters for Response Spectra", justify="left").grid(row=rr,column=0,sticky="w");rr+=1
    canvas.LabelDampin = Label(canvas, text="Damping Ratio = ", justify="right").grid(row=rr,column=0,sticky="e")
    canvas.entry_Damping  = Entry(canvas)
    canvas.entry_Damping.insert(0,str(0.05))
    canvas.entry_Damping.grid(row=rr,column=1,sticky="ew"); rr+=1

    canvas.LabelendPeriod = Label(canvas, text="End Period for Response Spectra (secs) = ", justify="right").grid(row=rr,column=0,sticky="e")
    canvas.entry_endPeriod  = Entry(canvas)
    canvas.entry_endPeriod.insert(0,str(6.0))
    canvas.entry_endPeriod.grid(row=rr,column=1,pady =7, sticky="ew"); rr+=1

    ttk.Separator(canvas, orient='horizontal').grid(row=rr, column=0, columnspan=2, pady =10, sticky="ew"); rr+=1

    canvas.label=Label(canvas, text="Parameters for ASCE7-22 Spectra", justify="left").grid(row=rr,column=0,sticky="w");rr+=1
    canvas.LabelLat = Label(canvas, text="Latitude = ", justify="right").grid(row=rr,column=0,sticky="e")
    canvas.entry_Lat  = Entry(canvas)
    canvas.entry_Lat.insert(0,str(latitude))
    canvas.entry_Lat.grid(row=rr,column=1,sticky="ew"); rr+=1

    canvas.LabelLong = Label(canvas, text="Longitude = ", justify="right").grid(row=rr,column=0,sticky="e")
    canvas.entry_Long = Entry(canvas)
    canvas.entry_Long.insert(0,str(longitude))
    canvas.entry_Long.grid(row=rr,column=1,pady =7, sticky="ew"); rr+=1

    SiteClassList=["A","B","BC","C","CD","D","DE","E","F"]
    canvas.SelectedSiteClass=StringVar()
    canvas.SelectedSiteClass.set("CD")
    canvas.label_SiteClass = Label(canvas,text="Site Class").grid(row=rr,column=0,sticky="e")
    canvas.list_SiteClass = OptionMenu(canvas, canvas.SelectedSiteClass,*SiteClassList)
    canvas.list_SiteClass.grid(row=rr,column=1,sticky="ew"); rr+=1

    RiskCategoryList=["I","II","III","IV"]
    canvas.SelectedRiskCategory =StringVar()
    canvas.SelectedRiskCategory.set("IV")
    canvas.label_RiskCategory = Label(canvas,text="Risk Category").grid(row=rr,column=0,sticky="e")
    canvas.list_RiskCategory = OptionMenu(canvas,canvas.SelectedRiskCategory,*RiskCategoryList)
    canvas.list_RiskCategory.grid(row=rr,column=1,sticky="ew"); rr+=1

    ttk.Separator(canvas, orient='horizontal').grid(row=rr, column=0, columnspan=2, pady =10, sticky="ew"); rr+=1


    ttk.Checkbutton(canvas, text="Create Response Spectra?", variable=canvas.createRS2).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="Create PSA vs Disp Spectra?", variable=canvas.createRS).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="Plot Arias Intensity?", variable=canvas.arias).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="Plot Velocity?", variable=canvas.plotVel).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="Plot Displacement?", variable=canvas.plotDisp).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="3D Orbit Plot?", variable=canvas.plotOrbit).grid(row=rr,column=0, sticky="w"); rr+=1

    ttk.Checkbutton(canvas, text="Include ASCE722 Spectra? (Only for Rotated Plots)", variable=canvas.includeASCE).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="Create Tripartite Response Spectra (only for Rotated Plots)?", variable=canvas.createTrip).grid(row=rr,column=0, sticky="w"); rr+=1

    ttk.Separator(canvas, orient='horizontal').grid(row=rr, column=0, columnspan=2, pady =10, sticky="ew"); rr+=1
    pb.grid(row=rr,column=0,pady=3,columnspan = 2); rr+=1
    ttk.Separator(canvas, orient='horizontal').grid(row=rr, column=0, columnspan=2, pady =10, sticky="ew"); rr+=1
    canvas.Labelplot= Label(canvas, text= "Click the button to plot (Will take time to complete with any Spectra option enabled)").grid(row=rr,column=0, sticky="w",columnspan = 2); rr+=1

    ttk.Button(canvas, text="Plot", command=on_click).grid(row=rr,column=0,padx=7,columnspan = 1)
    ttk.Button(canvas, text="Plot Rotated Max Resp", command=on_clickRotDisp).grid(row=rr,column=1,padx=7,columnspan = 1); rr+=1
    ttk.Button(canvas, text="Plot Rotated Max Velocity", command=on_clickRotVel).grid(row=rr,column=0,padx=7,columnspan = 1)
    ttk.Button(canvas, text="Plot Rotated Max Accel", command=on_clickRot).grid(row=rr,column=1,padx=7,columnspan = 1); rr+=1
    canvas.LabelAngle = Label(canvas, text="Rotate by Angle = ", justify="right").grid(row=rr,column=0,sticky="e")
    canvas.entry_Angle  = Entry(canvas); 
    canvas.entry_Angle.insert(0,str("45"))
    canvas.entry_Angle.grid(row=rr,column=1,pady =7, sticky="ew"); rr+=1
    ttk.Button(canvas, text="Plot Rotated Angle", command=on_clickRotAngle).grid(row=rr,column=1,padx=7,columnspan = 1); rr+=1
    ttk.Separator(canvas, orient='horizontal').grid(row=rr, column=0, columnspan=2, pady =10, sticky="ew"); rr+=1
    ttk.Checkbutton(canvas, text="Plot Azimuth angle vs Spectra, 3D plot?", variable=canvas.createAz3D).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Button(canvas, text="Plot RotD50 Spectra", command=on_clickRotD50).grid(row=rr,column=0,pady=3,columnspan = 2); rr+=1




ttk.Separator(canvas, orient='horizontal').grid(row=rr, column=0, columnspan=2, pady =10, sticky="ew"); rr+=1

canvas.Labelsave= Label(canvas, text= "Write data files (space separated time-acceleration values)").grid(row=rr,column=0); rr+=1
if EOF==1:

    ttk.Checkbutton(canvas, text="Save Channel 1?", variable=canvas.SaveCh1).grid(row=rr,column=0, sticky="w"); rr+=1
else:
    ttk.Checkbutton(canvas, text="Save Channel 1?", variable=canvas.SaveCh1).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="Save Channel 2?", variable=canvas.SaveCh2).grid(row=rr,column=0, sticky="w"); rr+=1
    ttk.Checkbutton(canvas, text="Save Channel 3?", variable=canvas.SaveCh3).grid(row=rr,column=0, sticky="w"); rr+=1


ttk.Button(canvas, text="Save File", command=saveFile).grid(row=rr,column=0,columnspan = 1)
ttk.Button(canvas, text="Quit", style = 'W.TButton', command=lambda:onclick2()).grid(row=rr,column=1,columnspan = 1); rr+=1


win.mainloop()