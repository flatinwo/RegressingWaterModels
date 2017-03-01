import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc  
from matplotlib.widgets import Slider, Button, RadioButtons
from load_data import *

Data1 = load_data(filename="../Data/a2_003_3.pkl")
Data2 = load_data(filename="../Data/a2_01_3.pkl")

rc('ytick',labelsize=20)
rc('xtick',labelsize=20)

x01,y01 = Data1["TMD"]
x02,y02 = Data2["TMD"]

x11,y11 = Data1["LLT"]
x12,y12 = Data2["LLT"]

x21,y21 = Data1["LLW"]
x22,y22 = Data2["LLW"]

x31,y31 = Data1["LVS"]
x32,y32 = Data2["LVS"]

fig, ax = plt.subplots(figsize=(10,10))
plt.subplots_adjust(left=0.25,bottom=0.30)

plt.plot(x01,y01,'-',color="red",lw=2, label="TMD")
plt.plot(x11,y11,'-',color="blue",lw=2,label="LLT")
plt.plot(x21,y21,'-',color="green",lw=2,label="LLW")
plt.plot(x31,y31,'-',color="black",lw=2,label="LVS")

l1,=plt.plot(x02,y02,'-.',color="red",lw=2)
l2,=plt.plot(x12,y12,'-.',color="blue",lw=2)
l3,=plt.plot(x22,y22,'-.',color="green",lw=2)
l4,=plt.plot(x32,y32,'-.',color="black",lw=2)

plt.xlabel('Temperature',fontsize=24)
plt.ylabel('Pressure',fontsize=24)
plt.legend(frameon=False,loc='upper right')

ax.set_ylim(-20,40)
ax.set_xlim(0.,1.0)

axcolor = 'lightgoldenrodyellow'
axPscale = plt.axes([0.25,0.15,0.65,0.03],axisbg=axcolor)
axTscale = plt.axes([0.25,0.05,0.65,0.03],axisbg=axcolor)

Pinit = 1.0
Tinit = 1.0


axPfreq = Slider(axPscale,'ps',-5,5.,valinit=Pinit)
axTfreq = Slider(axTscale,'ts', -2.,2.,valinit=Tinit)

def update(val):
	Ps, Ts = axPfreq.val, axTfreq.val
	l1.set_ydata(Ps*y02)
	l2.set_ydata(Ps*y12)
	l3.set_ydata(Ps*y22)

	l1.set_xdata(Ts*x02)
	l2.set_xdata(Ts*x12)
	l3.set_xdata(Ts*x22)

axPfreq.on_changed(update)
axTfreq.on_changed(update)

resetax = plt.axes([0.05,0.1,0.1,0.04])
button  = Button(resetax, 'Reset', color=axcolor)

def reset(event):
	axTfreq.reset()
	axPfreq.reset()

button.on_clicked(reset)


plt.show()







