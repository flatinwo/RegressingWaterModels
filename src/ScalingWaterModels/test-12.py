import GuggenWater as gw
import numpy as np

T,P,Feature,Src = np.loadtxt('Features.dat',dtype={'names':('T','P','F','S'),'formats':(np.float,np.float,'S9','S9')},
 usecols=(0,1,2,3),unpack=True,skiprows=1)

Property = np.array([str(s,'utf-8') for s in Feature])
Source = np.array([str(s,'utf-8') for s in Src])


mycb = "classicAndWidom"

model = gw.GuggenheimWater()
for t,p,ppty,src in zip(T,P,Property, Source): model.addPoint(t,p,ppty,src,callback=mycb)
model.plot(savefig=True,skip=['Widom'],style='ggplot',wait=True)
model.writetoFile('All.dat',callback=mycb)

st2 = gw.ST2I_mean_field()
#st2.spinodal[0] -= 250.
model2 = gw.GuggenheimWater(st2)
T,P,Feature,Src = np.loadtxt('FeaturesST2.dat',dtype={'names':('T','P','F','S'),'formats':(np.float,np.float,'S9','S9')},
 usecols=(0,1,2,3),unpack=True,skiprows=1)

Property = np.array([str(s,'utf-8') for s in Feature])
Source = np.array([str(s,'utf-8') for s in Src])

offset = 0.7044#/(gw.ST2I_mean_field().getSpinodalCriticalPointOffset() / gw.TIP4P2005().getSpinodalCriticalPointOffset())

for t,p,ppty,src in zip(T,P,Property,Source): model2.addPoint(t,p,ppty,src,offset,callback=mycb)
model2.writetoFile("All.dat",offsetCriticalPt=offset,isnew=False,callback=mycb)
#model2.plot(savefig=True,skip=['Spinodal'],style='ggplot',wait=False,marker='*')
