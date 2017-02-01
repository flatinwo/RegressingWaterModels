import GuggenWater as gw
import numpy as np

T,P,Feature,Src = np.loadtxt('Features.dat',dtype={'names':('T','P','F','S'),'formats':(np.float,np.float,'S9','S9')},
 usecols=(0,1,2,3),unpack=True,skiprows=1)

Property = np.array([str(s,'utf-8') for s in Feature])
Source = np.array([str(s,'utf-8') for s in Src])


mycb = "delineateSpinodalandWidomLine"
mycb2 = "delineateCriticalPoint"

K = 480.
model = gw.GuggenheimWater()
for t,p,ppty,src in zip(T,P,Property, Source): model.addPoint(t,p,ppty,src,callback=mycb)
# get critical transition data
PPst = 444. #model.PropertyDict['Tmax'][221][1] # where two curves cross each other
Pt =   model.PropertyDict['Tmax'][221][4] # where two curves cross each other
Tt =   model.PropertyDict['Tmax'][221][3]
model.clearAll()
model.findScalingConstantPerProperty(Pt,PPst,Tt)
model.K = K
model.l = 1.#250.916
for t,p,ppty,src in zip(T,P,Property, Source): model.addPoint(t,p,ppty,src,callback=mycb2)
model.plot(savefig=True,skip=['Widom'],style='ggplot',wait=True)
model.writetoFile('All.dat',callback=mycb2)



st2 = gw.ST2I_mean_field()
model2 = gw.GuggenheimWater(st2)
T,P,Feature,Src = np.loadtxt('FeaturesST2.dat',dtype={'names':('T','P','F','S'),'formats':(np.float,np.float,'S9','S9')},
 usecols=(0,1,2,3),unpack=True,skiprows=1)

Property = np.array([str(s,'utf-8') for s in Feature])
Source = np.array([str(s,'utf-8') for s in Src])

offset = 0.7044#/(gw.ST2I_mean_field().getSpinodalCriticalPointOffset() / gw.TIP4P2005().getSpinodalCriticalPointOffset())

for t,p,ppty,src in zip(T,P,Property,Source): model2.addPoint(t,p,ppty,src,offset,callback=mycb)
PPst = 309 #model2.PropertyDict['Tmax'][21][1] #where the two curves interesect
Pt =  model2.PropertyDict['Tmax'][21][4] #where the two curves cross
Tt = model2.PropertyDict['Tmax'][21][3]
model2.clearAll()
model2.findScalingConstantPerProperty(Pt,PPst,Tt)
model2.K = K
model2.l = 253.5/182.#308.3
for t,p,ppty,src in zip(T,P,Property,Source): model2.addPoint(t,p,ppty,src,offset,callback=mycb2)
model2.writetoFile("All.dat",offsetCriticalPt=offset,isnew=False,callback=mycb2)
#model2.plot(savefig=True,skip=['Spinodal'],style='ggplot',wait=False,marker='*')

