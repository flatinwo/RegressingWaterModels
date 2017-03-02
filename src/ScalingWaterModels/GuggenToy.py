from GuggenWater import TIP4P2005,GuggenheimWater
from TData import CriticalParams
from load_data import *
import numpy as np
from numpy.polynomial.polynomial import polyval2d,polyval
from numpy import poly1d

class ToyModel(TIP4P2005):
	def __init__(self,Tc,Pc,a2=-0.3):
		self.CriticalParams = CriticalParams(Tc,Pc,1.)
		L0, a, b, d, f, e = 0.5, a2, 0., -3.,0., 0.
		self.widom = np.array([[L0, a*L0, b*L0],[d*L0, f*L0, 0.], [e*L0, 0., 0.]])
		self.model_name = "LatticeGasRx"
		self.Tc = Tc
		self.Pc = Pc

	def WidomLine(self, T, P):
		return polyval2d(T,P,self.widom)

	def SpinodalPressure(self,T):
		p = np.poly1d([1., -1., T/4])
		rho = np.max(p.r)
		#print((T,rho,self.widom[0][1]))
		np.seterr(all='print')
		try:
			P = (1./(0.5+np.log(0.5)))*(T*np.log(1.-rho) + 2*rho*rho)
		except:
			print((T,rho,self.widom[0][1]))
			P = -10.5
		return P


	def WidomTemperature(self, P):
		Troots = ((np.poly1d(polyval(P,np.transpose(self.widom))[::-1])).r) #roots for temperature

		#check for good root
		T = np.extract(np.isreal(Troots),Troots)

		if len(T) == 0:
			raise ValueError("No real roots found for Pressure as %f" % P)
		elif len(T) == 1:
			return float(T)
		else:
			newT = np.extract(np.array([x >= 0. and x <= 1. for x in T]),T)
			if len(newT) > 1:
				raise ValueError("Still too many real roots found %f" % newT)
			else:
				return float(newT)


class GuggenModel(GuggenheimWater):
	def __init__(self,model):
		assert isinstance(model,ToyModel)
		self.model = model
		self.havefigure = False
		self.Trescale = []
		self.Prescale = []
		self.Property = []
		self.Source = []
		self.PropertyDict = {}
		self.Rescalings = {}
		self.l = 1.


if __name__ == "__main__":
	ty = ToyModel(1.,1.)
	gw = GuggenModel(ty)
	Data = load_data(filename="../Data/a2_03_3.pkl")

	for ppty in Data:
		T,P = Data[ppty]
		try: 
			for t,p in zip(T,P): gw.addPoint(t,p,ppty,Source="Unknown",callback="identity")
		except:
			gw.addPoint(T,P,ppty,Source="Unknown",callback="identity")
			

	gw.plot(marker='-',
		xlabel="Temperature",
		ylabel="Pressure",
		xlim=(0.0,1.0),
		ylim=(-20,60),
		style="seaborn-colorblind")


	ty = ToyModel(1.,1.)
	gw = GuggenModel(ty)
	Data = load_data(filename="../Data/a2_03_3.pkl")

	for ppty in Data:
		T,P = Data[ppty]
		try: 
			for t,p in zip(T,P): gw.addPoint(t,p,ppty,Source="Unknown",callback="delineateSpinodalandWidomLine")
		except:
			pass
			#gw.addPoint(T,P,ppty,Source="Unknown",callback="delineateSpinodalandWidomLine")
	
	ppty="LLCP"
	T,P = Data[ppty]
	gw.addPoint(T,P,ppty,Source="Unknown",callback="delineateSpinodalandWidomLine")
			

	gw.plot(marker='-',
		xlabel="Temperature",
		ylabel="Pressure",
		skip=["LLCP"],
		xlim=(-1.0,1.0),
		ylim=(-10,20),
		style="seaborn-colorblind", wait=False)

	ty = ToyModel(1.,1.,a2=-0.1)
	gw = GuggenModel(ty)
	Data = load_data(filename="../Data/a2_01_3.pkl")

	for ppty in Data:
		T,P = Data[ppty]
		try: 
			for t,p in zip(T,P): gw.addPoint(t,p,ppty,Source="Unknown",callback="delineateSpinodalandWidomLine")
		except:
			pass
	
	ppty="LLCP"
	T,P = Data[ppty]
	gw.addPoint(T,P,ppty,Source="Unknown",callback="delineateSpinodalandWidomLine")
			

	gw.plot(marker='-',
		skip=["LLCP"],
		xlabel="Temperature",
		ylabel="Pressure",
		xlim=(-1.0,1.0),
		ylim=(-10,20),
		style="seaborn-colorblind", wait=False)



