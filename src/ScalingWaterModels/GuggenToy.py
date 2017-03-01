from GuggenWater import TIP4P2005,GuggenheimWater
from TData import CriticalParams
from load_data import *
import numpy as np

class ToyModel(TIP4P2005):
	def __init__(self,Tc,Pc):
		self.CriticalParams = CriticalParams(Tc,Pc,1.)

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
	Data = load_data(filename="../Data/a2_003_3.pkl")

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