import numpy as np
import sys

class CriticalParams:
	def __init__(self,Tc=215,Pc=338.,Rhoc=0.05939):
		self.Tc = Tc
		self.Pc = Pc
		self.Rhoc = Rhoc
		self.R = 8.3144621
		self.MH2O = 0.018015268
		self.RhocRTc = self.Rhoc*self.R*self.Tc

	def reduceT(self,T):
		return (T-self.Tc)/self.Tc

	def reduceP(self,P):
		return (P-self.Pc)/self.RhocRTc

	def reduceTandP(self,T,P):
		return self.reduceT(T),self.reduceP(P)

	def invertT(self,That):
		return self.Tc*That + self.Tc

	def invertP(self,Phat):
		return self.RhocRTc*Phat + self.Pc

	def invertTandP(self,That,Phat):
		return self.invertT(That), self.invertP(Phat)

	def __repr__(self):
		return "<Critical Params:>\nT is %f K\nP is %f MPa\nDensity is %f kg/m^3" \
		%(self.Tc, self.Pc,self.Rhoc*self.MH2O*1e6)

	def __str__(self):
		print("Now printing...")
		return self.__repr__()


class DataPt:
	def __init__(self, T, P, Rho, CriticalPt=CriticalParams()):
		self.T = T
		self.P = P
		self.Rho = Rho
		self.CriticalParams = CriticalPt

	def set_critical_pt(self,critical_pt):
		self.CriticalParams = critical_pt

	def reducedUnits(self):
		return self.CriticalParams.reduceTandP(self.T, self.P)

	def __repr__(self):
		return "<Real Units:>\nT is %f K\nP is %f MPa\nDensity is %f kg/m^3" \
		%(self.T, self.P,self.Rho)


class RawData:
	def __init__(self, fn="compiledEOSAll"):
		x_T = np.loadtxt(fn,unpack=True,usecols=[0],skiprows=1)
		x_P = np.loadtxt(fn,unpack=True,usecols=[1],skiprows=1)
		x_Rho = np.loadtxt(fn,unpack=True,usecols=[2],skiprows=1)
		self.AllData = np.stack([x_T, x_P, x_Rho],axis=1)

	def getRawData(self):
		return self.AllData

	def trimTemperatureBounds(self,Tl=sys.float_info.min,Tu=sys.float_info.max):
		self.AllData = self.AllData[np.logical_not(self.AllData[:,0] > Tu)]
		self.AllData = self.AllData[np.logical_not(self.AllData[:,0] < Tl)]

	def trimPressureBounds(self,Pl=sys.float_info.min,Pu=sys.float_info.max):
		self.AllData = self.AllData[np.logical_not(self.AllData[:,1] > Pu)]
		self.AllData = self.AllData[np.logical_not(self.AllData[:,1] < Pl)]

	def trimPT(self,funcdict):
		return 0.

	def getTemperatures(self,copy=False):
		i=0
		return self.AllData[:,i].copy() if copy else self.AllData[:,i]

	def getPressures(self,copy=False):
		i=1
		return self.AllData[:,i].copy() if copy else self.AllData[:,i]

	def getDensities(self,copy=False):
		i=2
		return self.AllData[:,i].copy() if copy else self.AllData[:,i]

	def plotAll(self):
		return 0.



if __name__ == "__main__":
	myrawdata = RawData()
	cpt = CriticalParams(180,330,0.00545)
	print("Your move...")
