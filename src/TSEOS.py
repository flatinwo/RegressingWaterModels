import numpy as np
import pandas as pd
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
	def __init__(self, fn="./Data/compiledEOSAll"):
		x_T = np.loadtxt(fn,unpack=True,usecols=[0],skiprows=1)
		x_P = np.loadtxt(fn,unpack=True,usecols=[1],skiprows=1)
		x_Rho = np.loadtxt(fn,unpack=True,usecols=[2],skiprows=1)
		self.AllData = np.stack([x_T, x_P, x_Rho],axis=1)
		self.IsochoreDict = {}
		self.buildDict = True

	def getRawData(self):
		return self.AllData

	def trimTemperatureBounds(self,Tl=sys.float_info.min,Tu=sys.float_info.max):
		self.AllData = self.AllData[np.logical_not(self.AllData[:,0] > Tu)]
		self.AllData = self.AllData[np.logical_not(self.AllData[:,0] < Tl)]
		self.buildDict = True

	def trimPressureBounds(self,Pl=sys.float_info.min,Pu=sys.float_info.max):
		self.AllData = self.AllData[np.logical_not(self.AllData[:,1] > Pu)]
		self.AllData = self.AllData[np.logical_not(self.AllData[:,1] < Pl)]
		self.buildDict = True

	def trimPT(self,Pl=0,Tl=0):
		self.AllData = self.AllData[np.logical_and(self.AllData[:,0] > Tl, self.AllData[:,1] > Pl)]
		self.buildDict = True

	def trimAlongIsochore(self,Density=860,Tl=0,complement=False,densitymin=860):
		assert(Density>=densitymin)
		if complement:
			self.AllData = self.AllData[np.logical_and(self.AllData[:,0] > Tl, self.AllData[:,2] == Density)] # an associated array
		#be better for this
		else:
			self.AllData = self.AllData[np.logical_not(np.logical_and(self.AllData[:,0] < Tl, self.AllData[:,2] == Density))]
		self.buildDict = True

	def getTemperatures(self,copy=False):
		i=0
		return self.AllData[:,i].copy() if copy else self.AllData[:,i]

	def getPressures(self,copy=False):
		i=1
		return self.AllData[:,i].copy() if copy else self.AllData[:,i]

	def getDensities(self,copy=False):
		i=2
		return self.AllData[:,i].copy() if copy else self.AllData[:,i]

	def updateIsochoreDict(self):
		self.IsochoreDict = {}
		for i in range(len(self.AllData)):
			if self.AllData[i,2] not in self.IsochoreDict: self.IsochoreDict[self.AllData[i,2]] = [];
			self.IsochoreDict[self.AllData[i,2]].append([ self.AllData[i,0],self.AllData[i,1] ])
		self.buildDict = False

	def writeData(self,filename=None):
		if filename is None: print("Please specify filename")
		else:
			df = pd.DataFrame(self.AllData,columns=["Temperature(K)",
												"Pressure (MPa)",
												"Density (kg/m^3)"])
			df.to_csv(filename,header=True,index=False,sep="\t")




	def plotAll(self,wait=False,skip=[],style="seaborn-colorblind",savefig=False,ncol=4,figsize=(10,10),
		filename="./Data/RawEOS.pdf"):
		try:
			if self.buildDict: self.updateIsochoreDict()
			import matplotlib.pyplot as plt
			from matplotlib import rc
			from matplotlib.ticker import AutoMinorLocator
			rc('font',family="Times New Roman")
			rc('text',usetex=True)
			rc('ytick',labelsize=24)
			rc('xtick',labelsize=24)
			rc('font',size=30)
			rc('ytick.major',size=5)
			#import seaborn as sns #pick bigger color cycle fancier method
			hsv = plt.get_cmap('viridis')
			colors = hsv(np.linspace(0,1.0,len(self.IsochoreDict)))
			plt.figure(figsize=figsize)
			minorLocator = AutoMinorLocator(2)

			with plt.style.context(style):
				for isochore,color in zip(sorted(self.IsochoreDict),colors):
					T = np.array(self.IsochoreDict[isochore])[:,0]
					P = np.array(self.IsochoreDict[isochore])[:,1]
					if isochore not in skip: 
						plt.plot(T,P,'o',label=str(isochore),color=color, markersize=12)
					else:
						plt.plot(T,P,'-',label=(isochore), color=color, linewidth=2)

				plt.legend(numpoints=1,loc='lower right',ncol=ncol,fontsize=10,markerscale=1.,frameon=False)
				plt.xlabel('Temperature (K)')
				plt.axes().xaxis.set_minor_locator(minorLocator)
				plt.axes().yaxis.set_minor_locator(minorLocator)
				plt.ylabel('Pressure (MPa)')
				plt.tick_params(which='both', width=1)
				plt.tick_params(which='major', length=7)
				plt.tick_params(which='minor', length=4)
				if not wait: plt.show() if not savefig else plt.savefig(filename,bbox_inches='tight',dpi=600)

		except:
			print("Plotting error")



if __name__ == "__main__":
	myrawdata = RawData(fn="./Data/pruned025compiledEOSAll")
	myrawdata.plotAll(savefig=False)
	cpt = CriticalParams(180,330,0.00545)
	print("Your move...")
