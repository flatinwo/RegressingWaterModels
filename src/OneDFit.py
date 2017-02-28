import numpy as np
import TSEOS as ts
import pandas as pd 

class OneDFit:
	def __init__(self,name, RawData=None):
		self.__name = name
		self.fitParams = {}
		self.checkData = False
		if RawData is not None:
			self.Rawdata = RawData
			self.checkData = True

	def setName(self,name):
		self._name = name

	def getParams(self):
		if not self.checkData:
			print("No data, so no parameters")
			return
		else:
			return self.fitParams

class SpinodalFit(OneDFit):
	def __init__(self,fn='./Data/compiledEOSAll',Tl=215):
		OneDFit.__init__(self,name='spinodal',RawData=ts.RawData(fn))
		self.datapt = ts.DataPt(293,101.325,1000.)
		self.Rawdata.trimTemperatureBounds(Tl=Tl)
		self.fitParams['P(T) params'] = 0.
		self.fitParams['A(T) params'] = 0.
		self.setUpLog = False
		self.order0=2 # this is a polynomial for P vs V, change at your own risk
		self.order1=2 # analysis of anova shows that quadratic fit here is stat sig.
		self.countsl=5 # optimized number of observations to select

	def setFitInfo(self):
		return 0.

	def setCriticalPt(self,Tc,Pc,Rhoc):
		self.datapt.set_critical_pt(CriticalParams(Tc,Pc,Rhoc))

	def fit(self,Pmax=0.,writeData=False):
		SpinodalData=self.Rawdata.AllData.copy()
		
		SpinodalData=SpinodalData[np.less_equal(SpinodalData[:,1],Pmax)]
		T_unique = np.unique(SpinodalData[:,0])
		
		CalculatedSpinodal = np.zeros_like(SpinodalData[0:len(T_unique),:])
		DPDV = np.zeros_like(CalculatedSpinodal)

		i=0
		for Tl in T_unique:
			Sp=SpinodalData[np.equal(SpinodalData[:,0],Tl)]
			Sp=Sp[Sp[:,1].argsort()]
			pl, rhol = Sp[0:self.countsl,1],Sp[0:self.countsl,2]
			vl=1./rhol
			p = np.poly1d(np.polyfit(vl,pl,self.order0))
			vls = np.real(p.deriv().r)[0]
			rhols = 1./vls
			CalculatedSpinodal[i,:] = Tl, p(vls), rhols
			DPDV[i,:] = Tl, p.deriv()(vls), p.deriv().deriv()(vls)
			i = i+1


		A = np.sqrt(8.)/(3.*np.sqrt(DPDV[:,2]))
		Tshat, Pshat = self.datapt.CriticalParams.reduceTandP(CalculatedSpinodal[:,0],CalculatedSpinodal[:,1])

		# fit A vs T
		try: 
			p = np.polyfit(Tshat,A,self.order1,full=True)
		except: 
			print("Warning:\nError in A vs T fit")
			print(np.array([Tshat,A]))
		self.fitParams['A(T) params'] = np.poly1d(p[0])
		res1 = float(p[1])

		#fit Ps vs Ts
		polyorder = 2 #switching to quadratic fit after ANOVA analysis with p-value < 0.001 
		p = np.polyfit(Tshat,CalculatedSpinodal[:,1],polyorder,full=True) 
		self.fitParams['P(T) params']=np.poly1d(p[0]) 
		res2 = float(p[1])

		print(("Residuals for A(T) params and P(T) params are %e and %e respectively")%(res1,res2))

		if writeData:
			PT = np.transpose(np.array([Tshat,CalculatedSpinodal[:,1]]))
			df = pd.DataFrame(PT,columns=['reduced_temperature','pressure'])
			df.to_csv(r'./Data/PTspinodal.txt',sep='\t',index=False,header=True)
			print("Written spinodal to PTspiondal.txt")

			AT = np.transpose(np.array([Tshat,A]))
			df = pd.DataFrame(AT,columns=['reduced_temperature','Afactor'])
			df.to_csv(r'./Data/ATspinodal.txt',sep='\t',index=False,header=True)
			print("Written spinodal to ATspiondal.txt")


	def getFigure(self):
		return 0



if __name__ == "__main__":
	mysf = SpinodalFit(r'./Data/pruned025compiledEOSAll',Tl=230)
	#mysf = SpinodalFit()
	mysf.fit(writeData=True)
	print("Your move...")
