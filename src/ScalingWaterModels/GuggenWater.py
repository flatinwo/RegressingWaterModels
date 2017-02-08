import TData as td
import numpy as np
from numpy.polynomial.polynomial import polyval2d,polyval
from numpy import poly1d

class TIP4P2005:
	def __init__(self,CP=td.CriticalParams(182.,170.,0.05645)):
		assert isinstance(CP,td.CriticalParams)
		self.CriticalParams = CP
		self.spinodal = np.array([-462.,475.02,-215.306])
		L0, a, b, d, f = 1.55607, 0.154014, 0.125093, 0.0085418, 1.14576
		self.widom = np.array([[0., a*L0, d*L0],[L0, b*L0, 0.], [f*L0, 0., 0.]])
		self.model_name = "TIP4P2005"

	def getSpinodalCriticalPointOffset(self):
		return self.CriticalParams.Pc - self.spinodal[0]

	def WidomLine(self,T,P):
		That, Phat = self.CriticalParams.reduceTandP(T,P)
		return polyval2d(That,Phat,self.widom)

	def WidomLine_T(self,T,P):
		"""Derivative of Widom Line with respect to temperature"""
		That,Phat = self.CriticalParams.reduceTandP(T,P)
		return polyval2d(That,Phat,polyder(self.widom,axis=0))

#	@classmethod
	def SpinodalPressure(self,T):
		That = self.CriticalParams.reduceT(T)
		return polyval(That,self.spinodal)

#	@classmethod
	def WidomTemperature(self,P):
		Phat = self.CriticalParams.reduceP(P)
		Thatroots = ((poly1d(polyval(Phat,np.transpose(self.widom))[::-1])).r) #roots for temperature
		Troots = self.CriticalParams.invertT(Thatroots)

		#check for good root
		T = np.extract(np.isreal(Troots),Troots)

		if len(T) == 0:
			raise ValueError("No real roots found for Pressure as %f" % P)
		elif len(T) == 1:
			return float(T)
		else:
			newT = np.extract(np.array([x >= 100. and x <= 600. for x in T]),T)
			if len(newT) > 1:
				raise ValueError("Still too many real roots found %f" % newT)
			else:
				return float(newT)


class ST2I_mean_field(TIP4P2005):
	def __init__(self,CP=td.CriticalParams(253.5,160.,0.052478)):
		assert isinstance(CP,td.CriticalParams)
		self.CriticalParams = CP
		self.spinodal = np.array([-234.38-50.,628.41,-720.12])
		L0, a, b, d, f = 3.4915, 0.085811, 0., 0., 0.
		self.widom = np.array([[0., a*L0, d*L0],[L0, b*L0, 0.], [f*L0, 0., 0.]])
		self.model_name = "ST2"



class GuggenheimWater:
	def __init__(self,model=TIP4P2005()):
		self.model = model
		self.havefigure = False
		self.Trescale = []
		self.Prescale = []
		self.Property = []
		self.Source = []
		self.PropertyDict = {}
		self.Rescalings = {}
		self.l = 1.

	def rescale(self,T,P,offsetCriticalPt=0.8676):
		Tw, Ps =  self.model.WidomTemperature(P),self.model.SpinodalPressure(T)
		Tc, Pc = self.model.CriticalParams.Tc, self.model.CriticalParams.Pc
		if (P-Ps) < 0.75*(Pc-Ps):
			Tr, Pr = (T-Tw)/100.,(P-Ps)/644
		else:
			Tr, Pr = (T-Tw)/100., P/Pc
		return Tr, Pr

	def delineateSpinodal(self,T,P):
		return T, P - self.model.SpinodalPressure(T)

	def delineateWidomLine(self,T,P):
		return T-self.model.WidomTemperature(P), P

	def delineateSpinodalandWidomLine(self, T, P):
		return T-self.model.WidomTemperature(P), P - self.model.SpinodalPressure(T)

	def delineateCriticalPoint(self,T,P,l=1.):
		Ps = self.model.SpinodalPressure(T)
		PPs = P - Ps
		TTw = T - self.model.WidomTemperature(P)
		Tc = self.model.CriticalParams.Tc
		Pc = self.model.CriticalParams.Pc
		PsTc = self.model.SpinodalPressure(Tc)
		if PPs <= self.PPst : return T-Tc, P-Pc #TTw/Tc, PPs/Pc 
		else: 
			return T-Tc, P-Pc #TTw/Tc, PPs/(Pc) #

	def identity(self,T,P):
		return T,P

	def findScalingConstantPerProperty(self, Pt, PPst,Tt):
		Pc,Tc = self.model.CriticalParams.Pc, self.model.CriticalParams.Tc
		PsTc = self.model.SpinodalPressure(Tc)
		PsTt = self.model.SpinodalPressure(Tt)
		print("My (Pt, PPst, Tt) is %f, %f, %f"%(Pt,PPst,Tt))
		self.K = (Pc - PsTc)/((Pt-PsTc)*PPst)
		self.PPst = PPst
		self.l = (Pt - 0.5*Pc)/(PsTt*(1+0.5))

	def clearAll(self):
		mymodel = self.model
		self.__init__(model=mymodel)


	def addPoint(self,T,P,Property="Unknown",Source="Unknown",offsetCriticalPt=1.,callback="rescale"):
		try:
			func = getattr(self,str(callback))
		except:
			print("Callback function error, valid options are: delineateSpiondal\ndelineateWidomLine\n"\
				"delineateSpinodalandWidomLine\ndelineateCriticalPoint")
			exit(-1)
		x,y = func(T,P)
		self.Trescale.append(x)
		self.Prescale.append(y)
		self.Property.append(Property)
		self.Source.append(Source)

		if Property not in self.PropertyDict: self.PropertyDict[Property] = []
		self.PropertyDict[Property].append([x,y,Source,T,P])




	def writetoFile(self,filename="RescaledTIP4P2005.dat",isnew=True, add_critical_pt=True,offsetCriticalPt=1.,callback="rescale"):
		if add_critical_pt: self.addPoint(self.model.CriticalParams.Tc, self.model.CriticalParams.Pc,"LLCP","TSEOS",offsetCriticalPt,callback)
		if not self.PropertyDict:
			print("No entries have been added to dictionary. Drops mic!")
		else:
			model_name = self.model.model_name
			with open(filename,'a') as fn:
				if isnew: fn.write('Temperature\tPressure\tProperty\tSource\tModel\n')
				for myproperty in self.PropertyDict:
					T,P,S = np.transpose(np.array(self.PropertyDict[myproperty]))[0:3]
					for i in np.arange(len(T)):
						fn.write('%f\t%f\t%10s\t%10s\t%10s\n'%(float(T[i]),float(P[i]),myproperty,S[i],model_name))
			fn.close()



	def plot(self,savefig=False,skip=[],style='fivethirtyeight',wait=False,marker='o'):
		try:
			import matplotlib.pyplot as plt
			with plt.style.context(style):
				for myproperty in self.PropertyDict:
					T = np.array(self.PropertyDict[myproperty])[:,0]
					P = np.array(self.PropertyDict[myproperty])[:,1]
					if myproperty not in skip: 
						plt.plot(T,P,marker,label=myproperty, markersize=12)
					else:
						plt.plot(T,P,'-',label=myproperty, linewidth=2)

				plt.legend(numpoints=1)
				plt.xlabel('Spinodal Line',fontsize=20)
				plt.ylabel('Widom Line',fontsize=20)
				#plt.axes().set_aspect('equal','datalim')
				if not wait: plt.show() if not savefig else plt.savefig('RescaledPhaseDiagram.png',bbox_inches='tight',dpi=300)

		except:
			print("Plotting error")

if __name__ == "__main__":
	GW = GuggenheimWater()
	GW.addPoint(200,140,'Kt')
	GW.addPoint(230,120,'Kt')
	GW.addPoint(320,110,'Cp')
	GW.plot(savefig=True)
