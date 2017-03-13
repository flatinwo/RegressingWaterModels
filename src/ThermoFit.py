import numpy as np
import OneDFit as od
import TSEOS as ts
from scipy.optimize import brentq, minimize_scalar, leastsq
from numpy.polynomial.polynomial import polyval,polyval2d,polyder

def getXs(cxy,x,y):
	f   = polyval2d(x,y,cxy)
	fx  = polyval2d(x,y, polyder(cxy,axis=0))
	fy  = polyval2d(x,y, polyder(cxy,axis=1))
	fxx = polyval2d(x,y, polyder(polyder(cxy,axis=0),axis=0))
	fxy = polyval2d(x,y, polyder(polyder(cxy,axis=0),axis=1))
	fyy = polyval2d(x,y, polyder(polyder(cxy,axis=1),axis=1))
	return f,fx,fy,fxx,fxy,fyy


class ThermoPPty:
	def __init__(self,fn='./Data/compiledEOSAll', Tc=200, Pc=200, Rhoc=1.,Tl=230.):
		sf = od.SpinodalFit(fn,Tl) #Maybe i don't need this as a member fxn
		sf.setCriticalPt(Tc,Pc,Rhoc) #But ensures I have the same criticalpt
		sf.fit()

		self.compute = False
		self.criticalpt = ts.CriticalParams(Tc,Pc,Rhoc)
		
		self.cA = np.array(sf.fitParams['A(T) params'])[::-1]
		self.cS = np.array(sf.fitParams['P(T) params'])[::-1]
		self.Bg = np.zeros([5,5])
		self.Ws = np.zeros([5,5])
		self.Ls = np.zeros([5,5])
		self.delw = 1.

		self.fitparams = np.zeros(21)

		self.UniqueDensities = np.array(list(sf.Rawdata.IsochoreDict.keys()))
		self.IsochoreDict = {}
		for density in sf.Rawdata.IsochoreDict:
			if density not in self.IsochoreDict: self.IsochoreDict[density] = {}
			for T, P in sf.Rawdata.IsochoreDict[density]:
				self.IsochoreDict[density][T] = P
		
		self.f = lambda x, p : p[0]*x + p[1]*(x*np.log(x) + (1-x)*np.log(1-x)) + p[2]*x*(1-x)

	def getPressure(self,T=200.,Density=1000, guess=True, Pguess=0., dRho=20.):
		my_nearest = lambda array, value : array[np.abs(array-value).argmin()]
		RhoMax = my_nearest(self.UniqueDensities,Density+dRho)
		
		Tn = my_nearest(np.array(list(self.IsochoreDict[RhoMax].keys())),T)
		Pmax = self.IsochoreDict[RhoMax][Tn]
		myres = lambda P: Rho - self.getDensity(T,P,False)
		if not guess: Pmax = Pguess
		plsq = leastsq(myres,Pmax)
		return plsq[0]

	def getDensity(self,T=200.,P=100.,computed=False):
		if not computed: self.__compute(T,P)
		return self.criticalpt.Rhoc*self.criticalpt.MH2O*1e6/self.v # result in kg/m^3

	def getCompressibility(self,T=200,P=100.,computed=False):
		if not computed: self.__compute(T,P)
		return self.kt/(self.criticalpt.RhocRTc) # result in reciprocal MPa

	def getCp(self,T=200.,P=100.,computed=False):
		if not computed: self.__compute(T,P)
		return self.cp*self.criticalpt.R/self.criticalpt.MH2O #result in J kg^-1 K^-1

	def getThermalExpansivity(self,T=200,P=100.,computed=False):
		if not computed: self.__compute(T,P)
		return self.ap/self.criticalpt.Tc # result in inverse temperature

	def getCv(self,T=200.,P=100.,computed=False):
		cp = self.getCp(T,P,computed)
		ap = self.getThermalExpansivity(T,P,True)
		rho = self.getDensity(T,P,True)
		kt = self.getCompressibility(T,P,True)
		return cp - 1e6*T*ap*ap/(rho*kt) #result in J kg^-1 K^-1

	def getSpeedofSound(self,T=200.,P=100.,computed=False):
		cp = self.getCp(T,P,computed)
		ap = self.getThermalExpansivity(T,P,True)
		rho = self.getDensity(T,P,True)
		kt = self.getCompressibility(T,P,True)
		cv = cp - 1e6*T*ap*ap/(rho*kt) #result in J kg^-1 K^-1
		return np.sqrt(cp/(rho*cv*kt)) if cp/(rho*cv*kt) > 0 else print("Imaginary speed of sound") #what units

	def __str__(self,T=200.,P=100.):
		cp = self.getCp(T,P,self.compute)
		ap = self.getThermalExpansivity(T,P,True)
		rho = self.getDensity(T,P,True)
		kt = self.getCompressibility(T,P,True)
		cv = cp - 1e6*T*ap*ap/(rho*kt) #result in J kg^-1 K^-1
		cs = np.sqrt(cp/(rho*cv*kt))
		return "<At this given state with T=%f and P=%f, the response functions are the following>\
		\nDensity is %f kg/m^3\nc_P is %f J.(kgK)^-1\nc_V is %f J.(kgK)^-1\nk_T is %f MPa^-1\nalpha_p\
		is %f K^-1\nc_s is %f\n"%(T,P,rho,cp,cv,kt,ap,cs)


	def __compute(self, T=200,P=100,ind=0):
		That, Phat = self.criticalpt.reduceTandP(T,P)
		tau = 1. + That

		A = polyval(That,self.cA)
		Ps = polyval(That,self.cS)
		Pshat,S0 = self.criticalpt.reduceP(np.array([Ps,self.cS[0]]))
		PminusPhat = Phat - Pshat

		if PminusPhat < 0: PminusPhat=0.
		if S0 > 0: print("there was a numerical error with spinodal")
		
		Vs = 1.5*A*np.sqrt(PminusPhat) #spinodal volume... write for case without spinodal

		# get reaction term, that is -lnK or in this case L
		L, Lt, Lp, Ltt, Ltp, Lpp = getXs(self.Ls,That,Phat)
		W, Wt, Wp, Wtt, Wtp, Wpp = getXs(self.Ws,tau,Phat)

		params = np.array([L,W,tau,Phat])
		x = self.x = self.getEquibConc(params)
		dxdt = self.dxdt = (Lt + np.log(x/(1-x)) + Wt*(1.-2*x))/(2*W - tau/(x*(1.0-x)))
		dxdp = self.dxdp = (Lp + Wp*(1-2*x))/(2*W - tau/(x*(1-x)))

		#dGA/dP
		Bp = polyval2d(That,Phat,polyder(self.Bg,axis=1))
		w0 = 2.0*(self.Ws[0][1]/self.Ws[0][0])
		Bp += (1. - 0.5*self.Ls[0][1] - 0.25*w0 - 1.5*self.cA[0]*np.sqrt(-1*S0))
		Bpp = polyval2d(That,Phat,polyder(polyder(self.Bg,axis=1),axis=1))
		Btt = polyval2d(That,Phat,polyder(polyder(self.Bg,axis=0),axis=0))
		Btp = polyval2d(That,Phat,polyder(polyder(self.Bg,axis=0),axis=1))

		#response functions
		v = self.v = Bp + Vs + Wp*x*(1-x) # volume in dim units
		kt = self.kt = (-1./v)*(Bpp + x*tau*Lpp + tau*Wpp*x*(1-x) + (tau*Lp + tau*Wp*(1.-2*x))*dxdp) #compressibility
		cp = self.cp = -tau*(Btt + 2*x*Lt + 2*x*(1.-x)*Wt + tau*x*Ltt + tau*Wtt*x*(1-x) + (tau*Lt + tau*(1.-2.*x)*Wt)*dxdt) #heat capacity
		ap = self.ap = (1./v)*(Btp + x*Lp + x*tau*Ltp + Wp*x*(1.-x) + tau*Wtp*x*(1.-x) + (tau*Lp + tau*Wp*(1.-2.*x))*dxdt)

		self.compute = False


	def getEquibConc(self,params):
		"""params contains L, W, tau, and deltaP"""
		if params[0] == 0 and params[3] < 0:
			return 0.5
		try:
			res = minimize_scalar(self.f,args=params[0:3],bounds=(0.0001,0.9999),method='bounded')
			return res.x
		except:
			print("Unable to converge calculation for equilibirum conc.")
			exit(0)

	def setValues(self, x):
		L0, a, b, d, f = x[0:5]
		self.Ls[0][1], self.Ls[0][2] = a*L0, d*L0
		self.Ls[1][0], self.Ls[1][2] = L0, b*L0
		self.Ls[2][0] = f*L0

		w0 = x[5:6]
		delw = self.delw
		self.Ws[0][0], self.Ws[0][1] = 2.*delw, w0*delw
		self.Ws[1][0], self.Ws[1][1] = 2.*(1-delw), (1-delw)*w0

		bg = x[6:]
		self.Bg[0][2], self.Bg[1][1], self.Bg[2][0] = bg[0], bg[1], bg[2]
		self.Bg[0][3], self.Bg[1][2], self.Bg[2][1] = bg[3], bg[4], bg[5]
		self.Bg[3][0], self.Bg[1][3], self.Bg[2][2] = bg[6], bg[7], bg[8]
		self.Bg[3][1], self.Bg[4][0], self.Bg[0][4] = bg[9], bg[10],bg[11]


	def residual(self,x,data):
		res = 0
		self.setValues(x)
		for t,p,rho in data:
			res += (rho - self.getDensity(t,p))**2
		print(res)
		return res



from scipy.optimize import minimize, minimize_scalar
#implement constraints
# constraints : dict or sequence of dict, optional
#Constraints definition (only for COBYLA and SLSQP). Each constraint is defined in a dictionary with fields:
class ThermoFit:
	def __init__(self,model,minimizer=minimize):
		assert(isinstance(model,ThermoPPty))
		self.model = model
		self.minimizer = minimizer

	def runAll(self):
		pass

	def runPerIsochore(self):
		pass

	def runIsochore(self,Density,x0,method="Nelder-Mead"):
		if Density not in self.model.IsochoreDict:
			print("Isochore chosen is unavailable.")
		else:
			data = []
			for t in self.model.IsochoreDict[Density]:
				data.append([t,self.model.IsochoreDict[Density][t],Density])
			res = self.minimizer(self.model.residual,x0,
								args=(data),
								method=method,
								options={'disp':True,
								'maxfev':1000,
								'disp':True})
			np.savetxt('current_fit',res.x)



if __name__ == "__main__":
	myThermoPPty = ThermoPPty(fn='./Data/pruned025compiledEOSAll',Tc=213,Pc=338.,Rhoc=0.0595)
	start = False
	nLs, nWs , nBg, delw, f = 5, 1, 12, 1., 0.
	x=np.zeros(nLs+nWs+nBg)
	if start:
		x[0:nLs]=1.55607,0.154014,0.125093,0.00854418, 1.14576
		x[nLs:nLs+nWs]=0.03
		x[nLs+nWs:nLs+nWs+nBg] = [-0.00261876, 0.257249, -6.30589, 0.000605678,
		0.0248091, -0.0400033, 2.18819, -0.000994166, -0.00840543, 0.0719058, -0.256674, 0.]
	else:
		x=np.loadtxt('./Data/current_fit_x2_06')

	myThermoPPty.setValues(x)

	print(myThermoPPty.__str__(T=273.,P=101.))
	print(myThermoPPty.UniqueDensities)

	tfit = ThermoFit(myThermoPPty)
	tfit.runIsochore(Density=880.,x0=np.copy(x))

	print(myThermoPPty.__str__(T=273.,P=101.))







