import numpy as np
import OneDFit as od
import TSEOS as ts
from scipy.optimize import brentq, minimize_scalar, leastsq
from numpy.polynomial.polynomial import polyval,polyval2d,polyder
import pandas as pd

def getXs(cxy,x,y):
	f   = polyval2d(x,y,cxy)
	fx  = polyval2d(x,y, polyder(cxy,axis=0))
	fy  = polyval2d(x,y, polyder(cxy,axis=1))
	fxx = polyval2d(x,y, polyder(polyder(cxy,axis=0),axis=0))
	fxy = polyval2d(x,y, polyder(polyder(cxy,axis=0),axis=1))
	fyy = polyval2d(x,y, polyder(polyder(cxy,axis=1),axis=1))
	return f,fx,fy,fxx,fxy,fyy


class ThermoPPty:
	def __init__(self,fn='./Data/compiledEOSAll', Tc=200, Pc=200, Rhoc=1.,Tl=230.,spinodal=None):
		if not spinodal:
			sf = od.SpinodalFit(fn,Tl) #Maybe i don't need this as a member fxn
			sf.setCriticalPt(Tc,Pc,Rhoc) #But ensures I have the same criticalpt
			sf.fit()
		else:
			assert isinstance(spinodal,od.SpinodalFit)
			sf = spinodal

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
		self.IsochoreDictCv = {}


		self.residual_scale = {}

		for density in sf.Rawdata.IsochoreDict:
			if density not in self.IsochoreDict: self.IsochoreDict[density] = {}
			for T, P in sf.Rawdata.IsochoreDict[density]:
				self.IsochoreDict[density][T] = P
		
		self.f = lambda x, p : p[0]*x + p[1]*(x*np.log(x) + (1-x)*np.log(1-x)) + p[2]*x*(1-x)


	def addCvData(self,fn):
		# only add Cv data in ThermoSet
		df = pd.read_csv(fn,sep='\t')
		grouped_df = df.groupby('Density')
		for key,item in grouped_df:
			if key not in self.IsochoreDictCv: self.IsochoreDictCv[key] = {}
			for T, Cv in zip(item['Temperature'],item['Cv']):
				self.IsochoreDictCv[key][T] = Cv

		self.residual_scale['Cv'] = 1.




	def getPressure(self,T=200.,Density=1000, override_guess=False, Pguess=0., dRho=0.):
		my_nearest = lambda array, value : array[np.abs(array-value).argmin()]
		RhoMax = my_nearest(self.UniqueDensities,Density+dRho)
		
		Tn = my_nearest(np.array(list(self.IsochoreDict[RhoMax].keys())),T+10.)
		Pmax = self.IsochoreDict[RhoMax][Tn]
		myres = lambda P: Density - self.getDensity(T,P,False)
		if override_guess: Pmax = Pguess
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
		Pshat,S0 = self.criticalpt.reduceP(np.array([Ps-20.,self.cS[0]]))
		PminusPhat = Phat - Pshat

		At, Att = polyval(That,polyder(self.cA)), polyval(That,polyder(polyder(self.cA)))
		Pst, Pstt = polyval(That,polyder(self.cS)), polyval(That,polyder(polyder(self.cS)))
		K2 = self.criticalpt.RhocRTc

		#print(PminusPhat)


		if PminusPhat < 0: PminusPhat=0.
		if S0 > 0: print("there was a numerical error with spinodal")
		
		Vs = 1.5*A*np.sqrt(PminusPhat) #spinodal volume... write for case without spinodal
		kts = 0.75*A/np.sqrt(PminusPhat) #spinodal contribution to isothermal compressibility


		cps = Pst*(-2.*Vs*At/(K2*A) + kts*Pst/(K2*K2)) + Att*(PminusPhat)**(1.5) - Vs*Pstt/K2 # d2GsdTs
		aps = 1.5*At*np.sqrt(PminusPhat) - kts*Pst/K2 # correction to thermal expansivity

		# get reaction term, that is -lnK or in this case L
		L, Lt, Lp, Ltt, Ltp, Lpp = getXs(self.Ls,That,Phat)
		W, Wt, Wp, Wtt, Wtp, Wpp = getXs(self.Ws,tau,Phat)

		params = np.array([L,tau,W,Phat])
		x = self.x = self.getEquibConc(params)


		#x = 0.93127849869435763
		dxdt = self.dxdt = (Lt + np.log(x/(1-x)) + Wt*(1.-2*x))/(2*W - tau/(x*(1.0-x)))
		dxdp = self.dxdp = (Lp + Wp*(1-2*x))/(2*W - tau/(x*(1-x)))


		#dGA/dP
		Bp = polyval2d(That,Phat,polyder(self.Bg,axis=1))
		w0 = 2.0*(self.Ws[0][1]/self.Ws[0][0])
		Bp += (1. - 0.5*self.Ls[0][1] - 0.25*w0 - 1.5*self.cA[0]*np.sqrt(-1*S0))
		Bpp = polyval2d(That,Phat,polyder(polyder(self.Bg,axis=1),axis=1))
		Btt = polyval2d(That,Phat,polyder(polyder(self.Bg,axis=0),axis=0))
		Btp = polyval2d(That,Phat,polyder(polyder(self.Bg,axis=0),axis=1))


		#print((Bpp,Lpp,Wpp,Lp,Wp,kts))
		#print((cps,Btt,Ltt,Wtt,Lt,Wt,dxdt))

		#response functions
		v = self.v = Bp + Vs + x*Lp + Wp*x*(1-x) # volume in dim units
		kt = self.kt = (-1./v)*(Bpp + x*Lpp + Wpp*x*(1-x) + (Lp + Wp*(1.-2*x))*dxdp + kts) #compressibility
		cp = self.cp = -tau*(cps + Btt + x*Ltt + x*(1.-x)*Wtt + (Lt + np.log(x/(1.-x)) + (1.-2.*x)*Wt)*dxdt) #heat capacity
		ap = self.ap = (1./v)*(aps + Btp + x*Ltp  + Wtp*x*(1.-x) + (Lp + Wp*(1.-2.*x))*dxdt)



		#print((T,P,cp,v))
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

	def setValues(self, x,interaction=True,index=None):
		if interaction:
			L0, a, b, d, f = x[0:5]
			self.Ls[0][1], self.Ls[0][2] = a*L0, d*L0
			self.Ls[1][0], self.Ls[1][1] = L0, b*L0
			self.Ls[2][0] = f*L0

			w0 = x[5:6]
			delw = self.delw
			self.Ws[0][0], self.Ws[0][1] = 2.*delw, w0*delw
			self.Ws[1][0], self.Ws[1][1] = 2.*(1-delw), (1-delw)*w0

		bg = x[6:]
		if index is not None:
			self.Bg[0][2], self.Bg[1][1], self.Bg[2][0] = bg[0], bg[1], bg[2]
			self.Bg[0][3], self.Bg[1][2], self.Bg[2][1] = bg[3], bg[4], bg[5]
			self.Bg[3][0], self.Bg[1][3], self.Bg[2][2] = bg[6], bg[7], bg[8]
			self.Bg[3][1], self.Bg[4][0], self.Bg[0][4] = bg[9], bg[10],bg[11]
		else:
			self.Bg[0][2],self.Bg[1][1],self.Bg[1][2] = bg[0],bg[1],bg[4]
			self.Bg[2][1],self.Bg[1][3]	= bg[3],bg[7]


	def residual(self,x,data,cv_flag=False):
		self.setValues(x)
		if not cv_flag:
			res1 = 0.
			norm1 = 0.
			for t,p,rho in data:
				res1 += (rho - self.getDensity(t,p))**2
				norm1 += rho**2
			res = res1/norm1
			print(res)
			return res
		else:
			res1,res2 = 0.,0.
			norm1,norm2 = 0.,0.
			for t,p,rho,cv in data:
				res1 += (rho - self.getDensity(t,p))**2
				if abs(cv) > 0.00001: 
					res2 += (cv - self.getCv(t,self.IsochoreDict[rho][t],computed=False))**2
				norm1 += rho**2
				norm2 += cv**2
				#print((t,p,rho,cv,res1,res2))
				
			res = 2.*(res1/norm1) + (res2/norm2)
			print((res,res1/norm1,res2/norm2))
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

	def runIsochore(self,Densities,x0,method="Nelder-Mead",figure=True,holdFigure=False,predict=False,
		cv_flag=False):
		if Densities[0] not in self.model.IsochoreDict:
			print("Isochore chosen is unavailable.")
		else:
			data = []
			T_obs = []
			P_obs = []
			for Density in Densities:
				for t in self.model.IsochoreDict[Density]:
					if not cv_flag: data.append([t,self.model.IsochoreDict[Density][t],Density]) 
					else: data.append([t,self.model.IsochoreDict[Density][t],Density,
						self.model.IsochoreDictCv[Density][t]])
					T_obs.append(t)
					P_obs.append(self.model.IsochoreDict[Density][t])
			
			if not predict:
				def compressibility_constraint(x,data):
					kts = np.zeros(len(data))
					self.model.setValues(x)
					for i in range(len(kts)):
						kts[i] = self.model.getCompressibility(data[i][0],data[i][1])
					#print("ktc called.")
					return kts

				def cv_constraint(x,data):
					cvs = np.zeros(len(data))
					self.model.setValues(x)
					for i in range(len(cvs)):
						t,p,rho = data[i][0],data[i][1],data[i][2]
						cvs[i] = self.model.getCv(t,p)
					#print("cvc called")
					return cvs

				def rho_constraint(x,data):
					return np.array([self.model.getDensity(t,p) for t,p,rho,cv in data])



				if method == "SLSQP" or method == "COBYLA":
					rho_cons = ({'type':'ineq',
							'fun': rho_constraint, 
							'args': (data,)
						})

					kt_cons = ({'type':'ineq',
							'fun': compressibility_constraint, 
							'args': (data,)
						})
					cv_cons = ({'type':'ineq',
							'fun': cv_constraint, 
							'args': (data,)
						})
					if cv_flag: 
						cons = (rho_cons,kt_cons, cv_cons)
					else: 
						cons = (rho_cons,kt_cons)
					res = self.minimizer(self.model.residual,x0,
								args=(data,cv_flag),
								method=method,
								constraints=cons,
								options={'disp':True,'tol':1e-8,
								'maxiter':1000})
				else:
					res = self.minimizer(self.model.residual,x0,
								args=(data,cv_flag),
								method=method,
								options={'disp':True,
								'maxiter':1000})

				np.savetxt('./params/current_fit_sub',res.x)
				self.model.setValues(res.x)
			if figure:
				import matplotlib.pyplot as plt
				if predict: self.model.setValues(x0)
				T = np.linspace(213.,350.,200)
				for Density in Densities:
					P_pred = []
					for t in T:
						P_pred.append(self.model.getPressure(T=t,Density=Density,dRho=0.))
					plt.plot(T,P_pred,'-')
				plt.plot(T_obs,P_obs,'o')
				if not holdFigure: plt.show()
				if cv_flag:
					T = np.linspace(213.,350.,200)
					for Density in Densities:
						print(Density)
						cv_pred = []
						for t in T:
							p = self.model.getPressure(t,Density)
							cv_pred.append(self.model.getCv(t,p))
						plt.plot(T,cv_pred,'_')
						plt.plot(list(self.model.IsochoreDictCv[Density].keys()),
							list(self.model.IsochoreDictCv[Density].values()),'o')
						plt.show()





if __name__ == "__main__":
	tc, pc, rhoc = 213., 340., 0.059
	mys = od.SpinodalFit("./Data/compiledEOSAll",230.)
	mys.setCriticalPt(tc,pc,rhoc)
	mys.fitParams["A(T) params"] = np.array([  5.62586691e-06,  -2.43215926e-05,   1.70903857e-04,
        -2.05070945e-04])[::-1]
	mys.fitParams["P(T) params"] = np.array([ -496.40292951,  1503.16997588, -1632.66291085])[::-1]
	myThermoPPty = ThermoPPty(fn='./Data/pruned025compiledEOSAll',Tc=tc,Pc=pc,Rhoc=rhoc)#,spinodal=mys)
	#myThermoPPty.addCvData(fn='./Data/compiledEOSCV50')
	myThermoPPty.addCvData(fn='./Data/collatedCV')

	start = False
	nLs, nWs , nBg, delw, f = 5, 1, 12, 1., 0.
	x=np.zeros(nLs+nWs+nBg)
	if start:
		x[0:nLs]=1.55607,0.154014,0.125093,0.00854418, 1.14576
		x[nLs:nLs+nWs]=0.03
		x[nLs+nWs:nLs+nWs+nBg] = [-0.00261876, 0.257249, -6.30589, 0.000605678,
		0.0248091, -0.0400033, 2.18819, -0.000994166, -0.00840543, 0.0719058, -0.256674, 0.]
	else:
		#x=np.loadtxt('./Data/current_fit_x2_06')
		x=np.loadtxt('./params/current_fit')
		#x=np.loadtxt('./best_current_fit_3')

	#x = np.zeros_like(x)


	print(myThermoPPty.criticalpt)
	myThermoPPty.setValues(x,interaction=True,index=True)
	y = np.loadtxt('./params/current_fit_sub')
	myThermoPPty.setValues(y,interaction=True)

#290.0, -177.96899999999999, 860.0, 7996.25, 3458.3573345851223, nan

	print(myThermoPPty.__str__(T=273.,P=101.))
	print(myThermoPPty.__str__(T=290.,P=-177.97))

	densities = myThermoPPty.UniqueDensities
	#densities = [1120.]

	tfit = ThermoFit(myThermoPPty)
	#for density in densities:
	#	if density != 1120.:
	#		tfit.runIsochore(Density=density,x0=np.copy(x),holdFigure=True,predict=True)

	#tfit.runIsochore(Density=920.,x0=np.copy(x),holdFigure=True,predict=True)
	#tfit.runIsochore(Density=1020.,x0=np.copy(x),holdFigure=True,predict=True)
	tfit.runIsochore(Densities=densities,x0=np.copy(y),holdFigure=False,predict=True,method="Nelder-Mead",cv_flag=True,figure=True) #SLSQP is great for 
	#rapid search



	print(myThermoPPty.__str__(T=273.,P=101.))







