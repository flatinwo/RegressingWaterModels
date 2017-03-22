import pickle
import numpy as np
from GuggenToy import ToyModel, GuggenModel

def load_data(filename="../Data/a2_003_3.pkl"):
	Data = open (filename, 'rb')
	([[Tde, Pde], 
	[Tko, Pko], 
	[Tcp, Pcp], 
	[Tvl, Pxvl], 
	[Tsp1, Pxsp1], 
	[Tsp2, Pxsp2], 
	[Tllt, Pllt], 
	[Tllw, Pllw], 
	[Tc, Pc]]) = pickle.load(Data,encoding="latin1")

	AllData = {}
	AllData["TMD"] = np.array([Tde,Pde])
	AllData["KtM"] = np.array([Tko,Pko])
	AllData["CpM"] = np.array([Tcp,Pcp])
	AllData["VLT"] = np.array([Tvl,Pxvl])
	AllData["VLS"] = np.array([Tsp1,Pxsp1])
	AllData["LVS"] = np.array([Tsp2,Pxsp2])
	AllData["LLT"] = np.array([Tllt,Pllt])
	AllData["LLW"] = np.array([Tllw,Pllw])
	AllData["LLCP"] = np.array([Tc,Pc])
	AllData["VLCP"] = np.array([1.0,1.0])


	return AllData


def load_data_and_rescale(filename="../Data/a2_003_3.pkl",a2=-0.3,callback="delineateSpinodalandWidomLine",writeFile=None):
	Data = load_data(filename)
	ty = ToyModel(1.,1.,a2)
	gw = GuggenModel(ty)

	for ppty in Data:
		T,P = Data[ppty]
		try: 
			for t,p in zip(T,P): gw.addPoint(t,p,ppty,Source="Unknown",callback=callback)
		except:
			pass
			#gw.addPoint(T,P,ppty,Source="Unknown",callback="delineateSpinodalandWidomLine")
	
	ppties=["LLCP","VLCP"]
	for ppty in ppties:
		T,P = Data[ppty]
		gw.addPoint(T,P,ppty,Source="Unknown",callback=callback)

	AllData = {}
	for property in gw.PropertyDict:
		rT = [gw.PropertyDict[property][i][0] for i in range(len(gw.PropertyDict[property]))]
		rP = [gw.PropertyDict[property][i][1] for i in range(len(gw.PropertyDict[property]))]
		AllData[property] = np.array([rT,rP])

	if writeFile is not None:
		gw.writetoFile(filename=writeFile,isnew=True, add_critical_pt=False,offsetCriticalPt=1.,callback="rescale")

	return AllData


if __name__ == "__main__":
	Data1 = load_data()
	Data2 = load_data_and_rescale() 
	Data3 = load_data_and_rescale(callback="identity",writeFile="RescaledModel.dat")


