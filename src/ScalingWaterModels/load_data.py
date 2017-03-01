import pickle
import numpy as np

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


	return AllData


if __name__ == "__main__":
	Data = load_data()

