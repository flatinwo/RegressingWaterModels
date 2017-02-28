import numpy as np 

""" A simple function to rougly determine
if a data point is relaxed based on fskts
author: Folarin Latinwo
date: 02/17/17
"""

def pruneFskts(mydir="/Users/Folarin/Research/UMD/Analysis/EvanAllModels/TIP5P/fsktandk/",
				prefix="compiledFskt", suffix=".dat", datarange=[], cutofftime=200000.,
				cutoffF=np.exp(-1),scale=1000.):
	prunelist = []
	for density in datarange:
		filename = mydir+prefix+str(density)+suffix
		x,y,z = np.loadtxt(filename,usecols=(0,1,2),unpack=True)
		FsktDict = {}
		for t,F,T in zip(x,y,z):
			if T not in FsktDict: FsktDict[T] = []
			FsktDict[T].append([t,F])
		for T in FsktDict:
			if min(FsktDict[T],key= lambda x: abs(x[0]-cutofftime))[1] > cutoffF:
				prunelist.append([float(density)*scale,T])
	return prunelist


if __name__ == "__main__":
	import TSEOS as ts
	rd = ts.RawData()
	#rd.plotAll()

	df = ["%.2f"%(0.86 + 0.02*i) for i in range(10)]
	prunelist = pruneFskts(datarange=df,cutofftime=200000)
	for density,T in prunelist: rd.trimAlongIsochore(Density=density,Tl=T)
	#rd.plotAll(filename="./Data/PrunedEOS_050ns.pdf",savefig=True)
	rd.writeData(filename=r"./Data/pruned200CompiledEOSAll")