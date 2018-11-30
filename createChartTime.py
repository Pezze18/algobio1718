import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#['o', '.', ',', 'x', '+', 'v', '^', '<', '>', 's', 'd']

#names=["oldVersion","noboundMatrix","noboundLevelsVec","nobound_migliorato"]
#names=["nobound_migliorato","bound_min","bound_min_migliorato","bound_min_migliorato_optimized"]
#names=["bound_min_migliorato", "bound_min_migliorato_optimized","bound_min_migliorato_iterations", "bound_min_migliorato_iterations_optimized"]
#names=["bound_order","bound_order_optimized","bound_order_improved_iterations_percentiles"]
#names=["bound_means","bound_means_optimized","bound_means_iterations","bound_means_iterations_optimized"]

#names=["nobound_migliorato", "bound_min_migliorato_iterations_optimized", "bound_order_improved_iterations_percentiles", "bound_means_iterations_optimized"   ]
#translator={"nobound_migliorato":"Enumerazione Completa", "bound_min_migliorato_iterations_optimized":"Limite dei minimi", "bound_order_improved_iterations_percentiles":"Limite dell'ordinamento", "bound_means_iterations_optimized":"Limite della media"}

#names=["det_old","det_numpy","det_LevelsVec"]
#translator={"det_old":"oldVersion", "det_numpy":"noboundMatrix","det_LevelsVec":"noboundLevelsVec"}

names=["oldVersion","onlyBFS","onlyBFSAndProbCover","BFSAndProbCoverAndComplement","BFSAndLevelsVec"]
translator={"oldVersion":"Versione Originale","onlyBFS":"BFS","onlyBFSAndProbCover":"BFS e Numpy", "BFSAndProbCoverAndComplement":"BFS e Numpy e Complemento","BFSAndLevelsVec":"BFS e Numpy e Complemento e Suddivisione"}

#names=["det_oldVersion","det_onlyBFS","det_onlyBFS","det_onlyBFSAndNumpy","det_BFSAndNumpyAndComplement","det_BFSAndLevelsVec"]
#translator={"det_oldVersion":"oldVersion","det_onlyBFS":"onlyBFS","det_onlyBFSAndNumpy":"onlyBFSAndNumpy",           "det_BFSAndNumpyAndComplement":"BFSAndNumpyAndComplement","det_BFSAndLevelsVec":"BFSAndLevelsVec"}

#names=["OriginalVersion BDDE", "BestVersion BDDE", "BestVersion Combinatorial", "OriginalVersion Combinatorial"]
#names=["BDDE senza tagli","BestVersion BDDE","BestVersion Combinatorial"]

strategy="combinatorial"#enumertae, combinatorial
translate=True
marker='--x'
columns = ["method"] + ["k=" + str(i) for i in range(1, 11)]

filepath = "data/summary_tables/times_"+strategy+".csv"

frame = pd.read_csv(filepath, sep=",")
#print(frame)

title="Compare Time"
plt.title(title)
for index, row in frame.iterrows():
    if row["method"] in names:
        x=[]
        y=[]
        for k in range(1,11):
            if( isinstance(row["k="+str(k)],type("str"))):
                #print(type(row["k="+str(k)]))
                #print(row["k=" + str(k)])
                (time, seconds) = row["k="+str(k)].split("m")
                time = int(time)
                seconds = (float(seconds.split("s")[0]))
                ms = seconds - int(seconds)
                seconds = int(seconds)
                h = int(time / 60)
                m = time - 60 * h

                seconds=seconds+m*60+h*3600

                y.append(seconds)
                x.append(k)
        if(translate):
            plt.plot(x, y, marker, label=translator[row["method"]])
        else:
            plt.plot(x,y,marker,label=row["method"])

plt.ylabel('Time (m)')
if strategy=="enumerate" or strategy=="combinatorial" or strategy=="compare":
    plt.yscale("log")
plt.xlabel("k")
plt.legend()
plt.savefig(title+".png")
plt.show()



