import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#['o', '.', ',', 'x', '+', 'v', '^', '<', '>', 's', 'd']
strategy = "combinatorial"
#names=["oldVersion","noboundMatrix","noboundLevelsVec","nobound_migliorato"]
#names=["nobound_migliorato","bound_min","bound_min_migliorato"]
#names=["bound_min_migliorato", "bound_min_migliorato_optimized","bound_min_migliorato_iterations", "bound_min_migliorato_iterations_optimized"]
#names=["bound_order","bound_order_optimized","bound_order_improved_iterations_percentiles"]
#names=["bound_means","bound_means_optimized","bound_means_iterations","bound_means_iterations_optimized"]
names=["nobound_migliorato", "bound_min_migliorato_iterations_optimized", "bound_order_improved_iterations_percentiles", "bound_means_iterations_optimized"   ]

marker='--x'

columns = ["method"] + ["k=" + str(i) for i in range(1, 11)]
filepath = "data/summary_tables/times_enumerate.csv"

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
        plt.plot(x,y,marker,label=row["method"])

plt.ylabel('Time (m)')
plt.yscale("log")
plt.xlabel("k")
plt.legend()
plt.savefig(title+".png")
plt.show()



