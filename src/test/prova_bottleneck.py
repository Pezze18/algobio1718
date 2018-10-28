import time
import numpy as np
import bottleneck as bottle

a=np.asarray([0,1,2,4,5,6,7,8,9,10,13,34])
thresholds=np.asarray([0,0.000001,3,7,29,50])
indices=np.digitize(a,thresholds)-1
print(indices)
lista=[]
print(indices==1)

lista.append(a[indices==1])
print(lista)
raise ValueError


a=np.random.rand(100000000)*10
#a=list(range(90,100))+list(range(40,50))+list(range(0,5))+list(range(10,19))
#a=np.asarray(a)
#print(a)
print()
t_start = time.time()
b=bottle.partition(a, 10)[:10]
b=np.sort(b)
t_end = time.time()
print("Elapsed time: ", time.strftime("%H:%M:%S", time.gmtime(t_end - t_start)))
print(b)
print()
t_start = time.time()
#print(np.sort(a)[:1])
print(np.min(a))
t_end = time.time()
print("Elapsed time: ", time.strftime("%H:%M:%S", time.gmtime(t_end - t_start)))
