import time
import numpy as np
import bottleneck as bottle

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
