from numba import jit
from numpy import arange
import time

# jit decorator tells Numba to compile this function.
# The argument types will be inferred by Numba when function is called.
@jit
def sum2d_gpu(arr):
    M, N = arr.shape
    result = 0.0
    for i in range(M):
        for j in range(N):
            result += arr[i,j]
    return result

def sum2d_cpu(arr):
    M, N = arr.shape
    result = 0.0
    for i in range(M):
        for j in range(N):
            result += arr[i,j]
    return result

start=time.time()
a = arange(90000).reshape(300,300)
print(sum2d_gpu(a))
end=time.time()
print("Time GPU: "+str(end-start))



start=time.time()
a = arange(9000000).reshape(3000,3000)
print(sum2d_cpu(a))
end=time.time()
print("Time CPU: "+str(end-start))

