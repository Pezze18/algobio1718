import numpy as np

a=[]
my_values=[]
b=np.searchsorted(a, my_values)
print(b)



"""
vec=np.random.rand(100)*100
print(vec)
vec=np.sort(vec)
print(vec)
vec=vec[::-1]
print(vec)
"""

"""
x = np.array([6, 6.4, 6.0, 1.6])
bins = np.array([0.0, 1.0, 2.5, 4.0, 10.0])
inds = np.histogram(x,bins)[0]
print(inds)
"""

"""
a=[6, 6.4, 6.0, 1.6]
print(sorted(a))
print(a)
print(a.sort())
print(a)
"""