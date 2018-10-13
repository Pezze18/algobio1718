import sys

time="745m11.742s"
(time,seconds)=time.split("m")
time=int(time)
seconds=int(float(seconds.split("s")[0]))

h=int(time/60)
m=time-60*h 
print(str(h)+":"+str(m)+":"+str(seconds))
