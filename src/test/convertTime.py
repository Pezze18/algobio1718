import sys

times=["0m11.588s", " 0m15.902s", "4m48.703s", "403m10.665s"]
times_trad=[]
for time in times:
    (time,seconds)=time.split("m")
    time=int(time)
    seconds=(float(seconds.split("s")[0]))

    ms=seconds-int(seconds)
    seconds=int(seconds)

    h=int(time/60)
    m=time-60*h
    if h==0:
        h="00"
    if m==0:
        m="00"
    if seconds==0:
        seconds="00"
    ms=str(ms*1000)[0:3]

    times_trad.append(str(h)+"h:"+str(m)+"m:"+str(seconds)+"s:"+str(ms)+"ms")
print(times_trad)
