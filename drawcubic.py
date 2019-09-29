import matplotlib.pyplot as plt 
import numpy as np
import sys
with open('data.txt', 'r') as f:
    data=f.readlines()
lists=[[float(num) for num in line.split(' ')]for line in data]

with open('dat.in','r') as f:
    data=f.readlines()
xy=[[float(num)for num in line.split(' ')]for line in data]

l = int(xy[0][0])
def fixedpointplot(origin,formula,x,title):
    y1=origin(x)
    plt.plot(x,y1)
    y2=formula(x)
    plt.plot(x,y2)
    plt.grid(b=True,which='major',color='gray',linestyle='--')
    plt.savefig(title+'.png')
    plt.show()
    return
def cubicFunc(ls,x):
    ret = 0
    for i in reversed(ls):
        ret=ret*x+i
    return ret
def func2a(x):
   newli=[]
   for xs in x:
        print(xs)
        cnt = 2
        if xs < xy[cnt][0]:
           newli.append(cubicFunc(lists[0],xs))
           print(0)
        else:
           cnt=cnt+1
           while cnt<l:
               if xs < xy[cnt][0]:
                   print(cnt-2)
                   newli.append(cubicFunc(lists[cnt-2],xs))
                   break
               cnt=cnt+1
           if(cnt==l):
               print(cnt-2)
               newli.append(cubicFunc(lists[cnt-2],xs))
   return newli
def func(x):
    return 1/(25*x*x+1)

x=np.arange(-1.0,1.0,0.01)


fixedpointplot(func,func2a,x,'paint')

