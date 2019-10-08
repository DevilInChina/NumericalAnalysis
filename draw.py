import matplotlib.pyplot as plt 
import numpy as np
import sys
with open('data.txt', 'r') as f:
    data=f.readlines()
lists=[float(line)for line in data]

def fixedpointplot(origin,formula,x,title):
    y1=origin(x)
    plt.plot(x,y1)
    plt.grid(b=True,which='major',color='gray',linestyle='--')
    for funcs in formula:
        if(sys.argv[1]=='lag' or sys.argv[1]=='newtown'):
            y1=np.array(lists)
        else:
            y1 = funcs(x)
        y1=np.resize(y1,x.size)
        print(y1.size)
        plt.plot(x,y1)
        plt.grid(b=True,which='major',color='gray',linestyle='--')
    plt.savefig(title+'.png')
    plt.show()
    return 
def func2a(x):
	ret=0.0
	for i in reversed(lists):
		ret=x*ret+i
	return ret
def func(x):
    return 1/(25*x*x+1)

x=np.arange(-1.0,1.0,0.01)

funclis=[func2a]

fixedpointplot(func,funclis,x,'paint')

