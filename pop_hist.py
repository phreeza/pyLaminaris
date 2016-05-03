import numpy as np
from matplotlib import pyplot as plt

y_sep = [0.,500.,500.,1000.]
x_sep = [20,10,5]
h = np.zeros(40000)
for k in range(300000):
    p = [[[5*(k-1),0.] for m in range(2**n)] for n in range(4)]
    p[0][0][1] = 4000.+ 200*np.random.randn()
    #plt.plot([5*(k-1),5*(k-1)],[2000,p[0][0][1]],color=c[k],lw=lw)
    h[int(77*np.random.random()):int(p[0][0][1]):77] += 1
    for n in range(1,4):
        for m in range(2**n):
            if m%2 == 0:
                p[n][m+1][1] = p[n][m][1] = p[n-1][0][1] + y_sep[n-1] + 200*np.random.randn()
                #plt.plot([p[n-1][m/2][0],p[n-1][m/2][0]],[p[n][m][1],p[n-1][m/2][1]],color=c[k],lw=lw)
                h[int(p[n-1][m/2][1]):int(p[n][m][1]):77] += 1
            if n == 3:
                lm = y_sep[n] + 200*np.random.randn()
                #plt.plot([p[n][m][0],p[n][m][0]],[p[n][m][1],p[n][m][1]-lm],color=c[k],lw=lw)
                h[int(p[n][m][1]):int(p[n][m][1]+lm):77] += 1
                #plt.plot([p[n][m][0]],[p[n][m][1]-lm],color=c[k],lw=lw,marker='o',mec='none',ms=3.5)

