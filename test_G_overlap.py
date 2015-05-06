import glob
from sys import stderr
import numpy as np

import time
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mCols
import matplotlib.cm as cm
import matplotlib.mlab as mlab

from scipy.stats import norm

def pw_GMM_overlap(self, gmm):
    
    overlaps = []       
    mus = gmm.means[:,0]
    vars = np.array([v[0][0] for v in gmm.covars])
    weights = np.array(gmm.weights)
    order = np.argsort(mus)
    mus = mus[order]
    weights = weights[order]
    vars = vars[order]
    
    for i in xrange(mus.shape[0]-1):
        mu1 = mus[i]
        mu2 = mus[i+1]
        v1 = vars[i]
        v2 = vars[i+1]
        sd_max = np.sqrt(max(v1, v2))
        mn = min(mu1, mu2) - 10*sd_max
        mx = max(mu1, mu2) + 10*sd_max
        xs = np.arange(mn,mx,0.01)
        o = np.sum(np.min(np.c_[norm.pdf(xs,loc=mu1,scale=v1)*weights[i], 
                                norm.pdf(xs,loc=mu2,scale=v2)*weights[i+1]], 1))  * 0.01
        overlaps.append([o,o/weights[i],o/weights[i+1]])
    return overlaps

def eval_G(G, x):

    u, v = G
    s = np.sqrt(v)
    sq2pi = np.power(2*np.pi,0.5)

    y = (1/(s*sq2pi)) * np.exp( -1*((x-u)*(x-u))/(2*s*s) )
    #y = mlab.normpdf(x, mu, s)
    return y

def get_intersection(G1, G2, ws, tol=0.01):
    #sort so G1.mu < G2.mu
    #ui < uj
    oGs = [G1, G2] 
    ows = ws
    Gs, ws = [], []
    args = np.argsort([G1[0],G2[0]])
    
    for i in args:
        Gs.append(oGs[i])
        ws.append(ows[i])
    ui, vi = Gs[0]
    uj, vj = Gs[1]
    si, sj = np.sqrt(vi), np.sqrt(vj)
    al, be = ws
    print ui, si, uj, sj
    
    if si == sj:
        x=(ui+uj)/2.0
    else:
        sq2pi = np.power(2*np.pi,0.5)
        c = (2*si*si*sj*sj) * ( np.log( al/(si*sq2pi) ) - np.log( be/(sj*sq2pi) ) )
        c = c  + (si*si*uj*uj)-(sj*sj*ui*ui)
        b = -((2*uj*si*si)-(2*ui*sj*sj))
        a = (si*si)-(sj*sj)
        
        q=(b**2 - 4*a*c)
        if q<0: 
            x=None
        else:
            x1 = (-b + np.sqrt(q)) / (2*a)
            x2 = (-b - np.sqrt(q)) / (2*a)
            
            x=x1
            if (x1 < ui and x1 < uj) or (x1 > ui and x1 > uj):
                x=x2
    
    if x==None:
        return None, None, None, None

    y = al*eval_G(G1, x) 

    mn = ui - 5*si
    mx = uj + 5*sj
    xis = np.arange(x,mx, tol)
    xjs = np.arange(mn,x, tol)

    i_integral = np.sum(mlab.normpdf(xis, ui, si)*al)*tol
    j_integral = np.sum(mlab.normpdf(xjs, uj, sj)*be)*tol
    overlap = i_integral+j_integral

    return x, y, overlap/al, overlap/be

def plot_G(ax, Gs, weights, intersect):
    
    G_x=np.arange(0,5,.001)
    l = len(Gs)
    G_ys = []
    for i in xrange(l):
        c = cm.hsv(float(i)/l,1)
        mu = Gs[i][0]
        var = Gs[i][1]
        
        G_y = mlab.normpdf(G_x, mu, var**.5)*weights[i]
        G_ys.append(G_y)
        ax.plot(G_x,G_y,color=c)
        ax.plot(mu,-.001,"^",ms=10,alpha=.7,color=c)
    
    #ax.plot(G_x,np.power(G_ys[1]-G_ys[0],1),color='k')
    if intersect[0] !=None:
        ax.plot(intersect[0],intersect[1],"|",ms=10,alpha=.7,color='k')
    ax.plot([0,5],[0,0],color='k')

if __name__=="__main__":

    plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')

    fig, axarr = plt.subplots(2, 1)
    fig.set_figwidth(11)
    fig.set_figheight(8.5)
    axescolor  = '#f6f6f6'
        
    #axarr[0].plot(cps, sunk_cps, 'ro', alpha=0.2)
    #axarr[0].set_xlim(-0.10,max(cps)+1)
    #axarr[0].set_ylim(-0.10,max(sunk_cps)+1)
    Gs = [[0.8046310969121766,0.025287736617867537], [1.7864071995597282,0.10703174214076343]]
    weights = [.1,.9]
    Gs = [[2,.1], [3,.1]]
    weights = [.5,.5]
    
    Gs = [[1.96,.3*.3], [2.8,.17*.17]]
    weights = [.85,.15]

    Gs = [[0,1.0000000000000001e-05*1.0000000000000001e-05], [0.015190999999999967,1.000000002439455e-05*1.000000002439455e-05]]
    weights = [0.94444444444444153,1-0.94444444444444153]
    
    Gs = [[2.9674184163362765,0.19997407419512075*0.19997407419512075], [3.5563537847950277,0.011010508240214663*0.011010508240214663]]
    weights = [0.98735838860058134,1-0.98735838860058134]
    
    ix,iy, fi, fj = get_intersection(Gs[0], Gs[1], weights, tol=0.001)
    print ix, iy, fi, fj

    plot_G(axarr[0], Gs, weights, [ix,iy])   
    #axarr[0].set_xlim(0,.1)
    fig.savefig("./test_overlap.png")
    plt.close()



