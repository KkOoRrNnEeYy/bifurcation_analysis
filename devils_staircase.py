import numpy as np
from itertools import product
import matplotlib.pyplot as plt


class PSI:
    def __init__(self, n=2, k=4, g=10):
        self.n = n
        self.k = k
        self.g = g
        self.Dk = self.get_Dk(1)
        self.psi_vals = np.array([self.psi(dk) for dk in self.Dk])
        
    def get_alpha(self, p):
        alpha = 0
        r = 1
        while True:
            inc = self.g**(-(p-1) * (self.n**r-1)/(self.n-1))
            if inc < 10**(-4):
                return alpha + inc
            alpha += inc
            r += 1

    def dev_g(self, arr):
        return [arr[i]/self.g**(i+1) for i in range(len(arr))]
    def get_Dk(self, extend=False):
        Is = list(map(np.array, [i for i in product(range(0, self.g), repeat = self.k)]))
        Dk = np.array([np.sum(self.dev_g(I)) for I in Is])
        if extend:
            Dk = np.array([*Dk, *(Dk+1)])
        return Dk
    def get_ir(self, dk):
        formatter = f":.{self.k}f"
        formatter = ("{" + formatter + "}").format
        dk = formatter(dk)
        I = [int(i) for i in dk[2:]]
        while len(I) < self.k:
            I.append(0)
        return I
    
    def i1(self, I, r):
        #<i>
        if r == 1:
            return 0
        else:
            if I[r-1] < self.g-1:
                return 0
            else:
                return 1
    def i2(self, I, r):
        #[i]
        if r == 1:
            return 0
        else:
            if I[r-1] < self.g-2:
                return 0
            else:
                return 1
    def i_tilde(self, I, r):
        return I[r-1] - (self.g-2)*self.i1(I, r)
    def p(self, r, mr):
        return (self.n**(r-mr)-1)/(self.n-1)
    def m(self, I, r):    
        S = 0
        for s in range(1,r):
            p = 1
            for l in range(s, r):
                p *= self.i2(I, l)
            S += p
        mr = self.i1(I, r) * (1+S)
        return mr
    
    def psi(self, dk):
        add = 0
        if dk >= 1: add = 1
        I = self.get_ir(dk)
        S = 0
        for r in range(1, self.k+1):
            it = self.i_tilde(I, r)
            mr = self.m(I, r)
            S += it * 2**(-mr) * self.g**(-self.p(r, mr))
        return np.round(S, 20) + add
    
    def der_psi(self, dk, der_degree=1):
        h = 10**(-self.k)
        if der_degree == 1:
            y1 = self.psi(dk)
            y2 = self.psi(dk+h)
            return (y2-y1)/h
        else:
            y1 = self.der_psi(dk, 1)
            y2 = self.der_psi(dk+h, 1)
            return (y2-y1)/h
        
    def inv_psi(self, psi_val):
        ind = np.abs(self.psi_vals - psi_val).argmin()
        return self.Dk[ind]
    
    def der_inv_psi(self, psi_val, der_degree=1):
        ind = np.abs(self.psi_vals - psi_val).argmin()
        return (self.Dk[ind+1]-self.Dk[ind])/(self.psi_vals[ind+1] - self.psi_vals[ind])
  
    def plot_psi(self, Dk, y_psi):
        n = len(y_psi)
        fig, axs = plt.subplots(1,n,figsize=(7*n, 6))
        fig.suptitle(f'k={self.k}, g={self.g}')
        if n == 1:
            y = y_psi[0]
            axs.plot(Dk[:len(y)], y, c='black')
            axs.set(ylabel=r"$\psi$", xlim=(0, 1), ylim=(0, 1), xticks=np.arange(11)/10)
        else:
            for i in range(n):
                y = y_psi[i]
                axs[i].plot(Dk[:len(y)], y, c='black')
                axs[i].set(ylabel=r"$\psi$" + "'"*i, xlim=(0, 1), xticks=np.arange(11)/10)
        plt.savefig(f'psi_plots\k={self.k}_g={self.g}_psi.png')
        plt.show()
        
        

