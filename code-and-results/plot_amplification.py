## plotfile for amplification factors in 1D
import numpy as np
import matplotlib.pyplot as plt

## functions for amplification factors

def A_be(F,rho): # backward euler
    return 1/(1+4*F*np.sin(rho)**2)

def A_fe(F,rho): # forward euler
    return 1-4*F*np.sin(rho)**2

def A_cn(F,rho): # crank nicolson
    return (1-2*F*np.sin(rho)**2)/(1+2*F*np.sin(rho)**2)

def A_e(F,rho): # analytical amplification factor
     return np.exp(-4*F*rho**2)

# rho = k*dx
rho = np.linspace(0,1.5,1000)

F = [2,0.1]

#subplot or subplots
fig, ax = plt.subplots(1,2, figsize = (10,5))

markers_on = []
for i in range(10):
    markers_on.append(i*100)




# Plot for different F values
for i in range(len(F)):
    ax[i].set_title("F:{:.1f}".format(F[i]),fontsize=14)
    ax[i].plot(rho, A_be(F[i],rho),"-*",label = "BE", markevery = markers_on)
    ax[i].plot(rho, A_fe(F[i],rho),'-1', label = "FE",markevery = markers_on)
    ax[i].plot(rho, A_cn(F[i],rho), "-2",label = "CN",markevery = markers_on)
    ax[i].plot(rho, A_e(F[i],rho), label = "Exact")
    ax[i].set(xlabel = r"$\rho = k \Delta x$", ylabel = r"A($\rho$)")
    ax[i].tick_params(axis="x", labelsize=12)

ax[0].set(ylim =(-1,1))
ax[1].legend(fontsize=13)

plt.tight_layout()

## save

plt.show()
