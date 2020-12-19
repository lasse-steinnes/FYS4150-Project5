## Python file for time and error analysis ##
import numpy as np
import matplotlib.pyplot as plt

# time 1D cases
timeFE = np.log10(np.array([0.511504,2.93544,2457.85]))
timeBE = np.log10(np.array([0.538446,3.9505,2786.51]))
timeCN = np.log10(np.array([0.327743,4.0914,3019.07]))

dx = np.array([0.2,0.1,0.01])


# plot time on a log scale
plt.figure()
plt.plot(dx,timeFE,'-x', label = "FE")
plt.plot(dx,timeBE,'-x', label = "BE")
plt.plot(dx,timeCN,'-x', label = "CN")
plt.xlabel('dx',fontsize = 13)
plt.ylabel('log10(runtime) [ms]',fontsize = 13)
plt.legend()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("./results/figures/analysis/time-analysis.png")
plt.show()


# read off files to check error by
# comparing with analytical expressions (1D case)
Nx1 = 5; Nt1 =  75

Nx2 = 10; Nt2 = 300

Nx3 = 100; Nt3 = 30000

t = 1 # total time numerical methods

x1 = np.linspace(0,1,Nx1+1); t_arr1 = np.linspace(0,t,Nt1+1)
x2 = np.linspace(0,1,Nx2+1); t_arr2 = np.linspace(0,t,Nt2+1)
x3 = np.linspace(0,1,Nx3+1); t_arr3 = np.linspace(0,t,Nt3+1)



infile_FE = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx1) + "-Nt-" + str(Nt1) +  "-FE.txt", "r")
infile_BE = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx1) + "-Nt-" + str(Nt1) +  "-BE.txt", "r")
infile_CN = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx1) + "-Nt-" + str(Nt1) +  "-CN.txt", "r")

infile_FE2 = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx2) + "-Nt-" + str(Nt2) +  "-FE.txt", "r")
infile_BE2 = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx2) + "-Nt-" + str(Nt2) +  "-BE.txt", "r")
infile_CN2 = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx2) + "-Nt-" + str(Nt2) +  "-CN.txt", "r")

infile_FE3 = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx3) + "-Nt-" + str(Nt3) +  "-FE.txt", "r")
infile_BE3 = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx3) + "-Nt-" + str(Nt3) +  "-BE.txt", "r")
infile_CN3 = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx3) + "-Nt-" + str(Nt3) +  "-CN.txt", "r")

infile_FE.readline(); infile_BE.readline(); infile_CN.readline()
infile_FE2.readline(); infile_BE2.readline(); infile_CN2.readline()
infile_FE3.readline(); infile_BE3.readline(); infile_CN3.readline()

u_FE = np.loadtxt(infile_FE); u_BE = np.loadtxt(infile_BE); u_CN = np.loadtxt(infile_CN)
u_FE2 = np.loadtxt(infile_FE2); u_BE2 = np.loadtxt(infile_BE2); u_CN2 = np.loadtxt(infile_CN2)
u_FE3 = np.loadtxt(infile_FE3); u_BE3 = np.loadtxt(infile_BE3); u_CN3 = np.loadtxt(infile_CN3)


def u(x):
    sum = 0
    T = 1
    N = int(1e5)
    for n in range(1,N):
        sum += (((-1)**(n))/n)*np.sin(n*np.pi*x)*np.exp(-T*(n*np.pi)**2)
    uxt = x + (2/np.pi)*sum
    return uxt

uxt_an1 = u(x1)
uxt_an2 = u(x2)
uxt_an3 = u(x3)

# get absolute error by using supremum norm (for simplicity)
eps1FE = abs(np.max(uxt_an1 - u_FE[-1,:]))
eps1BE = abs(np.max(uxt_an1 - u_BE[-1,:]))
eps1CN = abs(np.max(uxt_an1 - u_CN[-1,:]))

eps2FE = abs(np.max(uxt_an2 - u_FE2[-1,:]))
eps2BE = abs(np.max(uxt_an2 - u_BE2[-1,:]))
eps2CN = abs(np.max(uxt_an2 - u_CN2[-1,:]))

eps3FE = abs(np.max(uxt_an3 - u_FE3[-1,:]))
eps3BE = abs(np.max(uxt_an3 - u_BE3[-1,:]))
eps3CN = abs(np.max(uxt_an3 - u_CN3[-1,:]))

epsFE = np.array([eps1FE,eps2FE,eps3FE])
epsBE = np.array([eps1BE,eps2BE,eps3BE])
epsCN = np.array([eps1CN,eps2CN,eps3CN])

# plotting errors
plt.figure()
plt.plot(dx,epsFE,'-x', label = "FE")
plt.plot(dx,epsBE,'-x', label = "BE")
plt.plot(dx,epsCN,'-x', label = "CN")
plt.xlabel('dx',fontsize = 13)
plt.ylabel(r'supremum norm of error $|\varepsilon|$',fontsize = 13)
plt.legend()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig("./results/figures/analysis/error-analysis.png")
plt.show()
