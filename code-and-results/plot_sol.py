import numpy as np
import matplotlib.pyplot as plt

dim = int(input("Want to plot in 1D or 2D (1/2)?: "))

if dim == 1:
    Nx = 100
    Nt = 3000
    T = 0.1
    infile_FE = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx) + "-Nt-" + str(Nt) +  "-FE.txt", "r")
    infile_BE = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx) + "-Nt-" + str(Nt) +  "-BE.txt", "r")
    infile_CN = open("./results/1D-solutions/" + "1Dsol-Nx-" + str(Nx) + "-Nt-" + str(Nt) +  "-CN.txt", "r")

    infile_FE.readline()
    infile_BE.readline()
    infile_CN.readline()

    u_FE = np.loadtxt(infile_FE)
    u_BE = np.loadtxt(infile_BE)
    u_CN = np.loadtxt(infile_CN)

    x = np.linspace(0,1,Nx+1)
    t_array = np.linspace(0,0.1,Nt+1)


    def u(x):
        sum = 0
        N = int(1e5)
        for n in range(1,N):
            sum += (((-1)**(n))/n)*np.sin(n*np.pi*x)*np.exp(-T*(n*np.pi)**2)
        uxt = x + (2/np.pi)*sum
        return uxt

    uxt_an = u(x)



    plt.figure()
    plt.plot(x,u_FE[int(Nt/5),:], label = 't1')
    plt.plot(x,u_FE[-1,:], label = 't2')
    plt.plot(x,uxt_an, label = 'analytical expression')
    plt.title('Forward Euler solution')
    plt.legend()

    plt.figure()
    plt.plot(x,u_BE[int(Nt/5),:], label = 't1')
    plt.plot(x,u_BE[-1,:], label = 't2')
    plt.plot(x,uxt_an, label = 'analytical expression')
    plt.title('Implisit backward solution')
    plt.legend()

    plt.figure()
    plt.plot(x,u_CN[int(Nt/5),:], label = 't1')
    plt.plot(x,u_CN[-1,:], label = 't2')
    plt.plot(x,uxt_an, label = 'analytical expression')
    plt.title('Crank-Nicholson solution')
    plt.legend()
    plt.show()

if dim == 2:
    #plott 2D
    d = 2
