import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter

dim = int(input("Want to plot 1D or 2D solution (1/2)?: "))

if dim == 1:

    stab_crit = str(input("Plot with stability criterion not fullfilled (yes/no)?: "))
    if stab_crit == 'yes':
        Nx = int(input("Set Nx (int):")) # 100
        Nt =  int(input("Set Nt (int):"))# 3000
        t2 = 0.1
        dx = 1/Nx

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
        t_array = np.linspace(0,t2,Nt+1)
        t1 = t_array[int(Nt/500)]
        print(t1)

        def u(x):
            sum = 0
            N = int(1e5)
            for n in range(1,N):
                sum += (((-1)**(n))/n)*np.sin(n*np.pi*x)*np.exp(-t1*(n*np.pi)**2)
            uxt = x + (2/np.pi)*sum
            return uxt

        uxt_an = u(x)


        plt.plot(x,u_FE[int(Nt/500),:], label = 'Forward Euler solution at $t_1$')
        plt.plot(x,u_BE[int(Nt/500),:], label = 'Implicit Backward solution at $t_1$')
        plt.plot(x,u_CN[int(Nt/500),:], label = 'Crank-Nicholson solution at $t_1$')
        plt.plot(x,uxt_an, label = 'Analytical solution $u(x,t)$, at $t_1$')
        plt.xlabel('x',fontsize = 13)
        plt.ylabel('$u(x,t)$',fontsize = 13)
        plt.title('Forward Euler unstable at $t_1 = 0.0002$ and $\Delta x$ = 0.01',fontsize = 14)
        plt.legend()
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.savefig("./results/figures/1d-sol/stab_crit_fail_t={:2f}dx={:2f}.png".format(t1,dx))

        plt.show()



    if stab_crit == 'no':
        Nx = int(input("Set Nx (int):")) # 100
        Nt =  int(input("Set Nt (int):"))# 3000
        t2 = 1.0
        dx = 1/Nx

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
        t_array = np.linspace(0,t2,Nt+1)
        t1 = t_array[int(Nt/5)]

        def u(x):
            sum = 0
            T = 10
            N = int(1e5)
            for n in range(1,N):
                sum += (((-1)**(n))/n)*np.sin(n*np.pi*x)*np.exp(-T*(n*np.pi)**2)
            uxt = x + (2/np.pi)*sum
            return uxt

        uxt_an = u(x)

        plt.figure()
        plt.plot(x,u_FE[int(Nt/5),:], label = 'Forward Euler solution at $t_1$')
        plt.plot(x,u_BE[int(Nt/5),:], label = 'Implicit Backward solution at $t_1$')
        plt.plot(x,u_CN[int(Nt/5),:], label = 'Crank-Nicholson solution at $t_1$')
        plt.plot(x,uxt_an, label = 'Analytical solution $u(x,t)$, at t = 10')
        plt.xlabel('x',fontsize = 13)
        plt.ylabel('$u(x,t)$',fontsize = 13)
        plt.title('Numerical solutions at time $t_1 = 0.2$ and $\Delta x$ = 0.01',fontsize = 14)
        plt.legend()
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.savefig("./results/figures/1d-sol/FE_BE_CN_t1={:2f}dx={:2f}.png".format(t1,dx))


        plt.figure()
        plt.plot(x,u_FE[-1,:], label = 'Forward Euler solution at $t_2$')
        plt.plot(x,u_BE[-1,:], label = 'Implicit Backward solution at $t_2$')
        plt.plot(x,u_CN[-1,:], label = 'Crank-Nicholson solution at $t_2$')
        plt.plot(x,uxt_an, label = 'Analytical solution $u(x,t)$, at t = 10')
        plt.xlabel('x',fontsize = 13)
        plt.ylabel('$u(x,t)$',fontsize = 13)
        plt.title('Numerical solutions at time $t_2 = 1.0$ and $\Delta x$ = 0.01',fontsize = 14)
        plt.legend()
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.savefig("./results/figures/1d-sol/FE_BE_CN_t2={:2f}dx={:2f}.png".format(t2,dx))

        plt.show()



if dim == 2:
    # make a 2D plot
    Nx = int(input("Set Nx (int):")) # 10
    Ny = int(input("Set Ny (int):")) # 10
    Nt = int(input("Set Nt (int):")) # 6000


    infile_2DBE = open("./results/2D-solutions/" + "2Dsol-Nx-" + str(Nx) + "-Ny-" + str(Ny) + "-Nt-" + str(Nt) +  "-BE.txt", "r")

    # read from file
    line = infile_2DBE.readline()
    T = line.split()[0]

    # get vector
    u2D = np.loadtxt(infile_2DBE)

    # Create mesh for surface plotting
    x = np.linspace(0,1,Nx+1)
    y = np.linspace(0,1,Ny+1)
    X,Y = np.meshgrid(x,y)

    # initialize plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot the surface.
    surf = ax.plot_surface(X, Y, u2D,cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(0, 0.8)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    ## chose colour of axis
    ax.tick_params(axis='x',colors='grey')
    ax.tick_params(axis='y',colors='grey')
    ax.tick_params(axis='z',colors='grey')

    ax.w_xaxis.line.set_color("grey")
    ax.w_yaxis.line.set_color("grey")
    ax.w_zaxis.line.set_color("grey")

    # formatting
    plt.title('{:}'.format(T),fontsize = 14)
    ax.set_xlabel('x',fontsize = 13)
    ax.set_ylabel('y',fontsize = 13)
    ax.set_zlabel("u(x,y)", fontsize = 13)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.tight_layout()

    plt.savefig("./results/figures/2d-sol/2d-Nx-{:d}-Ny-{:d}-Nt-{:d}.png".format(Nx,Ny,Nt))
    plt.show()
