import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

dim = int(input("Want to plot in 1D or 2D (1/2)?: "))

if dim == 1:
    Nx = int(input("Set Nx (int):")) # 100
    Nt =  int(input("Set Nt (int):"))# 3000
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
    plt.plot(x,u_FE[int(Nt/5),:], label = 'Forward Euler t1')
    plt.plot(x,u_BE[int(Nt/5),:], label = 'Implicit Backward t1')
    plt.plot(x,u_CN[int(Nt/5),:], label = 'Crank-Nicholson t1')
    plt.plot(x,uxt_an, label = 'analytical expression')
    plt.title('Method comparison for time t1')
    plt.legend()


    plt.figure()
    plt.plot(x,u_FE[-1,:], label = 'Forward Euler $t_2$')
    plt.plot(x,u_BE[-1,:], label = 'Implicit Backward t2')
    plt.plot(x,u_CN[-1,:], label = 'Crank-Nicholson t2')
    plt.plot(x,uxt_an, label = 'analytical expression')
    plt.title('Method comparison for time t2')
    plt.legend()

    plt.show()


if dim == 2:
    # make a 2D plot
    #T = 0.1
    Nx = int(input("Set Nx (int):")) # 10
    Ny = int(input("Set Ny (int):")) # 10
    Nt = int(input("Set Nt (int):")) # 6000


    infile_2DBE = open("./results/2D-solutions/" + "2Dsol-Nx-" + str(Nx) + "-Ny-" + str(Ny) + "-Nt-" + str(Nt) +  "-BE.txt", "r")

    # read from file
    line = infile_2DBE.readline()
    T = line.split()[1]

    # get vector
    u2D = np.loadtxt(infile_2DBE)

    #print(u2D)

    # Create mesh for surface plotting
    x = np.linspace(0,1,Nx+1)
    y = np.linspace(0,1,Ny+1)
    X,Y = np.meshgrid(x,y)

    # initialize plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    print("u2d shape:",np.shape(u2D))
    print("\n")
    print("xy shape:",np.shape(X))
    # Plot the surface.
    surf = ax.plot_surface(X, Y, u2D,cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

    # Customize the z axis.
    #ax.set_xlim(-1.01, 1.01)
    #ax.zaxis.set_major_locator(LinearLocator(10))
    #ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()
