import numpy as np
import matplotlib.pyplot as plt
from qcd import CF, CA, TR

# Generate values for kinematics, limited to [0.1,0.9] to avoid divergences
y = np.linspace(0.1, 0.9, 200)
z = np.linspace(0.1, 0.9, 200)

def kernels2d():
    '''
    This function calculates the Catani-Seymour dipole splitting functions for y=1, and plots them
    on the same figure to show their relative contributions.
    '''

    # Define splitting functions for y=1
    Pqq = CF*(2/(1-z)-(1+z))
    Pgq = TR/2*(1-2*z*(1-z))
    Pgg = CA/2*(2/(1-z)-2+z*(1-z))

    plt.plot(z, Pqq, label = "q->qg")
    plt.plot(z, Pgq, label = "g->qq")
    plt.plot(z, Pgg, label = "g->gg")
    plt.xlabel("z")
    plt.ylabel("P(z)")
    plt.yscale("log")
    plt.title("Relative contribution of splitting kernels")
    plt.legend()
    plt.savefig("kernels.png")
    plt.show()

def kernels3d():
    '''
    This function calculates the Catani-Seymour dipole splitting functions, and plots them
    on the same figure to show their relative contributions.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    # Create a grid of (y,z) values
    y_arr, z_arr = np.meshgrid(y, z)

    # Define splitting functions
    Pqq = CF*(2/(1-z_arr*(1-y_arr))-(1+z_arr))
    Pgq = TR/2*(1-2*z_arr*(1-z_arr))
    Pgg = CA/2*(2/(1-z_arr*(1-y_arr))-2+z_arr*(1-z_arr))

    ax.plot_surface(y_arr, z_arr, Pqq, label = "q->qg")
    ax.plot_surface(y_arr, z_arr, Pgq, label = "g->qq")
    ax.plot_surface(y_arr, z_arr, Pgg, label = "g->gg")

    ax.set_xlabel("y")
    ax.set_ylabel("z")
    ax.set_zlabel("P(y,z)")
    ax.set_zscale("log")
    ax.set_title("Relative contribution of splitting kernels")
    plt.legend()
    plt.savefig("kernels3d.png")
    plt.show()

kernels2d()