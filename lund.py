import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines
import random as r

from shower import *
from particle import CheckEvent
from qcd import AlphaS
from matrix import eetojj

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Set the energy scale according to the Z mass
ecom = 91.1876

# Set up classes for running coupling, hard scattering and shower
alphas = AlphaS(ecom,0.118,cmw=True)
hardxs = eetojj(alphas)
shower = Shower(alphas,t0=1.)

def Magnitude3(mom):
    '''
    This function calculates the magnitude of the spacelike components of a given 4-momentum

    Parameters
    ----------
    mom (arr) : The 4-momentum of the parton

    Returns
    -------
    (float) : The magnitude of the particle's 3-momentum
    '''

    return m.sqrt(mom[1]**2 + mom[2]**2 + mom[3]**2)
    
def EmissionAngle(mom1, mom2):
    '''
    This function determines the angle of emission between the emitting parton and the emitted parton

    Parameters
    ----------
    mom1 (arr) : The 4-momentum of the emitting parton
    mom2 (arr) : The 4-momentum of the emitted parton

    Returns
    -------
    (float) : The emission angle between the emitting parton and the emitted parton
    '''

    dotprod = (mom1[1]*mom2[1])+(mom1[2]*mom2[2])+(mom1[3]*mom2[3])
    return m.acos(dotprod / (Magnitude3(mom1) * Magnitude3(mom2)))

def constantkt(kt_values):
    '''
    This function plots a series of lines of constant transverse momentum.

    Parameters
    ----------
    kt_values (arr) : The values of transverse momentum
    '''

    # Generate angles distributed geometrically across the curve
    theta = np.geomspace(0.001, np.pi, 250)

    for kt in kt_values:
        # Define energy based on relation to transverse momentum kt
        E = kt / np.sin(theta)

        # Ensure plotted lines stay within allowed phase space
        mask = 2 * E / ecom <= 1 - theta
        theta_f = theta[mask]
        E_f = E[mask]

        plt.plot(theta_f, 2 * E_f / ecom, color="k", linestyle="-.", alpha=0.5,
                     label=rf"$k_\perp$={kt}")

    # Add labels to the lines of constant kt
    labelLines(plt.gca().get_lines(), align=True, fontsize=10, color="grey", zorder = 2.5)
    plt.legend([])

    return

def emissions1(evno):
    '''
    This function generates a Lund plane for a series of events with one emission, parameterised by the
    normalised angle of emission and normalised energy of the emitted parton.

    Parameters
    ----------
    evno (int) : The number of events to generate
    '''

    fig, ax = plt.subplots()

    # Set the seed for random number generation
    r.seed(123456)

    em_angle = np.zeros(evno)
    em_energy = np.zeros(evno)

    for ev in range(evno):
        # Generate hard scattering and define maximum value of kt
        event, weight = hardxs.GenerateLOPoint()
        t = (event[0].mom+event[1].mom).M2()

        # Ensure that only one event is generated
        while len(event) < 5:
            shower.c = 1
            shower.t = t
            moms = shower.GeneratePoint(event)
        
        # Ensure momentum conservation
        if not CheckEvent(event):
            print("Something went wrong:")
            for p in event:
                print(p)
        
        # Define parameters for the emission
        pi = moms[0][0]
        pj = moms[0][1]
        em_angle[ev] = EmissionAngle((pi+pj), pj) / np.pi
        em_energy[ev] = 2*pj[0] / ecom

        # Check generated point lies in allowed phase space
        if em_energy[ev] >= 1 - em_angle[ev]:
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan
        elif em_energy[ev] <= 0 or em_angle[ev] <= 0:
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan

        # Increment progress bar as a visual aid
        progress_bar(ev + 1, evno)

    # Define offset close to zero to avoid log errors
    offset = 0.005

    # Plot phase space boundaries
    x1, y1 = np.linspace(offset, 1, 100), np.full(100, offset)
    x2, y2 = np.linspace(offset, 1, 100), 1 - np.linspace(offset, 1, 100)
    x3, y3 = np.full(100, offset), np.linspace(offset, 1, 100)
    ax.plot(x1, y1, color="k", linestyle="--", label="Phase space boundary")
    ax.plot(x2, y2, color="k", linestyle="--")
    ax.plot(x3, y3, color="k", linestyle="--")

    # Plot emissions
    ax.scatter(em_angle, em_energy, s=1, label = "Emissions")

    ax.set_xlabel(r"$\theta/\pi$")
    ax.set_ylabel(r"$2E/Q$")
    ax.set_title("Lund jet plane for first emission")

    # Create inset plot for log scale
    inset = ax.inset_axes([0.6, 0.45, 0.38, 0.38])

    # Plot phase space boundaries for inset
    inset.plot(x1, y1, color="k", linestyle="--", label="Phase space boundary")
    inset.plot(x2, y2, color="k", linestyle="--")
    inset.plot(x3, y3, color="k", linestyle="--")

    # Plot emission for inset
    inset.scatter(em_angle, em_energy, s=1, label = "Emissions")

    inset.set_xlabel(r"$\theta/\pi$")
    inset.set_ylabel(r"$2E/Q$")
    inset.set_xscale("log")
    inset.set_yscale("log")

    plt.legend(loc = "upper right")
    plt.savefig("lund1em_10000_CMW.pdf")
    plt.show()

    return

def emissionsN(evno, log=False):
    '''
    This function generates a Lund plane for a specified number of events, allowed to evolve until all
    emissions are vetoed. The plane is parameterised in terms of the normalised angle of emission and 
    normalised energy of the emitted parton.

    Parameters
    ----------
    evno (int) : The number of events to generate
    log (bool) : Whether the graph is to be displayed in a log scale or not
    '''

    plt.figure()

    # Set the seed for random number generation
    r.seed(123456789)
    
    kt_list = []
    for ev in range(0, evno):
        em_angle = []
        em_energy = []
        moms_list = []

        # Generate hard scattering event and define highest transverse momentum 
        event, weight = hardxs.GenerateLOPoint()
        t = (event[0].mom+event[1].mom).M2()
        
        # Run the shower and return all intermediate momenta and t values
        gen = shower.Run_Return(event, t)

        # Extract momenta and t values
        for moms, tlist in gen:
            moms_list.append(moms)
            kt_list.append(m.sqrt(tlist[0]))
        kt_list = [x for x in kt_list if x != 0]

        # Check momentum conservation
        if not CheckEvent(event):
            print("Something went wrong:")
            for p in event:
                print(p)

        # Check for zero emissions
        if moms_list == [] or moms_list[0] == None: 
            print("Event {0} terminated due to invalid phase space conditions".format(ev))
            continue

        for em in range(0, len(moms_list)-1):
            # Define parameters for emission
            pi = moms_list[em][0]
            pj = moms_list[em][1]
            em_angle.append(EmissionAngle(pi+pj, pj) / np.pi)
            em_energy.append(2*pj[0] / (ecom))

            # Check if emissions lie in allowed phase space
            if em_angle[em] <= 0 or em_energy[em] <= 0 or em_energy[em] >= 1 - em_angle[em]:
                print("Event {0} terminated due to invalid phase space conditions".format(ev))
                em_angle = em_angle[0:em]
                em_energy = em_energy[0:em]
                break
        if em_angle != [] and em_energy != []:
            plt.scatter(em_angle[0], em_energy[0], s = 10, label = "Event")
            plt.plot(em_angle, em_energy)

    # Set log scale
    if log:
        plt.xscale("log")
        plt.yscale("log")
    else:
        plt.ylim((-0.1,1.1))

    # Plot lines of constant kt
    kt_values = [1, 3, 5, 7, 10]
    constantkt(kt_values)

    # Define offset close to zero to avoid log errors
    offset = 0.005

    # Plot phase space boundaries
    x1, y1 = np.linspace(offset, 1, 100), np.full(100, offset)
    x2, y2 = np.linspace(offset, 1, 100), 1 - np.linspace(offset, 1, 100)
    x3, y3 = np.full(100, offset), np.linspace(offset, 1, 100)
    plt.plot(x1, y1, color="k", linestyle="--", label="Phase space boundary")
    plt.plot(x2, y2, color="k", linestyle="--")
    plt.plot(x3, y3, color="k", linestyle="--")

    plt.xlabel(r"$\theta/\pi$")
    plt.ylabel(r"$2E/Q$")
    plt.title("Lund jet plane for N emissions")
    plt.savefig("lundemNlog_ev1.pdf")
    plt.legend(loc = "best")
    plt.show()

    return

def emissions_temp(evno):
    '''
    This function generates a Lund plane for a series of events with one emission, parameterised by
    theta/pi and 2E/Q sin(theta/pi).

    Parameters
    ----------
    evno (int) : The number of events to generate
    '''

    fig, ax = plt.subplots()

    # Set the seed for random generation
    r.seed(123456)

    em_angle = np.zeros(evno)
    em_energy = np.zeros(evno)

    for ev in range(evno):
        # Generate hard scattering event and define highest transverse momentum
        event, weight = hardxs.GenerateLOPoint()
        t = (event[0].mom+event[1].mom).M2()

        # Ensure one emission is generated
        while len(event) < 5:
            shower.c = 1
            shower.t = t
            moms = shower.GeneratePoint(event)

        # Ensure momentum conservation
        if not CheckEvent(event):
            print("Something went wrong:")
            for p in event:
                print(p)

        # Define parameters for emission
        pi = moms[0][0]
        pj = moms[0][1]
        em_angle[ev] = EmissionAngle((pi+pj), pj)
        em_energy[ev] = np.sin(em_angle[ev])*pj[0]

        # Check if emission lies in allowed phase space
        if em_energy[ev] >= ecom/2*np.sin(em_angle[ev])*(1-em_angle[ev]/np.pi):
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan
        elif em_energy[ev] <= 0 or em_angle[ev] <= 0:
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan

        # Increment progress bar as a visual aid
        progress_bar(ev + 1, evno)

    # Define offset close to zero to avoid log errors
    offset = 0.005

    # Plot phase space boundaries
    angles = np.linspace(offset, np.pi, 100)
    x1, y1 = angles, np.full(100, offset)
    x2, y2 = angles, ecom/2*np.sin(angles)*(1-angles/np.pi)
    ax.plot(x1, y1, color="k", linestyle="--", label="Phase space boundary")
    ax.plot(x2, y2, color="k", linestyle="--")

    # Plot emissions
    ax.scatter(em_angle, em_energy, s=1, label = "Emissions")

    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$E\sin\theta$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Lund jet plane for first emission")

    plt.legend(loc = "upper right")
    plt.savefig("lund1em_Esintheta.png")
    plt.show()
    
    return

def constantkt_temp(kt_values):
    '''
    This function plots a series of lines of constant transverse momentum for the parameterisation of 
    the Lund plane in terms of theta/pi and 2E/Q sin(theta/pi).

    Parameters
    ----------
    kt_values (arr) : The values of transverse momentum
    '''

    # Generate angles distributed geometrically
    theta = np.geomspace(0.001, np.pi, 250)

    for kt in kt_values:
        # Define energy from kt relation
        E = kt / np.sin(theta)

        # Filter out valid phase space points
        mask = 2*E/ecom*np.sin(theta/np.pi) < (1-theta/np.pi)*(theta/np.pi*np.cos(theta/np.pi)+np.sin(theta/np.pi))
        theta_f = theta[mask]
        E_f = E[mask]

        plt.plot(theta_f/np.pi, 2*E_f/ecom*np.sin(theta_f/np.pi), color="k", linestyle="-.", alpha=0.5,
                     label=rf"$k_\perp$={kt}")

    # Label the lines of constant kt
    labelLines(plt.gca().get_lines(), align=True, fontsize=10, color="grey", zorder = 2.5)
    plt.legend([])

    return

def emissionsN_temp(evno, log=False):
    '''
    This function generates a Lund plane for a specified number of events, allowed to evolve until all
    emissions are vetoed. The plane is parameterised in terms of theta/pi and 2E/Q sin(theta/pi).

    Parameters
    ----------
    evno (int) : The number of events to generate
    log (bool) : Whether the graph is to be displayed in a log scale or not
    '''

    plt.figure()

    # Set the seed for random number generation
    r.seed(253)
    
    kt_list = []
    for ev in range(0, evno):
        em_angle = []
        em_energy = []
        moms_list = []

        # Generate hard scattering event and define highest value of transverse momentum
        event, weight = hardxs.GenerateLOPoint()
        t = (event[0].mom+event[1].mom).M2()

        # Return all intermediate momenta and t values
        gen = shower.Run_Return(event, t)

        # Extract momenta and t values
        for moms, tlist in gen:
            moms_list.append(moms)
            kt_list.append(m.sqrt(tlist[0]))
        kt_list = [x for x in kt_list if x != 0]

        # Ensure momentum conservation
        if not CheckEvent(event):
            print("Something went wrong:")
            for p in event:
                print(p)
        
        # Check for no emissions
        if moms_list == []: 
            print("Event {0} terminated due to invalid phase space conditions".format(ev))
            continue

        for em in range(0, len(moms_list)-1):
            # Define parameters for emission
            pi = moms_list[em][0]
            pj = moms_list[em][1]
            em_angle.append(EmissionAngle((pi+pj), pj)/np.pi)
            em_energy.append(2*np.sin(em_angle[em])*pj[0]/ecom)

            # Check emission lies in allowed phase space
            if em_angle[em] <= 0 or em_energy[em] <= 0 or em_energy[em] >= (1-em_angle[em])*(em_angle[em]*np.cos(em_angle[em])+np.sin(em_angle[em])):
                print("Event {0} terminated due to invalid phase space conditions".format(ev))
                em_angle = em_angle[0:em]
                em_energy = em_energy[0:em]
                break
        if em_angle != [] and em_energy != []:
            plt.scatter(em_angle[0], em_energy[0], s = 10, label = "Event")
            plt.plot(em_angle, em_energy)

    # Plot lines of constant kt
    kt_values = [1, 3, 5, 7, 10]
    constantkt_temp(kt_values)

    # Define offset close to zero to avoid log errors
    offset = 0.005
    
    # Plot phase space boundary
    angles = np.linspace(offset, 1, 100)
    x1, y1 = angles, (1-angles)*(angles*np.cos(angles)+np.sin(angles))

    plt.plot(x1, y1, color="k", linestyle="--")

    # Set log scale
    if log:
        plt.xscale("log")
        plt.yscale("log")

    plt.xlabel(r"$\theta/\pi$")
    plt.ylabel(r"$2E/Q\cdot\sin\left(\theta/\pi\right)$")
    plt.title("Lund jet plane for N emissions")
    plt.legend(loc = "best", ncols = 3)
    plt.savefig("lundNem_Esintheta2.png")
    plt.show()

    return

def emissions1_3d_histo(evno, bins, lower_lim =0, upper_lim=1, log=False):
    
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')

    r.seed(123456)
    em_angle = np.zeros(evno)
    em_energy = np.zeros(evno)

    for ev in range(evno):
        event, weight = hardxs.GenerateLOPoint()
        t = (event[0].mom+event[1].mom).M2()
        while len(event) < 5:
            shower.c = 1
            shower.t = t
            moms = shower.GeneratePoint(event)
        if not CheckEvent(event):
            print("Something went wrong:")
            for p in event:
                print(p)
        pi = moms[0][0]
        pj = moms[0][1]
        em_angle[ev] = EmissionAngle((pi+pj), pj) / np.pi
        em_energy[ev] = 2*pj[0] / ecom
        if em_energy[ev] >= 1 - em_angle[ev]:
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan
        elif em_energy[ev] <= 0 or em_angle[ev] <= 0:
            em_angle[ev] = np.nan
            em_energy[ev] = np.nan

        progress_bar(ev + 1, evno)
        
    em_angle = em_angle[~np.isnan(em_angle)]
    em_energy = em_energy[~np.isnan(em_energy)]

    # Create 3D histogram
    # Define the number of bins and the ranges for each axis
    x_bins = y_bins = np.linspace(lower_lim, upper_lim, bins)  # sets the bins

    if log:
        hist, xedges, yedges = np.histogram2d(np.log10(em_angle), np.log10(em_energy), bins=[x_bins, y_bins])
        ax.set_xlabel(r'Logged angle $\log_{10}{\theta}$')
        ax.set_ylabel(r'Logged energy $\log_{10}{E}$')  
    else:
        hist, xedges, yedges = np.histogram2d(em_angle, em_energy, bins=[x_bins, y_bins])
        ax.set_xlabel(r'Angle $\theta$ (in \pi units)')
        ax.set_ylabel(r'Energy $E$ (in GeV)')

    # Set labels
    ax.set_zlabel('Frequency')

    # Compute the bar heights
    dx = dy = (upper_lim-lower_lim)/bins
    dz = hist.flatten()

    offset = dx/2

    # Construct the grid for plotting
    xpos, ypos = np.meshgrid(xedges[:-1] + offset, yedges[:-1] + offset)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)

    #save_histogram_data(hist,xedges,yedges)
    
    # Normalize the histogram counts (dz) for colormap
    dz_normalized = (dz - dz.min()) / (dz.max() - dz.min())  # Normalize to [0, 1]

    # Apply a colormap (Red for high, Green for low)
    colors = plt.cm.YlOrRd(dz_normalized)
    """
    # Adjust alpha (transparency) for yellow values
    for i, color in enumerate(colors):
        if 0.3 < dz_normalized[i] < 0.7:  # Range for "yellow"
            colors[i][-1] = 1.0  # Set alpha to full opacity"""

    # Find bin centres
    x_centres = (xedges[:-1] + xedges[1:]) / 2
    y_centres = (yedges[:-1] + yedges[1:]) / 2

    # Plot surface
    x_arr, y_arr = np.meshgrid(x_centres, y_centres)
    ax.plot_surface(x_arr, y_arr, hist, facecolors=plt.cm.winter(hist / hist.max()), shade=False, 
                    edgecolor = "black")

    # Plot the 3D histogram
    #ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')

    plt.savefig("histsurfcmw.png")
    plt.show()
    
    return

def save_histogram_data(hist, xedges, yedges, filename="histogram_data.txt"):
    """
    Saves the histogram data to a text file.
    """
    with open(filename, 'w') as f:
        # Write bin edges for reference
        f.write("# X bin edges:\n")
        f.write(" ".join(map(str, xedges)) + "\n")
        f.write("# Y bin edges:\n")
        f.write(" ".join(map(str, yedges)) + "\n")
        
        # Write histogram counts
        f.write("# Histogram counts:\n")
        for row in hist:
            f.write(" ".join(map(str, row)) + "\n")
            
def read_histogram_data(filename="histogram_data.txt"):
    """
    Reads histogram data from a text file saved by the `save_histogram_data` function.
    Returns the histogram counts (as a 2D numpy array), x bin edges, and y bin edges.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Parse the file for bin edges and histogram counts
    xedges = []
    yedges = []
    hist = []

    is_reading_hist = False  # Flag to determine if we are reading the histogram counts

    for line in lines:
        if line.startswith("# X bin edges:"):
            xedges = np.fromstring(lines[lines.index(line) + 1], sep=" ")
        elif line.startswith("# Y bin edges:"):
            yedges = np.fromstring(lines[lines.index(line) + 1], sep=" ")
        elif line.startswith("# Histogram counts:"):
            is_reading_hist = True
        elif is_reading_hist:
            # Read histogram rows
            hist.append(list(map(float, line.split())))

    # Convert data to numpy arrays
    hist = np.array(hist)
    return hist, xedges, yedges


emissionsN(1, True)
