# MPHYS-Project---QCD-Parton-Showers
A shared repository to house the code needed for the MPHYS project on QCD Parton Showers.

# The following files were provided to us, and we have not made any changes to them
.matrix.py - Prewritten matrix element generation
.shower.py - Prewritten shower program
durham.py - Generates data from yij analysis
eshapes.py - Generates data from event shape analysis
histogram.py - Handles binning of data
matrix.py - Handles matrix element generation
particle.py - Defines some parameters of particles
plots.conf - Defines parameters for event shape plots
Sherpa.yoda - yij data from Sherpa
vector.py - Defines vector operations

# The following files were provided to us, and we have been making continuous changes to
shower.py - Modified shower class
qcd.py - Class for strong running coupling, with CMW scheme added

# The following files were created by us
kdetest.py - Testing the kernel density equation contour generation
kernels.py - Visualise relative contribution of splitting kernels
lund.py - Generates Lund planes
runshower.py - Runs shower and performs event shape analysis

# Access Rivet
docker pull hepstore/rivet:latest

# Running rivet with docker
docker run -it --rm hepstore/rivet:latest

# Command line for Sam
docker run -it --rm -v /c/Users/samue/OneDrive/Documents/Python_Scripts/MPHYS-Project---QCD-Parton-Showers:/app -w /app hepstore/rivet:latest /bin/bash

# Setting up alias to run rivet for Tom
alias rivet='docker run -it --rm -u `id -u tomas`:`id -g` -v "/Volumes/Macintosh HD/Users/tomas/python_scripts/MPHYS-Project---QCD-Parton-Showers/":/app -w /app hepstore/rivet:latest'

# Shortcut to run rivet - enter this in the command line
rivet

# Setting up shortcut to run rivet-mkhtml for Tom
alias rivet-mkhtml='docker run -it --rm -u `id -u tomas`:`id -g` -v "/Volumes/Macintosh HD/Users/tomas/python_scripts/MPHYS-Project---QCD-Parton-Showers/":/app -w /app hepstore/rivet:latest rivet-mkhtml'

# Shortcut make plots using rivet-mkhtml
rivet-mkhtml -c plots.conf myshower.yoda
exit