from shower import *

import sys
from matrix import eetojj
from durham import Analysis
from eshapes import EAnalysis

# Set up classes for strong coupling, hard scatter, shower and analysis
alphas = AlphaS(91.1876,0.118,cmw=False)
hardxs = eetojj(alphas)
shower = Shower(alphas,t0=1.)
jetrat = Analysis()
eshape = EAnalysis()

# Set the seed for random number generation
r.seed(123456)

# Define number of events
eventno = 100000
for i in range(eventno):
    # Generate hard scattering event and define highest value of ordering variables
    event, weight = hardxs.GenerateLOPoint()
    t = (event[0].mom+event[1].mom).M2()

    # Generate full evolution of the shower
    shower.Run(event,t)

    # Check momentum conservation
    if not CheckEvent(event):
        print("Something went wrong:")
        for p in event:
            print(p)

    # Increment progress bar for visual aid
    progress_bar(i, eventno)

    # Perform event shape analysis
    jetrat.Analyze(event, weight)
    eshape.Analyze(event, weight)

# Format data files for plotting
jetrat.Finalize("myshower")
eshape.Finalize("myshower")
print("")