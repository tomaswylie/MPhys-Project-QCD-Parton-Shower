import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

x = np.random.randn(10000)
y = np.random.randn(10000)

# This function creates a contour plot based on the density of points x,y
kde_plot = sns.kdeplot(x = x, y = y, fill = True, cmap = "viridis", levels = 10, thresh = 0.05)
cbar = plt.colorbar(kde_plot.collections[-1], label = 'Density')

plt.scatter(x, y, s = 1, alpha = 0.3)

plt.xlabel('x')
plt.ylabel('y')
plt.title('KDE contour plot of point density')
plt.savefig("kdetest.png")
plt.show()

