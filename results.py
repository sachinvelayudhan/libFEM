import numpy as np
from numpy import loadtxt
import matplotlib.pyplot as plt

f1 = loadtxt('node_file.txt')
f2 = loadtxt('solution.txt')

#plt.plot(f1[:,0], f1[:,1], linestyle='none',marker='o', markerfacecolor='blue', markersize=12)

nodex = 3;
nodey = 3;
it = 0;
y = np.zeros((3,1))
u = np.zeros((3,1))

for i in range(0,nodey):
    y[i] = f1[it,1]
    u[i] = f2[it]
    it = 3*i + 3
    
# plotting the points 
plt.plot(u, y, color='green', linewidth = 3, marker='o', markerfacecolor='blue', markersize = 8)

# naming the x axis
plt.xlabel('Displacement, u')
# naming the y axis
plt.ylabel('y-coordinate')
  
# Title
plt.title('Variation of displacement along the y-direction')