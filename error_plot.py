filename = "error.txt"
text = open(filename).read().split("\n")
import numpy as np
import matplotlib.pyplot as plt
n = len(text)
X = np.zeros(n)
Y = np.zeros(n)
G = np.zeros(n)

for i in range(n):
    X[i] = text[i].split()[0]
    Y[i] = text[i].split()[1]
    G[i] = text[i].split()[2]
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
plt.xlim([0, 2.0])
plt.ylim([-1.0, 1.0])
cntr = ax.tricontour(X, Y, G, levels=[-2.5, 2.5, 7.5])
ax.clabel(cntr)
ax.grid()
ax.set_title("NISL error plot")
plt.show()

filename1 = "error_equator.txt"
text1 = open(filename1).read().split("\n")
n = len(text1)
X = np.zeros(n)
Y = np.zeros(n)
for i in range(n):
    X[i] = text1[i].split()[0]
    Y[i] = text1[i].split()[1]
plt.plot(X, Y)
plt.show()