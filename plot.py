filename = "log.txt"
text = open(filename).read().split("\n")
import numpy as np
import matplotlib.pyplot as plt
n = len(text)
X = np.zeros(n)
Y = np.zeros(n)
G = np.zeros(n)

id = 0
ma = 0
for i in range(n-1):
    X[i] = text[i].split()[0]
    Y[i] = text[i].split()[1]
    G[i] = text[i].split()[2]
print(max(G))
fig, ax = plt.subplots(figsize=(10,5))
cntr = ax.tricontourf(X, Y, G, 14, cmap="jet")
fig.colorbar(cntr, ax=ax)
ax.grid()
plt.show()