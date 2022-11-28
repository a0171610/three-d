text = open("animation.txt").read().split("\n")
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from PIL import Image
import os

step = 20
n = len(text)
n = int(n / step)
X = np.zeros((step, n))
Y = np.zeros((step, n))
G = np.zeros((step, n))

for i in range(step):
    for j in range(n):
        X[i, j] = text[i * n + j].split()[0]
        Y[i, j] = text[i * n + j].split()[1]
        G[i, j] = text[i * n + j].split()[2]

ims = []
cdir = os.getcwd()

fig, ax = plt.subplots()

for i in range(step):
    im = ax.tricontourf(X[i, :], Y[i, :], G[i, :], 14, cmap="jet", levels=np.linspace(-0.05,1.05, 23))
    if i == 0:
        fig.colorbar(im)
    add_arts = im.collections
    filename = os.path.join(cdir, "No{}.png".format(i))
    plt.savefig(filename)
    im = Image.open(filename)
    ims.append(im)

ims[0].save(os.path.join(cdir, "NISL_animation.gif"), save_all=True, append_images=ims[1:], loop=0, duration=15)

for i in range(step):
    filename = os.path.join(cdir, "No{}.png".format(i))
    os.remove(filename)
