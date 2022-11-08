import matplotlib.pyplot as plt
import numpy as np

filename = "./log_cbell.txt"
filename1 = "./log_ccbell.txt"

text = open(filename).read().split("\n")
text1 = open(filename1).read().split("\n")

n = len(text)
G = np.zeros(n)
G1 = np.zeros(n)
G2 = np.zeros(n)
G3 = np.zeros(n)
C = np.zeros(n)
W = np.zeros(n)
L = np.zeros(n)
chi = np.zeros(n)
psi = np.zeros(n)

a = -0.8
b = 0.9

for i in range(n):
  G[i] = text[i].split()[0]

chi_min = 0.1
chi_max = 1.0

surf_x = np.linspace(0.1, 1.0, n)
surf_y = np.zeros_like(surf_x)
surf_y[:] = surf_x[:] * surf_x[:] * a + b

for i in range(n):
  W[i] = text1[i].split()[0]
  G1[i] = text1[i].split()[1]
  G2[i] = G[i] * G[i] * a + b
  C[i] = ((65340*G[i] + 12.0*np.sqrt(12.0*(125.0*G1[i] - 52.0)**3 + 29648025*G[i]**2)) ** (1.0/3.0)) / (60.0)
  chi[i] = C[i] + (13.0/75.0 - 5.0/12.0*G1[i]) / C[i]
  psi[i] = chi[i] * chi[i] * a + b

g2_diff = 0.792
g2_max = 0.892
g2_min = 0.1

for i in range(n):
  L[i] = np.sqrt(( ((G[i] - chi[i]) / (chi_max - chi_min)) **2) + (((G1[i] - psi[i]) / (g2_diff)) ** 2))
  G3[i] = g2_max + ((g2_min - g2_max) / (chi_max - chi_min)) * (G[i] - chi_min)


l2_real = 0.0
l2_shape = 0.0
l2_over = 0.0

for i in range(n):
  if ((0.1 > G[i] and G[i] < 1.0) or (0.1 > G1[i] and G1[i] > 0.892)):
    l2_over += W[i] * L[i] / np.sum(W)
  elif (0.1 < G[i] and G[i] < 1.0 and G3[i] < G1[i] and G1[i] < G2[i]):
    l2_real += W[i] * L[i] / np.sum(W)
  else:
    l2_shape += W[i] * L[i] / np.sum(W)

print(l2_real, l2_shape, l2_over)

plt.scatter(G, G1, color="red")
plt.plot(surf_x, surf_y, color="black")
plt.plot([0.1, 1.0], [0.892, 0.892], color="black")
plt.plot([1.0, 1.0], [0.892, 0.1], color="black")
plt.plot([0.1, 1.0], [0.892, 0.1], color="black")
plt.text(0.1, 0.4, "l2_over = {:.3E}".format(l2_over), fontsize = 10)
plt.text(0.1, 0.3, "l2_real = {:.3E}".format(l2_real), fontsize = 10)
plt.text(0.1, 0.2, "l2_shape = {:.3E}".format(l2_shape), fontsize = 10)
plt.savefig('mixing_fd.png')
