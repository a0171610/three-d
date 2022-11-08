import numpy as np
import matplotlib.pyplot as plt
import japanize_matplotlib

x = [4, 3, 2, 1]
y1 = [0.201, 0.111, 0.0232, 0.004]
y2 = [0.381, 0.241, 0.131, 0.017]
y4 = [0.150, 0.053, 0.006, 0.0019]
y5 = [0.033333, 0.033333, 0.033333, 0.033333]

fig = plt.figure(facecolor="white")
ax = fig.add_subplot(111)
ax.set_ylim(0.0001, 1.0)
plt.xticks([4, 3, 2, 1], ["3.0", "1.5", "0.75", "0.375"])
ax.set_yticks([1.0, 0.1, 0.01, 0.001, 0.0001])
ax.set_yscale('log')
ax.grid()
ax.invert_xaxis()
ax.plot(x, y1, marker="o", label="内挿なし(4点)", color="blue", markersize=8)
ax.plot(x, y2, marker="o", label="内挿なし(1点)", color="#ff7f00", markersize=8)
ax.plot(x, y4, marker="o", label="内挿あり", color="green", markersize=8)
ax.plot(x, y5, color="purple", lw=4)
ax.set_xlabel('$\Delta\lambda$', fontsize=18)
ax.set_ylabel('$l_2$', fontsize=18)
ax.legend()
ax.set_title("resolution dependency of $l_2$ norm")
plt.show()