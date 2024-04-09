import numpy as np
from numpy.ma import masked_array
import matplotlib.pyplot as plt

# fname = "0_1_Jans/Av_Hamiltonian.txt"
# fname = "1nm_Skin/Av_Hamiltonian.txt"
fname = "0_2_Skin/Av_Hamiltonian.txt"
fhand = open(fname)

# get data
data = fhand.read()
data = data.strip().split()
data = data[1:]

# format data into a square
size = int(np.sqrt(0.25 + 2 * len(data)) - 0.5)
ham = np.zeros((size, size))
ham[np.triu_indices(size)] = data
ham = np.triu(ham) + np.tril(ham.T, -1)

# to allow for the log-type plot!
minimum = 0.001
ham[abs(ham) < minimum] = 0
diag = np.diag(ham)
logs = np.log2(np.abs(ham))
logs -= np.log2(minimum)
logs[ham < 0] *= -1
logs[np.diag_indices(size)] = diag


hamtoplot = ham  # ham or logs

off_diag = masked_array(hamtoplot, hamtoplot > 100)
on_diag = masked_array(hamtoplot, hamtoplot <= 100)

fig, ax = plt.subplots()
pa = ax.imshow(
    on_diag, interpolation='nearest', cmap=plt.cm.get_cmap('coolwarm_r'))
# cax: list of 4 items:
# [0] = horizontal location
# [1] = vertical location
# [2] = width of colored part of the colorbar
# [3] = height of the colored part of the colorbar
cba = plt.colorbar(pa, shrink=0.4, cax=fig.add_axes([0.8, 0.5, 0.03, 0.3]))

pb = ax.imshow(
    off_diag, interpolation='nearest', cmap=plt.cm.PRGn
)
cbb = plt.colorbar(pb, shrink=0.4, cax=fig.add_axes([0.8, 0.1, 0.03, 0.3]))
fig.subplots_adjust(right=0.75)

cba.set_label('frequency', rotation=0, y=1.12, labelpad=-22)
cbb.set_label('coupling', rotation=0, y=1.12, labelpad=-22)
ax.set_title("Hamiltonian using Skinner and TCC maps")
plt.show()