import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = "This program accepts a path to a file")

parser.add_argument("filename", type = str, help = "Path.txt")

args = parser.parse_args()

data = np.loadtxt(args.filename, unpack = "True")

cut = data[0]

tau = data[1]

sigma_tau = data[2]

plt.errorbar(cut, tau, sigma_tau, fmt = ".", capsize = 3)

plt.xlabel("Cut[ns]")

plt.ylabel("Muon Mealife [us]")

plt.grid(True)

plt.savefig("Cut_vs_ML")

plt.show()