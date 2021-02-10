#!/usr/bin/python

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

if len(sys.argv) != 4:
    print("Please provide three arguments: folder, file name, and library ID!")
    sys.exit()

folder = sys.argv[1]
filename = sys.argv[2]
library = sys.argv[3]

file1 = "{}/{}.unique.ratio".format(folder, filename)
file2 = "{}/{}.multi.ratio".format(folder, filename)
file3 = "{}/{}.unique.score".format(folder, filename)
file4 = "{}/{}.multi.score".format(folder, filename)
file5 = "{}/{}.unique.mismatch".format(folder, filename)
file6 = "{}/{}.multi.mismatch".format(folder, filename)

pp = PdfPages("{}/{}_alignment_quality.pdf".format(folder, library))

# plot alignment ratio
df_x = np.loadtxt(file1, delimiter="\t", dtype="int", usecols=0)
df_y = np.loadtxt(file1, delimiter="\t", dtype="int", usecols=1)
df = pd.DataFrame({"x": df_x, "y": df_y})
fig, ax = plt.subplots(figsize=(8, 8))
plt.bar(df["x"], df["y"], width=0.7, color="lightskyblue", edgecolor="black")
plt.yticks(rotation=90)
plt.xlabel("alignment ratio %")
plt.title("Histogram of unique alignment ratio")
plt.savefig(pp, format="pdf")

df_x = np.loadtxt(file2, delimiter="\t", dtype="int", usecols=0)
df_y = np.loadtxt(file2, delimiter="\t", dtype="int", usecols=1)
df = pd.DataFrame({"x": df_x, "y": df_y})
fig, ax = plt.subplots(figsize=(8, 8))
plt.bar(df["x"], df["y"], width=0.7, color="lightskyblue", edgecolor="black")
plt.yticks(rotation=90)
plt.xlabel("alignment ratio %")
plt.title("Histogram of multiple alignment ratio")
plt.savefig(pp, format="pdf")

# plot alignment score
df_x = np.loadtxt(file3, delimiter="\t", dtype="int", usecols=0)
df_y = np.loadtxt(file3, delimiter="\t", dtype="int", usecols=1)
df = pd.DataFrame({"x": df_x, "y": df_y})
fig, ax = plt.subplots(figsize=(8, 8))
plt.bar(df["x"], df["y"], width=0.7, color="lightskyblue", edgecolor="black")
plt.yticks(rotation=90)
plt.xlabel("alignment score")
plt.title("Histogram of unique alignment score")
plt.savefig(pp, format="pdf")

df_x = np.loadtxt(file4, delimiter="\t", dtype="int", usecols=0)
df_y = np.loadtxt(file4, delimiter="\t", dtype="int", usecols=1)
df = pd.DataFrame({"x": df_x, "y": df_y})
fig, ax = plt.subplots(figsize=(8, 8))
plt.bar(df["x"], df["y"], width=0.7, color="lightskyblue", edgecolor="black")
plt.yticks(rotation=90)
plt.xlabel("alignment score")
plt.title("Histogram of multiple alignment score")
plt.savefig(pp, format="pdf")

# plot alignment mismatch
df_x = np.loadtxt(file5, delimiter="\t", dtype="int", usecols=0)
df_y = np.loadtxt(file5, delimiter="\t", dtype="int", usecols=1)
df = pd.DataFrame({"x": df_x, "y": df_y})
fig, ax = plt.subplots(figsize=(8, 8))
plt.bar(df["x"], df["y"], width=0.7, color="lightskyblue", edgecolor="black")
plt.yticks(rotation=90)
plt.xlabel("alignment mismatch")
plt.title("Histogram of unique alignment mismatch")
plt.savefig(pp, format="pdf")

df_x = np.loadtxt(file6, delimiter="\t", dtype="int", usecols=0)
df_y = np.loadtxt(file6, delimiter="\t", dtype="int", usecols=1)
df = pd.DataFrame({"x": df_x, "y": df_y})
fig, ax = plt.subplots(figsize=(8, 8))
plt.bar(df["x"], df["y"], width=0.7, color="lightskyblue", edgecolor="black")
plt.yticks(rotation=90)
plt.xlabel("alignment mismatch")
plt.title("Histogram of multiple alignment mismatch")
plt.savefig(pp, format="pdf")

pp.close()
