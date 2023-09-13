from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

pd.options.display.width = 0

f = open('KM4_R1_001_viRNAs.txt', 'r')

#df = pd.DataFrame([f.readline().split()], columns= ["a", "b", "c", "d", "e", "f", "g", "h"])

df = pd.read_csv('KM4_R1_001_viRNAs.txt', sep="\t", header=None, names=["a", "b", "c", "d", "e", "f", "g", "h"])

df["sizes"] = df.apply(lambda row: len(row.e), axis=1)

sizes = df["sizes"].tolist()

sizes = np.array(sizes)

freqs = Counter(df.sizes)
#freqs = sorted(freqs.items(), key=lambda i: i[1])

df2 = pd.DataFrame({"Size":df.sizes, "Position":df.d, "+/-":df.b})

for i in range(len(df2)):
	if df2.loc[i].Size != 21:
		df2 = df2.drop(i)

positions = df2["Position"].tolist()

pos_sizes = []
neg_sizes = []

for i in range(df2["+/-"].size):
    if (df2["+/-"].iloc[i] == '+'):
        pos_sizes.append(df.sizes[i])
        neg_sizes.append(0)
    if (df2["+/-"].iloc[i] == '-'):
        neg_sizes.append(df.sizes[i])
        pos_sizes.append(0)

#sizes = df2["Size"].tolist()
pos_sizes = np.array(pos_sizes)
#sizes2 = -1 * sizes
neg_sizes = -1 * np.array(neg_sizes)

f, ax = plt.subplots()

ax.bar(positions, pos_sizes)
ax.bar(positions, neg_sizes)

# Formatting x labels
plt.xticks(rotation=90)
plt.tight_layout()
# Use absolute value for y-ticks
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])

plt.show()

#plt.xlabel("Sequence length (bp)")

#plt.ylabel("Count")

#plt.savefig("fig1d.png", dpi=300, facecolor='white', edgecolor='black', orientation='potrait', format='png', transparent='True')
