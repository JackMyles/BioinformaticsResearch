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

freqs = Counter(df.d)
#freqs = sorted(freqs.items(), key=lambda i: i[1])

df2 = pd.DataFrame({"Size":df.sizes, "Position":df.d, "+/-":df.b})

for i in range(len(df2)):
	if df2.loc[i].Size != 21:
		df2 = df2.drop(i)

positions = df2["Position"].tolist()

pos_reads = []
neg_reads = []

for i in range(df2["+/-"].size):
    if (df2["+/-"].iloc[i] == '+'):
        pos_reads.append(positions[i])
        neg_reads.append(0)
    if (df2["+/-"].iloc[i] == '-'):
        neg_reads.append(positions[i])
        pos_reads.append(0)

#sizes = df2["Size"].tolist()
pos_reads = np.array(pos_reads)
freq_pos = Counter(pos_reads)
pos_count = []
for i in range(len(positions)):
    pos_count.append(freq_pos[positions[i]])
pos_count = np.array(pos_count)
#sizes2 = -1 * sizes
neg_reads = -1 * np.array(neg_reads)
freq_neg = Counter(neg_reads)
neg_count = []
for i in range(len(positions)):
    neg_count.append(freq_neg[-positions[i]])
neg_count = -1 * np.array(neg_count)

f, ax = plt.subplots(figsize=(12,6))

ax.bar(positions, pos_count, color='red', edgecolor='red')
ax.bar(positions, neg_count, color='blue', edgecolor='blue')

ax.margins(x=0)

plt.ylim([-10, 40])

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
