from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

pd.options.display.width = 0

f = open('KM4_R1_001_viRNAs.txt', 'r')

#df = pd.DataFrame([f.readline().split()], columns= ["a", "b", "c", "d", "e", "f", "g", "h"])

df = pd.read_csv('KM4_R1_001_viRNAs.txt', sep="\t", header=None, names=["a", "b", "c", "d", "e", "f", "g", "h"])

df["sizes"] = df.apply(lambda row: len(row.e), axis=1)

sizes = df["sizes"].tolist()

freqs = Counter(df.sizes)
#freqs = sorted(freqs.items(), key=lambda i: i[1])

pos_sizes = []
neg_sizes = []

for i in range(df.b.size):
    if (df.b[i] == '+'):
        pos_sizes.append(df.sizes[i])
        neg_sizes.append(0)
    if (df.b[i] == '-'):
        neg_sizes.append(df.sizes[i])
        pos_sizes.append(0)

df2 = pd.DataFrame({"Positive":pos_sizes, "Negative":neg_sizes})

pos_freqs = Counter(df2.Positive)
neg_freqs = Counter(df2.Negative)

df3 = pd.DataFrame({"Positive":pos_freqs, "Negative":neg_freqs})
df3 = df3.sort_index()
df3 = df3.drop(index=0)
for i in df3.index:
    if i > 40 or i < 18:
        df3 = df3.drop(index=i)
df3 = df3.convert_dtypes(convert_integer=True)
ax = df3.plot.bar(color=["Red","Blue"], rot=0, title="Placeholder", width=0.8, figsize=(10,5))
ax.set_xlabel("X")
ax.set_ylabel("Y")
#ax.tick_params(rotation=45)

plt.show()

"""
width = 0.3
x = ["21", "22", "23", "24"]
X_axis
fig, ax = plt.subplots()

bar1 = ax.bar(x-width/2, pos_sizes, width, label="+")
bar2 = ax.bar(x+width/2, neg_sizes, width, label="-")

plt.show()
"""

"""
plt.hist(sizes, bins=25, color='#FFB6C1', edgecolor="black")

plt.title("%i vsiRNAs\nLengths %i to %i"%(len(sizes),min(sizes),max(sizes)))

plt.xlabel("Sequence length (bp)")

plt.ylabel("Count")

plt.savefig("fig1d.png", dpi=300, facecolor='white', edgecolor='black', orientation='potrait', format='png', transparent='True')
"""
