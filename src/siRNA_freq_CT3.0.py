from collections import Counter
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

plt.hist(sizes, bins=25, color='#FFB6C1', edgecolor="black")

plt.title("%i vsiRNAs\nLengths %i to %i"%(len(sizes),min(sizes),max(sizes)))

plt.xlabel("Sequence length (bp)")

plt.ylabel("Count")

plt.savefig("fig1d.png", dpi=300, facecolor='white', edgecolor='black', orientation='potrait', format='png', transparent='True')
