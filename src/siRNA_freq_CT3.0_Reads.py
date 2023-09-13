from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd

pd.options.display.width = 0

f = open('KM4_R1_001_viRNAs.txt', 'r')

df = pd.read_csv('KM4_R1_001_viRNAs.txt', sep="\t", header=None, names=["a", "b", "c", "d", "e", "f", "g", "h"])

df["sizes"] = df.apply(lambda row: len(row.e), axis=1)

sizes = df["sizes"].tolist()

freqs = Counter(df.d)

df2 = pd.DataFrame({"Size":df.sizes, "Position":df.d})

for i in range(len(df2)):
	if df2.loc[i].Size != 21:
		df2 = df2.drop(i)

plt.figure(figsize=(12, 6))

genome_length = 10836

# df2.sort_values("Position") for lowest and highest position of selected reads
# df.d.sort_values() for lowest and higest position of all reads
last_position = int(df.d.sort_values().tail(1))

# plt.hist(df2.Position, bins=genome_length, color='#FFB6C1', edgecolor="red")

plt.hist(df2.Position, bins=last_position, color='#FFB6C1', edgecolor="red")

plt.title("%i vsiRNAs\nLengths %i to %i"%(len(sizes),min(sizes),max(sizes)))

plt.xlabel("Position")

plt.ylabel("Reads")

plt.margins(x=0)

plt.savefig("fig1e.png", dpi=300, facecolor='white', edgecolor='black', orientation='potrait', format='png', transparent='True')

plt.show()

# while sorted.iloc[i] != 100 and i < len(sorted)-1:
#	i += 1
	    
