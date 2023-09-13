from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

pd.options.display.width = 0

f = open('KM4_R1_001_viRNAs.txt', 'r')

df = pd.read_csv('KM4_R1_001_viRNAs.txt', sep="\t", header=None, names=["a", "b", "c", "d", "e", "f", "g", "h"])

df["sizes"] = df.apply(lambda row: len(row.e), axis=1)

sizes = df["sizes"].tolist()

freqs = Counter(df.d)

df2 = pd.DataFrame({"Size":df.sizes, "Position":df.d, "String":df.e, "+/-":df.b})

#pos_df2 = df2 = df2[df2['+/-'] != '-']
#neg_df2 = df2 = df2[df2['+/-'] != '+']
df2.reset_index()

a = np.zeros(10)
c = np.zeros(10)
t = np.zeros(10)
g = np.zeros(10)

for i in range(len(df2)):
    if df2.loc[i].Size != 21:
        df2 = df2.drop(i)
    else:
        string = df2.loc[i].String
        for i in range(min(10, len(string))):
            if string[i] == "A":
                a[i] += 1
            elif string[i] == "C":
                c[i] += 1
            elif string[i] == "T":
                t[i] += 1
            elif string[i] == "G":
                g[i] += 1

for i in range(len(a)):
    total = a[i] + c[i] + t[i] + g[i]
    a[i] = a[i]/total * 100
    c[i] = c[i]/total * 100
    g[i] = g[i]/total * 100
    t[i] = t[i]/total * 100

plt.figure(figsize=(12, 6))

x = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
# y = ["1", "10"]

plt.bar(x, a, color='#e87e72', width=0.9)
plt.bar(x, c, bottom=a, color='#62b685', width=0.9)
plt.bar(x, t, bottom=a+c, color='#54bdc2', width=0.9)
plt.bar(x, g, bottom=a+c+t, color='#bc81f8', width=0.9)

# a2 = [a[0], a[9]]
# c2 = [c[0], c[9]]
# g2 = [g[0], g[9]]
# t2 = [t[0], t[9]]

genome_length = 10836

# df2.sort_values("Position") for lowest and highest position of selected reads
# df.d.sort_values() for lowest and higest position of all reads

plt.title("%i vsiRNAs\nLengths %i to %i"%(len(sizes),min(sizes),max(sizes)))

plt.xlabel("(nt)")

plt.ylabel("%")

plt.legend(["A", "C", "T", "G"], bbox_to_anchor=(1.0, 1.0), loc='upper left')

plt.margins(x=0.01)

plt.savefig("fig1f.png", dpi=300, facecolor='white', edgecolor='black', orientation='potrait', format='png', transparent='True')

plt.show()
