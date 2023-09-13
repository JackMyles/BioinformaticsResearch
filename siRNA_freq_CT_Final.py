from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def read_counter(file, save_file, len_range):
    pd.options.display.width = 0

    df = pd.read_csv(file, sep="\t", header=None, names=["a", "b", "c", "d", "e", "f", "g", "h"])

    df["sizes"] = df.apply(lambda row: len(row.e), axis=1)

    sizes = df["sizes"].tolist()

    freqs = Counter(df.sizes)
   
    plt.hist(sizes, bins=25, color='#FFB6C1', edgecolor="black")

    if len_range[0] != -1:
        plt.xlim(len_range[0], len_range[1])
        plt.title("%i vsiRNAs\nLengths %i to %i"%(len(sizes),min(len_range),max(len_range)))

    plt.title("%i vsiRNAs\nLengths %i to %i"%(len(sizes),min(sizes),max(sizes)))

    plt.xlabel("Sequence length (bp)")

    plt.ylabel("Count")

    plt.savefig(save_file, dpi=300, facecolor='white', edgecolor='black', orientation='potrait', format='png', transparent='True')

    plt.show()

def double_bar(file, save_file):
    pd.options.display.width = 0

    df = pd.read_csv(file, sep="\t", header=None, names=["a", "b", "c", "d", "e", "f", "g", "h"])

    df["sizes"] = df.apply(lambda row: len(row.e), axis=1)

    sizes = df["sizes"].tolist()

    freqs = Counter(df.sizes)

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

    plt.savefig(save_file, dpi=300, facecolor='white', edgecolor='black', orientation='potrait', format='png', transparent='True')

    plt.show()

def dual_bar(file, len_range, save_file):
    pd.options.display.width = 0

    df = pd.read_csv(file, sep="\t", header=None, names=["a", "b", "c", "d", "e", "f", "g", "h"])

    df["sizes"] = df.apply(lambda row: len(row.e), axis=1)

    sizes = df["sizes"].tolist()

    sizes = np.array(sizes)

    freqs = Counter(df.d)
    #freqs = sorted(freqs.items(), key=lambda i: i[1])

    df2 = pd.DataFrame({"Size":df.sizes, "Position":df.d, "+/-":df.b})

    for i in range(len(df2)):
            if df2.loc[i].Size not in range(len_range[0], len_range[1]+1):
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

    plt.savefig(save_file, dpi=300, facecolor='white', edgecolor='black', orientation='potrait', format='png', transparent='True')

    plt.show()

def stacked_bar(file, pos_neg, len_range, save_file):
    pd.options.display.width = 0
    
    df = pd.read_csv(file, sep="\t", header=None, names=["a", "b", "c", "d", "e", "f", "g", "h"])

    df["sizes"] = df.apply(lambda row: len(row.e), axis=1)

    sizes = df["sizes"].tolist()

    freqs = Counter(df.d)

    df2 = pd.DataFrame({"Size":df.sizes, "Position":df.d, "String":df.e, "+/-":df.b})

    if pos_neg == "+":
        pos_df2 = df2 = df2[df2['+/-'] != '-']
        df2 = pos_df2
    elif pos_neg == "-":
        neg_df2 = df2 = df2[df2['+/-'] != '+']
        df2 = neg_df2

    a = np.zeros(10)
    c = np.zeros(10)
    t = np.zeros(10)
    g = np.zeros(10)

    for i in range(len(df2)):
        if df2.loc[i].Size not in range(len_range[0], len_range[1]+1):
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

    plt.title("%i vsiRNAs\nLengths %i to %i"%(len(df2),min(len_range),max(len_range)))

    plt.xlabel("(nt)")

    plt.ylabel("%")

    plt.legend(["A", "C", "U", "G"], bbox_to_anchor=(1.0, 1.0), loc='upper left')

    plt.margins(x=0.01)

    plt.savefig(save_file, dpi=300, facecolor='white', edgecolor='black', orientation='potrait', format='png', transparent='True')

    plt.show()


file_found = False
while file_found != True:
    file = input("Enter the file name: ")
    try:
        f = open(file, 'r')
        file_found = True
    except FileNotFoundError:
        print("\nFile does not exist\n")
        file_found = False
    except:
        print("\nError Occured Opening File\n")
        file_found = False

choice = int(input(
'''\n1) Read Counter
2) Double Bar (+ vs -)
3) Split Bar
4) Stacked Bar

Enter the number of the graph you want to create: '''))

if choice == 1:
    len_range = input('''\nEnter the first and last number of the range of lengths you
want included seperated by a space (if you want the default values enter -1): ''')
    len_range = len_range.split()
    for i in range(len(len_range)):
        len_range[i] = int(len_range[i])

    save_file = input('''\nWhat do you want the saved file to be named: ''')

    read_counter(file, save_file, len_range)

elif choice == 2:
    save_file = input('''\nWhat do you want the saved file to be named: ''')

    double_bar(file, save_file)

elif choice == 3:
    len_range = input('''\nEnter the first and last number of the range of lengths you
want included seperated by a space (if only one number, enter only once): ''')
    len_range = len_range.split()
    for i in range(len(len_range)):
        len_range[i] = int(len_range[i])
    if len(len_range) == 1:
        len_range.append(len_range[0])

    save_file = input('''\nWhat do you want the saved file to be named: ''')
    
    dual_bar(file, len_range, save_file)
    
elif choice == 4:
    pos_neg = input('''\nShould the graph include positive, negative, or both strands?
(enter '+' for positive, '-' for negative, and '+-' for both): ''')

    len_range = input('''\nEnter the first and last number of the range of lengths you
want included seperated by a space (if only one number, enter only once): ''')
    len_range = len_range.split()
    for i in range(len(len_range)):
        len_range[i] = int(len_range[i])
    if len(len_range) == 1:
        len_range.append(len_range[0])

    save_file = input('''\nWhat do you want the saved file to be named: ''')
    
    stacked_bar(file, pos_neg, len_range, save_file)
    
    
