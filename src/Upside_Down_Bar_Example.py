from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

y = ['{} to {} years'.format(i, i+4) for i in range(0, 90, 4)]
d_1 = np.random.randint(0, 150, 23)
d_2 = -1 * d_1

fig, ax = plt.subplots()
ax.bar(y, d_1)
ax.bar(y, d_2)

# Formatting x labels
plt.xticks(rotation=90)
plt.tight_layout()
# Use absolute value for y-ticks
ticks =  ax.get_yticks()
ax.set_yticklabels([int(abs(tick)) for tick in ticks])

plt.show()
