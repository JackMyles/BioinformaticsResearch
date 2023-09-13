import matplotlib.pyplot as plt
import pandas as pd
times = pd.date_range('2018-09-01', periods=7, freq='5D')
yesSeries = pd.Series([1800,2000,3000,1000,2000,1500,1700], index=times)
nodSeries = pd.Series([200,500,700,600,300,50,0], index=times)

df = pd.DataFrame({"Example one":yesSeries,"Example two":nodSeries})
ax = df.plot.bar(color=["SkyBlue","IndianRed"], rot=0, title="Epic Graph\nAnother Line! Whoa")
ax.set_xlabel("date")
ax.set_ylabel("counts")
ax.xaxis.set_major_formatter(plt.FixedFormatter(times.strftime("%b %d %Y")))
plt.show()
