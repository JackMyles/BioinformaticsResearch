import scipy
from pylab import figure

def bihist(y1, y2, nbins=10, h=None):
        '''
        Bihistogram.
        h is an axis handle. If not present, a new figure is created.
        '''
        if h is None: h = figure().add_subplot(111)
        xmin = scipy.floor(scipy.minimum(y1.min(), y2.min()))
        xmax = scipy.ceil(scipy.maximum(y1.max(), y2.max()))
        bins = scipy.linspace(xmin, xmax, nbins)
        n1, bins1, patch1 = h.hist(y1, bins)
        n2, bins2, patch2 = h.hist(y2, bins)
        # set ymax:
        ymax = 0
        for i in patch1:
                height = i.get_height()
                if height > ymax: ymax = height
        # invert second histogram and set ymin:
        ymin = 0
        for i in patch2:
                height = i.get_height()
                height = -height
                i.set_height(height)
                if height < ymin: ymin = height
        h.set_ylim(ymin*1.1, ymax*1.1)          
        h.figure.canvas.draw()

y1 = [0, 1, 2, 3, 4, 5]
y2 = [2, 3, 4, 5, 6, 7]

bihist(y1, y2)
