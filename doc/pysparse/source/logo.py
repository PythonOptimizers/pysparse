import pylab

irow = [0, 1, 1, 1, 2, 2, 3, 3, 0]
jcol = [2, 0, 1, 3, 1, 2, 0, 2, 3]
nrow = ncol = 4
ms = 60

fig = pylab.figure()
ax = fig.gca()
ax.plot(jcol, irow, 'ks', markersize=ms, linestyle='None')
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(xmin=-1, xmax=ncol)
ax.set_ylim(ymin=nrow, ymax=-1)
ax.set_aspect('equal')
pylab.box(on=False)
pylab.show()
