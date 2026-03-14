import matplotlib.pyplot as plt
import numpy as np
import matplotlib.collections as mcoll
import matplotlib.path as mpath
from pylab import loadtxt

#these are my defaults that make pretty figures for posters
plot_params = {'axes.linewidth': 3,
               'xtick.labelsize': 'large',
               'ytick.labelsize': 'large',
               'xtick.major.size': 2,
               'xtick.minor.size': 1,
               'ytick.major.size': 2,
               'ytick.minor.size': 1,
               'xtick.direction': 'in',
               'ytick.direction': 'in',
               'patch.linewidth': 3,
               'font.size': 12,
               'font.family': 'serif',
               'font.serif': 'Times',
               'lines.linewidth': 2,
               'text.usetex': True,
	       'figure.subplot.right':0.9,
               'figure.subplot.left': 0.1,
               'figure.subplot.bottom': .1,
	       'figure.subplot.top':0.90,
               'figure.figsize': (6,6)}
plt.rcParams.update(plot_params)

def colorline(
    x, y, z=None, cmap=plt.get_cmap('copper'),
        linewidth=3, alpha=1.0):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap,
                              linewidth=linewidth, alpha=alpha)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

data3=loadtxt("PCyg2400.txt")
N = 10
np.random.seed(101)
x = data3[:,0]#
#x=np.random.rand(N)
y = data3[:,1]
#y=np.random.rand(N)
fig, ax = plt.subplots()

path = mpath.Path(np.column_stack([x, y]))
verts = path.interpolated(steps=3).vertices
x, y = verts[:, 0], verts[:, 1]
z = np.linspace(0, 1, len(x))
#plt.plot(x,z)
colorline(x, y, z, cmap=plt.get_cmap('rainbow'), linewidth=3)
ax.tick_params(right= True,top= True,which="both")
#PySplot
#vibgyor
base=15
shift=-5
plt.text(6550.5,base,"P",color='#7f4fc9',fontsize=100)
plt.text(6554,base,"y",color='#3e49bb',fontsize=80)
plt.text(6557,base+shift,"S",color='#526eff',fontsize=100)
plt.text(6560.5,base+shift,"p",color='#32c12c',fontsize=100)
plt.text(6564,base+shift,"l",color='#ffef00',fontsize=100)
plt.text(6566,base+shift,"o",color='#ff9a00',fontsize=100)
plt.text(6570,base+shift,"t",color='#d40c00',fontsize=100)
# plt.text(6550,19,"y",color='blue',fontsize=30)
#colorline(data3[:,0],data3[:,1],cmap='cubehelix', linewidth=3)
plt.axis((6550,6575,0,20))
ax.set_xticklabels([])
ax.set_yticklabels([])
# plt.axis('equal')
plt.savefig('icon.png')
plt.show()
