"""
This aims Python hook plots the electron density in the xy plane. It needs
to be run in parallel and requires matplotlib.
"""
import numpy as np
from scipy.interpolate import griddata
import matplotlib
matplotlib.use('Agg')  # no interactive plots, only PNG files
import matplotlib.pyplot as plt  # needs to be imported after matplotlib.use
from mpl_toolkits.mplot3d import axes3d  # 3D plots with projections
from matplotlib import cm  # color palettes

bohr = 0.52917721067


def run(ctx):
    # `ctx` is a context object which contains various quantities from within
    # AIMS. Here, we use a convenience function to gather the grid points and
    # electron density from all MPI tasks.
    points, rho, partition_tab = ctx.gather_all_grids(['rho', 'partition_tab'])
    if ctx.rank == 0:  # the grids are gather only on the root process
        points = points.T  # (3, npts) -> (npts, 3)
        rho = rho.sum(0)  # sum over spin
        in_plane = abs(points[:, 2]) < 1e-10  # points in xy plane
        points = bohr*points[in_plane, 0:2]  # filter & take x, y & scale
        rho = np.log10(1+rho[in_plane])  # filter & log scale
        X, Y = np.mgrid[-4:4:400j, -2:2:200j]  # get rectangular grid
        # interpolate density to rectangular grid
        rho = griddata(points, rho, (X, Y), method='cubic')
        fig = plt.figure()  # start figure
        ax = fig.gca(projection='3d')  # get axes
        ax.plot_surface(X, Y, rho, alpha=0.3)
        ax.contour(X, Y, rho, np.log10(1+10.0**np.array(range(-3, 3))),
                   zdir='z', offset=-1, cmap=cm.coolwarm)
        ax.contour(X, Y, rho, [-2, -1.8, -1.6, -1],
                   zdir='x', offset=-5, cmap=cm.coolwarm)
        ax.contour(X, Y, rho, [0, 0.2, 0.4, 1],
                   zdir='y', offset=3, cmap=cm.coolwarm)
        ax.set_xlim(-5, 4)
        ax.set_ylim(-2, 3)
        ax.set_zlim(-1, None)
        fig.savefig('density.png')
