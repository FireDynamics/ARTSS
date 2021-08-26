import matplotlib.pyplot as plt
import numpy as np


def plot_colormap_vel(data, cmap='seismic', cb_label='velocity [m/s]', title='', savefig=False, filetitle='',fileformat='pdf'):
    extrema = max(np.amax(data), abs(np.amin(data)))
    plt.pcolormesh(data, cmap=cmap, vmin=-extrema, vmax=extrema)
    cb = plt.colorbar()
    cb.set_label(cb_label)
    plt.title(title)
    plt.gca().set_aspect('equal', adjustable='box')
    if savefig:
        if filetitle == '':
            filetitle='colormap_velocity'
        if not filetitle.endswith(fileformat):
            filetitle=filetitle+'.'+fileformat
        plt.savefig(filetitle,format=fileformat)
    else:
        plt.show()
    plt.close()


def plot_colormap_temp(data, cmap='viridis', cb_label='Temperature [K]', title='', savefig=False, filetitle='',fileformat='pdf'):
    plt.pcolormesh(data, cmap=cmap)
    cb = plt.colorbar()
    cb.set_label(cb_label)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(title)
    if savefig:
        if filetitle == '':
            filetitle='colormap_velocity'
        if not filetitle.endswith(fileformat):
            filetitle=filetitle+'.'+fileformat
        plt.savefig(filetitle,format=fileformat)
    else:
        plt.show()
    plt.close()

