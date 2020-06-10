import numpy as np
import mdtraj as md
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
from mtools.gromacs.gromacs import make_comtrj


def calc_number_density(gro_file, trj_file, area,
        dim, box_range, n_bins, frame_range=None,maxs=None, mins=None):
    """
    Calculate a 1-dimensional number density profile for each residue

    Parameters
    ----------
    gro_file: str
        GROMACS '.gro' file to load 
    trj_file: str
        Trajectory to load
    area: int or float
        Area of dimensions not in number density profile
    dim: int
        Dimension to calculate number density profile (0,1 or 2)
    box_range: array
        Range of coordinates in 'dim' to evaluate
    frame_range: Python range() (optional)
        Range of frames to calculate number density function over
    maxs: array (optional)
        Maximum coordinate to evaluate 
    mins: array (optional)
        Minimum coordinate to evalute
    
    Attributes
    ----------
    """
    trj = md.load(trj_file, top=gro_file)
    resnames = np.unique([x.name for x in
                trj.topology.atoms])
    rho_list = list()
    res_list = list()
    
    for resname in resnames:
        sliced = trj.topology.select('name {}'.format(resname))
        trj_slice = trj.atom_slice(sliced)
        if frame_range:
            trj_slice = trj_slice[frame_range]
        for i,frame in enumerate(trj_slice):
            if maxs is None:
                indices = [[atom.index for atom in compound.atoms]
                          for compound in
                          list(frame.topology.residues)]
            else:
                indices = np.intersect1d(
                          np.intersect1d(np.where(frame.xyz[-1, :, 0]
                              > mins[0]),
                          np.where(frame.xyz[-1, :, 0] < maxs[0])),
                          np.intersect1d(np.where(frame.xyz[-1, :, 1] 
                              > box_range[0]),
                          np.where(frame.xyz[-1, :, 1] < box_range[1])))

            if frame_range:
                if i == 0:
                    x = np.histogram(frame.xyz[0,indices,dim].flatten(), 
                        bins=n_bins, range=(box_range[0], box_range[1]))
                    rho = x[0]
                    bins = x[1]
                else:
                    rho += np.histogram(frame.xyz[0, indices, dim].
                            flatten(),bins=n_bins, range=(box_range[0],
                                box_range[1]))[0]
            else:
                if i == 0:
                    x = np.histogram(frame.xyz[0,indices,dim].flatten(), 
                        bins=n_bins, range=(box_range[0], box_range[1]))
                    rho = x[0]
                    bins = x[1]
                else:
                    rho += np.histogram(frame.xyz[0, indices, dim].
                            flatten(),bins=n_bins, range=(box_range[0],
                                box_range[1]))[0]

        rho = np.divide(rho, trj_slice.n_frames*area*(bins[1]-bins[0]))
        rho_list.append(rho)
        res_list.append(resname)

    bin_list = bins[:-1]
    
    return (rho_list, bin_list, res_list)

# Set variabes for number density
area = 2.95 * 2.98
dim = 2
#box_range = [0.67, 2.17]
box_range = [10.83, 12.34]
n_bins = 100

# CASSANDRA
# Calculate number density function
rho, bins, residues = calc_number_density('coords.gro', 'wrapped.xyz', area, dim, box_range, n_bins)

# Plot number density
for i in [1,2]:
    plt.plot(bins, rho[i], label='Cassandra {}'.format(residues[i]))

# GROMACS
# Calculate number density function
rho, bins, residues = calc_number_density('gromacs/nvt_translate.gro',
    'gromacs/nvt_translate.trr', area, dim, box_range, n_bins)

# Plot number density
for i in [1,2]:
    plt.plot(bins, rho[i], label='GROMACS {}'.format(residues[i]))

plt.legend()
plt.ylabel(r'Number density ($nm^{-3}$)')
plt.ylim(0,110)
plt.xlabel('Position in slit pore (nm)')
plt.savefig('number-density.pdf')
