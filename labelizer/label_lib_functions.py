# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 23:22:17 2020

@author: gebha
"""
import LabelLib as ll
import numpy as np


def savePqr(fileName, grid):
    points = grid.points()
    template = "ATOM{0: 7}   AV  AV{1: 6}{2:12.1f}{3:8.1f}{4:8.1f}{5:8.2f}{6:7.3f}\n"
    r = grid.discStep * 0.5
    with open(fileName, "w") as out:
        for i, p in enumerate(points.T):
            x, y, z, w = p
            sz = template.format(i + 1, 1, x, y, z, w, r)
            out.write(sz)


def savePqrFromAtoms(fileName, atoms):
    sz = "ATOM{0: 7}    X   X{1: 6}{2:12.1f}{3:8.1f}{4:8.1f}{5:8.2f}{6:7.3f}\n"
    with open(fileName, "w") as out:
        for i, at in enumerate(atoms.T):
            out.write(sz.format(i, i, at[0], at[1], at[2], 1.0, at[3]))


def saveXYZ(fileName, av):
    points = av.points()
    template = "F {0:.3f} {1:.3f} {2:.3f}\n"
    with open(fileName, "w") as out:
        header1 = str(len(points.T) + 1) + "\n"
        out.write(header1)
        header2 = "FILE AND FLUOROPHORE SPECIFICATION - TBD\n"
        out.write(header2)
        for i, p in enumerate(points.T):
            x, y, z, w = p
            sz = template.format(x, y, z)
            out.write(sz)

        x, y, z = meanAV(av)
        sz = "MP {0:.3f} {1:.3f} {2:.3f}\n".format(x, y, z)
        out.write(sz)


def saveAV(fileName, av):

    points = av.points()
    # print("POINTS SAVE")
    # print(av.points().T[:3])
    str_shape = "S {0:d} {1:d} {2:d}\n".format(*av.shape)
    str_origin = "O {0:.3f} {1:.3f} {2:.3f}\n".format(*av.originXYZ)
    str_delta = "D {0:.3f}\n".format(av.discStep)
    template = "P {0:.3f} {1:.3f} {2:.3f}\n"
    # r = grid.discStep * 0.5
    with open(fileName, "w") as out:
        out.write(str_shape)
        out.write(str_origin)
        out.write(str_delta)
        for i, p in enumerate(points.T):
            x, y, z, w = p
            sz = template.format(x, y, z)
            out.write(sz)


def loadAV(fileName):

    with open(fileName, "r") as fp:
        lines = fp.readlines()
        shape = np.array([int(d) for d in lines[0].split()[1:]])
        origin = np.array([float(d) for d in lines[1].split()[1:]])
        delta = float(lines[2].split()[1])
        length = shape[0] * shape[1] * shape[2]
        grid = -1 * np.ones(length)
        # print("SHAPE ", shape)
        for point_string in lines[3:]:
            coords = np.array([float(d) for d in point_string.split()[1:]])
            [nx, ny, nz] = [int(round(d)) for d in (coords - origin) / delta]
            # idx = nx*shape[1]*shape[2] + ny*shape[2] + nz (wrong order)
            idx = nz * shape[1] * shape[2] + ny * shape[1] + nx
            grid[idx] = 1
        av = ll.Grid3D(shape, origin, delta)
        av.grid = grid
        # print("POINTS LOAD")
        # print(av.points().T[:3])
        return av


def removeWeights(av):
    av_new = ll.Grid3D(av.shape, av.originXYZ, av.discStep)
    grid_new = np.array(av.grid)
    grid_new[grid_new > 0] = 1
    av_new.grid = grid_new
    return av_new


def calcAV(
    coords,
    cb_position,
    linker_length=20.0,
    linker_width=2.0,
    dye_radii=None,
    grid_size=1.0,
    save=False,
    filename=None,
):
    if dye_radii is None:
        dye_radii = [8.5, 3.5, 1.5]
    atoms = np.array(coords).astype(float).T
    simulation_grid_spacing = grid_size
    source = np.array(cb_position).astype(float)
    radii = np.array(dye_radii).astype(float)
    av = ll.dyeDensityAV3(
        atoms, source, linker_length, linker_width, radii, simulation_grid_spacing
    )
    # print(av.points().T[:3])
    if save and filename is not None:
        saveAV(filename, av)

    return av


def meanAV(av):
    return np.mean(av.points()[:3].T, axis=0)


def rmpDistance(av1, av2):
    return np.linalg.norm(meanAV(av1) - meanAV(av2))


def rdaDistance(av1, av2, n=10000):
    return ll.meanDistance(av1, av2, n)


def effAV(av1, av2, R0, n=10000):
    return ll.meanEfficiency(av1, av2, R0, n)


def effDistance(av1, av2, R0, n=10000):
    meanE = effAV(av1, av2, R0, n=n)
    return R0 * (1 / meanE - 1) ** (1 / 6)


def calcDistances(av1, av2, foerster_R):
    # distances = ll.sampleDistanceDistInv(av1,av2,100000)
    Dmean = ll.meanDistance(av1, av2, 100000)
    Emean = ll.meanEfficiency(av1, av2, foerster_R, 100000)

    return 1
