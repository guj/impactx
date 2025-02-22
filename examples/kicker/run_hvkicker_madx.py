#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl, Marco Garten
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


import amrex.space3d as amr
from impactx import ImpactX, RefPart, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of  nm
bunch_charge_C = 1.0e-9  # used without space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle().load_file("hvkicker.madx")

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=4.0e-3,
    sigmaY=4.0e-3,
    sigmaT=1.0e-3,
    sigmaPx=3.0e-4,
    sigmaPy=3.0e-4,
    sigmaPt=2.0e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
sim.lattice.load_file("hvkicker.madx", nslice=1)

# run simulation
sim.evolve()

# clean shutdown
del sim
amr.finalize()
