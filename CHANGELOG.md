# Atomic Modules Refactoring Changelog

This document tracks the extraction of atomic physics logic from generalized loops into element-specific modules (`atom_H.f90`, `atom_He.f90`, `atom_heavy.f90`).

## Strategy
- Retain global variables (`HIRecLines`, `forbiddenLines`, etc.) in `common_mod.f90` to prevent breaking downstream modules (`output_mod.f90`, `update_mod.f90`).
- Extract specific heating, cooling, emission, and opacity calculations for elements into their own modules.
- Refactor `emission_mod.f90` and `ionization_mod.f90` drivers to delegate to these new modules.

## Changes

- Created `atom_H.f90` to handle hydrogen-specific physics (`atom_H_fb_ff`, `hydro2phot`, `atom_H_RecLinesEmission`, `atom_H_inOpacity`).
- Created `atom_He.f90` to handle helium-specific physics (`atom_He_fb_ff`, `atom_He_twoPhoton`, `atom_He_RecLinesEmission`, `atom_He_inOpacity`).
- Created `atom_heavy.f90` to handle heavy elements physics (`atom_heavy_fb_ff`, `atom_heavy_forLines`, `equilibrium`, `atom_heavy_putOpacity`).
- Cleaned up duplicated physics logic inside `output_mod.f90`, `mocassinPlot.f90`, and `update_mod.f90` to correctly use the new logic structures.
