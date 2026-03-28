# Feasibility Report: Refactoring MOCASSIN Atomic Modules

## Overview

This report evaluates the feasibility of reorganizing the MOCASSIN Monte Carlo radiative transfer code to separate the atomic physics calculations into element-specific modules. Specifically, the strategy entails:
1. Extracting Hydrogen and Helium into their own specific modules (`atom_H.f90`, `atom_He.f90`).
2. Creating a generic implementation module (`atom_heavy.f90`) that handles all heavy elements (Oxygen, Carbon, Nitrogen, Iron, etc.) using the standard matrix equilibrium solvers.

The constraints are that atomic data must continue to be loaded from external files, avoiding hard-coded physical data.

## 1. Current Architecture

Currently, MOCASSIN handles atomic physics via highly generalized, loop-driven routines located in `emission_mod.f90` and `ionization_mod.f90`.

- **Emission (`emission_mod.f90`)**: Calculates free-bound, free-free emission (`fb_ff`), two-photon emission (`twoPhoton`), recombination line emission (`RecLinesEmission`), and forbidden lines (`forLines`) by iterating over all elements (`elem = 1, nElements`) and ionization stages.
- **Opacities (`ionization_mod.f90`)**: Computes element opacities, continuum opacities, and solves for ionization states. Similar to emission, it uses generic loops over the `density` and `atomic_data_array` arrays.
- **Data Structures (`common_mod.f90` & `elements_mod.f90`)**: Key arrays like `HIRecLines`, `HeIRecLines`, `gammaHI`, and `forbiddenLines` are globally defined and used throughout the application, including inside `output_mod.f90` and `update_mod.f90`.

## 2. Evaluation of Proposed Strategy

The proposed strategy—extracting H, He, and a generic Heavy elements module—is **highly feasible and represents a solid architectural improvement**. It strikes a good balance between maintaining the generic data-driven approach and isolating element-specific edge cases.

### Pros of this Strategy:
*   **Separation of Concerns:** Hydrogen and Helium have unique physics in the codebase (e.g., `hydro2phot`, `HeI2photSub`, specific Lyman/Balmer logic, Benjamin, Skillman & Smits fits for He I). Isolating these into `atom_H.f90` and `atom_He.f90` removes a significant amount of specialized `if (elem == 1)` logic from the generic routines.
*   **Avoids Unnecessary Duplication:** By creating a single `atom_heavy.f90` module rather than individual modules per heavy element, the codebase avoids massive duplication of the generic `equilibrium` matrix solver and opacity additions.
*   **Data-Driven:** The `atom_heavy.f90` module can still seamlessly ingest the arrays populated in `elements_mod.f90` (from `fileNames.dat`).

### Required Code Changes

This refactor **would not require major structural changes to the core Monte Carlo simulation loop**, but it will require significant changes to the specific physics subroutines.

1.  **Module Creation:**
    *   Create `atom_H.f90`, `atom_He.f90`, and `atom_heavy.f90`.
2.  **Relocating Logic:**
    *   From `emission_mod.f90`, move H-specific logic from `fb_ff`, `twoPhoton`, and `RecLinesEmission` into `atom_H.f90`.
    *   Move He-specific logic into `atom_He.f90`.
    *   Move the `equilibrium` subroutine and `forLines` logic into `atom_heavy.f90`.
    *   From `ionization_mod.f90`, extract the H/He continuum opacity calculations and move the `putOpacity` heavy-element logic to `atom_heavy.f90`.
3.  **Delegation Strategy:**
    *   The `emissionDriver` and `ionizationDriver` routines would be transformed into delegators. Instead of containing the logic, they would call `atom_H_emission(...)`, `atom_He_emission(...)`, and then loop `atom_heavy_emission(elem, ...)` for the rest.
4.  **Handling Global Variables:**
    *   Because arrays like `HIRecLines` and `forbiddenLines` are used widely (e.g., in `output_mod.f90`), these arrays should either remain in `common_mod.f90` and be passed to or updated by the new modules, or the new modules must expose getter functions to provide this data to the output routines. Keeping them in `common_mod.f90` will minimize architectural disruption.

## 3. Conclusion

The strategy of creating `atom_H.f90`, `atom_He.f90`, and `atom_heavy.f90` is sound and feasible. It **can be done without major problems or disrupting the core architecture**. The changes are mostly limited to extracting tightly coupled subroutine logic from `emission_mod.f90` and `ionization_mod.f90` into the new modules and updating the drivers to delegate the calculations appropriately.