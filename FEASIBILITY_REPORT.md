# Feasibility Report: Refactoring MOCASSIN Atomic and Dust Modules

## Overview

This report evaluates the feasibility of reorganizing the MOCASSIN Monte Carlo radiative transfer code to separate the physics calculations into element-specific modules (e.g., `atom_H.f90`, `atom_O.f90`), and consolidating dust-related physics into a generic module (e.g., `dust.f90`). The primary constraints are that atomic data must continue to be loaded from external files, and the dust module must handle generalized physics (collisions, charge exchange, etc.) without hardcoded species data.

## 1. Feasibility of Element-Specific Modules (`atom_H.f90`, `atom_O.f90`)

### Current Architecture
Currently, MOCASSIN handles atomic physics via highly generalized, loop-driven routines located in `emission_mod.f90` and `ionization_mod.f90`.
- **`emission_mod.f90`**: Calculates free-bound, free-free emission (`fb_ff`), two-photon emission (`twoPhoton`), recombination line emission (`RecLinesEmission`), and forbidden lines (`forLines`) by iterating over all elements (`elem = 1, nElements`) and ionization stages.
- **`ionization_mod.f90`**: Computes element opacities, continuum opacities, and solves for ionization states using loops over generic data structures (`atomic_data_array`).
- **Data Model**: External files are read into global arrays (e.g., `xSecArray`, `atomic_data_array`) in `elements_mod.f90`.

### Feasibility: Partially Feasible but Requires Caution

It is technically feasible to break down these generic loops into element-specific modules, but there are significant architectural and performance trade-offs to consider.

**Pros:**
*   **Separation of Concerns:** Element-specific edge cases (like hydrogen/helium free-free calculations, or specific recombination line fits like Benjamin, Skillman & Smits for He I) can be isolated rather than cluttering generic routines with `if (elem == 1)` or `select case (z)` blocks.
*   **Maintainability:** Easier to locate and update the specific physics or fitting formulas for a single element.

**Cons & Challenges:**
*   **Code Duplication:** Much of the physics (e.g., calculating collision strengths, matrix inversions for forbidden lines) is identical across heavy elements (Oxygen, Carbon, Nitrogen, etc.). Creating an `atom_O.f90` and `atom_C.f90` could lead to massive code duplication unless an Object-Oriented or highly modular approach is taken.
*   **Performance Overhead:** The code relies heavily on vectorized operations and nested loops over flat arrays. Moving to polymorphic derived types or making individual subroutine calls per element/cell could introduce significant function-call overhead and break cache locality, which is crucial for a computationally expensive 3D Monte Carlo code.

### Proposed Architecture for Atomic Refactoring

Instead of hardcoding `atom_O.f90`, `atom_C.f90` from scratch, we should leverage Fortran 2003+ Object-Oriented features (or a module-based interface design) to achieve the desired separation without losing the data-driven flexibility.

1.  **Base Atomic Class (`atom_base.f90`)**: Create a base module with an abstract type `Atom` that defines the interfaces for heating, cooling, emission, and opacity calculations. It will contain the generalized routines currently found in `emission_mod.f90` and `ionization_mod.f90`.
2.  **Element-Specific Implementations**:
    *   **`atom_H.f90` & `atom_He.f90`**: These elements already have highly specific physics (e.g., `hydro2phot`, `HeI2photSub`, specific Lyman/Balmer routines). They would extend the `Atom` type and override the specific methods.
    *   **`atom_heavy.f90`**: A generic implementation that handles all heavy elements (O, C, N, Fe) using the standard matrix equilibrium solvers (`equilibrium` subroutine). This avoids creating `atom_O.f90` when oxygen's logic is identical to nitrogen's logic.
3.  **Data Loading**: The data will still be loaded dynamically from `data/fileNames.dat` into the instances of these types at runtime.

## 2. Feasibility of a Generic Dust Module (`dust.f90`)

### Current Architecture
The dust component is currently split. `dust_mod.f90` calculates basic dust opacities and handles the MPI reduction of these opacities. However, complex dust physics like Quantum Heating (`qHeat`), enthalpy calculations (`enthalpy`), and dust emission PDFs (`setDustPDF`) are currently embedded inside the gas-focused `emission_mod.f90`.

### Feasibility: Highly Feasible and Recommended

Consolidating the dust physics is very feasible and would greatly improve the codebase's organization.

**Proposed Architecture for Dust Refactoring**

1.  **Consolidated `dust_physics_mod.f90`**:
    *   Move the dust-specific routines out of `emission_mod.f90` (e.g., `qHeat`, `getTmin`, `enthalpy`, `clrate`, `cheat`, `setDustPDF`) into the new dust module.
    *   This module will act as the single source of truth for all general dust physics (heating, cooling, collisions, charge exchange).
2.  **Data-Driven Species**:
    *   The physics routines will remain completely generic. They will rely on the existing data arrays (like `dustAbsXsecP`, `dustScaXsecP`, `grainRadius`, `grainWeight`) which are populated from the external `dustData` directory.
    *   The existing hardcoded specific species logic (like the `select case(sorc)` checking for 'S' [Silicate] or 'C' [Carbon] in the `enthalpy` and `qHeat` routines) should be refactored. The properties that distinguish Silicate from Carbon (e.g., natom scaling factors, enthalpy polynomial coefficients) should ideally be read from a configuration file rather than hardcoded in a `case` statement. This ensures the module is truly species-agnostic.

## Conclusion

The proposed refactoring is feasible.

*   For the **Atomic components**, extracting Hydrogen and Helium into specific modules (`atom_H.f90`, `atom_He.f90`) is highly beneficial due to their unique physics. For heavier elements, a generic `atom_heavy.f90` module driven by external data is recommended over creating individual files for Oxygen, Carbon, etc., to prevent code duplication.
*   For the **Dust component**, consolidating all scattered dust physics from `emission_mod.f90` into a central, data-driven `dust.f90` module is both feasible and highly recommended for code clarity.
