# mocassin_HM_2025: Monte Carlo Simulations of Ionised Nebulae

**MOCASSIN** is a fully 3D or 2D photoionisation and dust radiative transfer code which employs a Monte Carlo approach to the transfer of radiation through media of arbitrary geometry and density distribution. It was originally developed for the modelling of photoionised regions like HII regions and planetary nebulae and has since expanded and been applied to a variety of astrophysical problems, including modelling clumpy dusty supernova envelopes, star forming galaxies, protoplanetary disks and inner shell fluorescence emission in the photospheres of stars and disk atmospheres.

The code can deal with arbitrary Cartesian grids of variable resolution, it has successfully been used to model complex density fields from SPH calculations and can deal with ionising radiation extending from Lyman edge to the X-ray. The dust and gas microphysics is fully coupled both in the radiation transfer and in the thermal balance.

## Overview

The project provides four main executables:

- `mocassin`: The main MOCASSIN driver that is used to start a new simulation.
- `mocassinWarm`: Resumes an interrupted simulation using existing grid files.
- `mocassinOutput`: Runs the output routines using the current grid files in the output subdirectory.
- `mocassinPlot`: Uses the current grid files in the output/ subdirectory to create 3d-emission maps.

## Building and Installation

Building the project requires MPI and Fortran compilers. On Debian/Ubuntu, install them via:
```bash
sudo apt-get install -y openmpi-bin libopenmpi-dev gfortran
```

The project is built using `make`:
```bash
make clean
make
```
Compiler flags can be overridden, for example: `make FCFLAGS="-Wall -Wextra"`.

## Running the Benchmarks

### Directory Setup
To run the models, the `data` and `dustData` directories must be copied or linked to `/usr/share/mocassin/` or available locally depending on your environment. Some missing benchmark data files may need manual symlinks (e.g., symlinking `forst-cryst-perp1-kyoto-uvx.nk` to `forsterite.nk` and copying `astroSilD2003.nk` to `astroSilD2003_base.nk`).

### Pure Photoionisation Benchmarks
You can try to run one or more of the available benchmark problems located in the `benchmarks/` directory.

To run the standard test problem (Meudon standard HII region), navigate to `benchmarks/test_Problem`:
```bash
cd benchmarks/test_Problem
mkdir -p output/
mpirun -np 1 ../../mocassin
```

Alternatively, to setup a benchmark manually:
1. Copy or link the benchmark `input.in` to your working `input` directory.
2. Ensure you have the corresponding abundance file.
3. Run using MPI: `mpirun -np <num_procs> ./mocassin`

### Pure Dust Benchmarks
There are also 1D and 2D benchmark models included under `benchmarks/dust/1D` and `benchmarks/dust/2D`. Copy the required input file to `input/input.in` and execute MOCASSIN.

## Input and Output Files

- **Input files**: Placed in the `input/` directory (e.g. `input.in`, abundance files, density distribution files).
- **Data files**: Located in the `data/` directory (atomic data files) and `dustData/` (dust optical data library).
- **Output files**: Saved in the `output/` directory (e.g. `ionratio.out`, `lineFlux.out`, `temperature.out`, `SED.out`, `grid0.out`, `grid1.out`, `grid2.out`, `grid3.out`, `dustGrid.out`).

See `man/mocassin.1` for a complete list of keywords that can be configured in the `input.in` file.

---

## Some Notes on Dust Species

| Grain Species | Type | Size Range (μm) | ρ (g/cm³) | Sublim. Temp (K) | Molecular Formula | Mol. Weight (amu) | ⟨amu/atom⟩ | Work Function (Ry) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Silicates (amorphous) | Mg/Fe silicates | 0.005–0.25 | 2.5–3.5 | 1200–1500 | Approx. MgFeSiO₄ | ~172 | ~21.5 | ~0.37 (5 eV) |
| Graphite | Carbonaceous | 0.005–0.25 | 2.2 | ~2000 | C | 12 | 12.0 | 0.33–0.37 (4.5–5 eV) |
| Olivine | (Mg,Fe)₂SiO₄ | 0.01–1.0 | 3.3–4.4 | 1400–1500 | Mg₂SiO₄ | 140.7 | 23.5 | ~0.37 (5 eV) |
| Forsterite | Mg₂SiO₄ (pure olivine) | 0.1–1.0 | 3.27 | 1400 | Mg₂SiO₄ | 140.7 | 23.5 | ~0.37 (5 eV) |
| Quartz | SiO₂ | 0.01–0.5 | 2.65 | 1200–1300 | SiO₂ | 60.08 | 20.0 | ~0.33 (4.5 eV) |
| Enstatite | MgSiO₃ | 0.1–1.0 | 3.2 | 1350–1400 | MgSiO₃ | 100.4 | 25.1 | ~0.37 (5 eV) |

Of course! Here's the **LaTeX version** of the table, formatted for use in an article or report (e.g., with `booktabs` for better spacing). Following the table, I’ve included **references and sources** for each set of values where applicable.

### References & Sources

Here are sources used to compile the physical parameters (need to double check these):

#### Density and Molecular Weights
- **Draine, B. T. (2003)**, *Annual Review of Astronomy and Astrophysics*, 41, 241  
- **Henning, T. (2010)**, *ARA&A*, 48, 21  
- **Lodders, K. (2003)**, *ApJ*, 591, 1220 (solar abundances)  
- **CRC Handbook of Chemistry and Physics** for molecular weights and densities

#### Sublimation Temperatures
- **Gail & Sedlmayr (1999)**, *A&A*, 347, 594  
- **Tielens, A. G. G. M. (2005)**, *The Physics and Chemistry of the Interstellar Medium*  
- **Duschl et al. (1996)**, *A&A*, 312, 624  
- Note: Temperatures vary with ambient pressure and gas density.

#### Work Function
- **Draine & Salpeter (1979)**, *ApJ*, 231, 77 (photoelectric emission from grains)  
- **Jenkins (2009)**, *ApJ*, 700, 1299 (dust charging)  
- **CRC Handbook** for graphite and silica values  
- For silicates, work functions are inferred from electron yield and ionization potentials (~5 eV = 0.37 Ry)
