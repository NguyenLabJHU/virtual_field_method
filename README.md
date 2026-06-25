# VFM_JHU
Ocular Biomechanics: Virtual Fields Method (VFM) Optimization
This repository contains the MATLAB implementation of the Virtual Fields Method (VFM) applied to ocular tissues. The project focuses on identifying constitutive material properties using finite element modeling and experimental data, specifically accounting for prestress conditions within the ocular structure.

📁 Repository Structure
The project is organized into a modular hierarchy to facilitate local development and high-performance computing (HPC):

src/: Core MATLAB scripts and functions, including optimization loops, stress computation routines, and prestress integration.

data/: Input files, such as FEBio model files (.feb) containing geometry and fiber definitions, and simulation log files. (Note: Large data files are ignored by Git to maintain repository efficiency).

results/: Output directory for optimization results and generated .mat files.

VF/: Cached virtual work components and deformation gradients stored as CSV and MAT files to optimize re-computation times.

resources/: MATLAB project metadata and environment settings.

🚀 Getting Started
Open Project: Open the VFM_Source_Git.prj file in MATLAB to automatically configure search paths and environment variables.

Main Script: Run src/Run_test_uniform_optimization_multistart_calc_Fpre.m to execute the optimization process.

🛠 Technical Details
Solver Compatibility: Designed for seamless integration with FEBio.

Prestress Modeling: The identification process explicitly incorporates prestress effects.

Optimization: Utilizes optimization strategies to identify material parameters. The effect of prestress is updated after the solver reaches a local minimum.

Virtual Fields: Implements integration-point-based virtual work pertubation.
