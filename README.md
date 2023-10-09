# Regularized-Stokesian-Dynamics-with-LE-boundary-condition
===

Abstract: We extended the regularized Stokesian dynamics (SD) method to theoretically analyze the rheology of suspended particles based on periodic spherical arrays and the Lees-Edwards periodic boundary approach. The proposed regularized Stokesian dynamics method is applicable to suspensions from dilute to dense regime moving with low Reynolds numbers, and has been verified by comparing previous studies on the suspended spheres under constant shear. Regularized stokeslet can be used as an effective tool for numerical simulation of fine sediment transport.
<!-- ## Papar Information
- Title:  `paper name`
- Authors:  `A`,`B`,`C`
- Preprint: [https://arxiv.org/abs/xx]()
- Full-preprint: [paper position]()
- Video: [video position]() -->
## Code Details
This code was originally developed for calculating the rheological behavior of granular materials and has now been expanded to handle rigid granular assemblies(conglomerate.f90) as well as flexible plant stems(filament.f90), primarily based on quaternion dynamics.
## Install & Dependence
- gcc/g++
- gfortran
- openmp & mpi
- blas-lapack

## Installation
  ```fortran
git clone https://github.com/sherlockyzd/Regularized-Stokesian-Dynamics
sudo apt install gcc
sudo apt install g++
sudo apt install gfortran
sudo apt install libmpich-dev
sudo apt install libblas-dev
sudo apt install liblapack-dev
make clean
make
  ```



## Use
- for train
  ```
  chmod a+x runTest.sh
  ./runTest.sh
  ```


## Directory Hierarchy
```
|—— examples
|    |—— control_file.dat
|    |—— control_file0..dat
|    |—— DEM
|        |—— Results
|            |—— PartVisuForStokes.xmf
|            |—— PartVisuForStokes0000005000
|            |—— PartVisuForStokes0000010000
|            |—— Stokes.log
|            |—— WallsForStokes.backup
|    |—— DEM.prm
|    |—— DIM_property.txt
|    |—— fila_control_file.dat
|    |—— gmmres_output.VTF
|    |—— initial_bdy.dat
|    |—— initial_config.dat
|    |—— lubdat
|        |—— lubscalar.txt
|        |—— lubscalarmf.txt
|        |—— r2babc.dat
|        |—— r2babcd.dat
|        |—— r2bgh.dat
|        |—— r2bghd.dat
|        |—— r2bm.dat
|        |—— scalars_general_resistance_text_d.txt
|    |—— output.VTF
|    |—— output_bdy.VTF
|    |—— runSimulation
|    |—— yeta_mu.txt
|—— GRPerY.exe
|—— Makefile
|—— runTest.sh
|—— src
|    |—— conglomerate.f90
|    |—— csd_post.f90
|    |—— DEM
|        |—— DEM_gz.f90
|        |—— Prtcl_CL_and_CF.f90
|        |—— Prtcl_Comm.f90
|        |—— Prtcl_ContactForce_inc.f90
|        |—— Prtcl_ContactSearch.f90
|        |—— Prtcl_ContactSearchPW.f90
|        |—— Prtcl_Decomp_2d.f90
|        |—— Prtcl_DEMSystem (copy).f90
|        |—— Prtcl_DEMSystem.f90
|        |—— Prtcl_DumpPrtcl.f90
|        |—— Prtcl_Geometry (copy).f90
|        |—— Prtcl_Geometry.f90
|        |—— Prtcl_Hrchl_Munjiza.f90
|        |—— Prtcl_Integration.f90
|        |—— Prtcl_IOAndVisu.f90
|        |—— Prtcl_LogInfo.f90
|        |—— Prtcl_NBS_Munjiza.f90
|        |—— Prtcl_Parameters.f90
|        |—— Prtcl_Property.f90
|        |—— Prtcl_Timer.f90
|        |—— Prtcl_TypeDef.f90
|        |—— Prtcl_Variables.f90
|    |—— Ewald_summation.f90
|    |—— filament.f90
|    |—— GRPY_YZD.f90
|    |—— hydro_lub.f90
|    |—— initial_conf.f90
|    |—— lubrication tabel
|        |—— FarfieldForce
|            |—— Farfield.f90
|            |—— Farfield.o
|            |—— lambda.dat
|            |—— lapblas_double_excerpts.f
|            |—— lapblas_double_excerpts.o
|            |—— mainf
|            |—— mainf.f90
|            |—— run.sh
|            |—— tablefar.txt
|            |—— tensors.mod
|            |—— tools.f90
|            |—— tools.o
|        |—— hydro_lub_AK.f90
|        |—— hydro_lub_AK0.f90
|        |—— hydro_lub_YZD.f90
|        |—— hydro_lub_YZD1.f90
|        |—— main_gz.f90
|        |—— near-mid
|            |—— AX.dat
|            |—— AY.dat
|            |—— BY.dat
|            |—— Calc_AX.m
|            |—— Calc_AY.m
|            |—— Calc_BY.m
|            |—— Calc_CX.m
|            |—— Calc_CY.m
|            |—— Calc_GX.m
|            |—— Calc_GY.m
|            |—— Calc_HY.m
|            |—— Calc_MX.m
|            |—— Calc_MY.m
|            |—— Calc_MZ.m
|            |—— CX.dat
|            |—— CY.dat
|            |—— GX.dat
|            |—— GY.dat
|            |—— HY.dat
|            |—— MX.dat
|            |—— MY.dat
|            |—— MZ.dat
|        |—— tabel
|            |—— AX_tabel.m
|            |—— AX_tb.mat
|            |—— AX_tbmf.mat
|            |—— AY_tabel.m
|            |—— AY_tb.mat
|            |—— AY_tbmf.mat
|            |—— BY_tabel.m
|            |—— BY_tb.mat
|            |—— BY_tbmf.mat
|            |—— CX_tabel.m
|            |—— CX_tb.mat
|            |—— CX_tbmf.mat
|            |—— CY_tabel.m
|            |—— CY_tb.mat
|            |—— CY_tbmf.mat
|            |—— FarfieldScalar.mat
|            |—— GX_tabel.m
|            |—— GX_tb.mat
|            |—— GX_tbmf.mat
|            |—— GY_tabel.m
|            |—— GY_tb.mat
|            |—— GY_tbmf.mat
|            |—— HY_tabel.m
|            |—— HY_tb.mat
|            |—— HY_tbmf.mat
|            |—— lub0.m
|            |—— lubscalar.txt
|            |—— lubscalarmf.txt
|            |—— lubtable.m
|            |—— lubtablemf.m
|            |—— MX_tabel.m
|            |—— MX_tb.mat
|            |—— MX_tbmf.mat
|            |—— MY_tabel.m
|            |—— MY_tb.mat
|            |—— MY_tbmf.mat
|            |—— MZ_tabel.m
|            |—— MZ_tb.mat
|            |—— MZ_tbmf.mat
|    |—— main.f90
|    |—— MakefileRSD
|    |—— matrices.f90
|    |—— mymodules.f90
|    |—— NBS
|        |—— NBSPrtclTypeDef.f90
|        |—— PrtclFill.f90
|        |—— PrtclNBS_Munjiza.f90
|    |—— per_GRPY_YZD.f90
|    |—— Regularization.f90
|    |—— source_force.f90
|    |—— stepper.f90
|    |—— stokesian_Green.f90
|    |—— tensors.f90
|    |—— tools
|        |—— abstract_interfaces.f90
|        |—— algebra_lin.f90
|        |—— basic_tools.f90
|        |—— cg_matrix_omp.f90
|        |—— cg_mpi.f90
|        |—— cg_omp.f90
|        |—— cg_rc.f90
|        |—— gmres_mod.f90
|        |—— gmres_mod_gmres.f90
|        |—— lapblas_double_excerpts.f
|        |—— R1d_mod.f90
|        |—— real_type.f90
|    |—— traj.f90
|    |—— Wall_interaction.f90
