# Robust Approaches to Handling Complex Geometries with Galerkin Difference Methods

Matlab functions and scripts used in the paper _Robust Approaches to Handling
Complex Geometries with Galerkin Difference Methods_ by Jeremy E. Kozdon, Lucas
C. Wilcox, Thomas Hagstrom, and Jeffrey W. Banks.

Note: These are research codes written for the above paper, and it will take
some work in order to understand and apply these codes to other problems.

In general files in directory `src` _may_ be useful for other applications
whereas those in `acoustic2D` are specific files for the acoustic wave equation.

## Weight Adjusted Projection Test
The weight adjusted projection tests can be run using the script
`wa_error_test.m` in subdirectory `wa_error`

## Acoustic2D run scripts

### GD only scripts
The following scripts are in the directory `acoustic2D`

  * `GDConservation` :: GD conservation test
  * `GDConservation_central` :: GD conservation test with central flux
  * `GDConstPres` :: GD constant-state preservation test
  * `GDBox_curved` :: GD ghost basis curved box test
  * `GDBox_curved_extrapolation` :: GD extrapolation basis curved box test
  * `GDDisk2D` :: Disk with curved GD elements using ghost basis
  * `GDDisk2D_extrapolation` :: Disk with curved GD elements using extrapolation
  * `GDDisk2D_dt` :: Compute the maximum stable time step for the curved ghost basis GD disk mesh
  * `GDDisk2D_extrapolation_dt` :: Compute the maximum stable time step for the curved extrapolation GD disk mesh

### GD and simplicial scripts
The following scripts require the installation of the nodal DG codes from
[HesthavenWarburton2008](http://dx.doi.org/10.1007/978-0-387-72067-8); bibtex
reference below, and please cite if you use these codes! Installation
instructions for these codes are given below.

The following scripts are in the directory `acoustic2D`

  * `box_run` :: script to check error and compute time step for simple box problem (runs simplicial, GD ghost basis, and GD extrapolation)
  * `box_spectrum` :: script to compute the eigenvalue spectrum for simplicial, GD ghost basis, and GD extrapolation (data cached on git LFS)
  * `compute_max_timestep` :: Compute the maximum stable time step for the Taylor time stepping method for each of the schemes (uses spectrums from `box_spectrum`)
  * `DGDisk2D` :: Disk with curved simplicial elements
  * `DGDisk2D_dt` :: compute the maximum stable time step with curved simplicial elements
  * `GDDGDisk2D` :: Disk with coupled ghost basis GD and curved simplicial elements
  * `GDDGDisk2D_dt` :: compute the maximum stable time step with coupled ghost basis GD and curved simplicial elements ### Installing nodal-dg codes
  * `GDDGDisk2D_extrapolation` :: Disk with coupled extrapolation GD and curved simplicial elements
  * `GDDGDisk2D_dt` :: compute the maximum stable time step with coupled extrapolation GD and curved simplicial elements ### Installing nodal-dg codes
  * `GDDGInclusion` :: Run inclusion test problem
  * `GDDGInclusion_error_interp` :: Do the convergence analysis to `GDDGInclusion` runs (see note below)
  * `GDDGInclusion_geometry` :: Check to make sure the geometry interpolations are working

Note on `GDDGInclusion_error_interp`: The data needed by this script must be
generated prior to running. The data is stored in the directory
`acoustic2D/data/` using [GIT LFS](https://git-lfs.github.com/).

If you do not want to install GIT LFS, you can download the data directly with
either the following wget or curl. If you use this method you may wish to
locally untrack these files since they will make your respository larger if they
raw binary is added. This can be done with the commands:

```
git rm --cached -r acoustic2D/data/
git commit -m "untracking data files"
```

#### `wget` data download
```
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n5_r0_t15.mat -O acoustic2D/data/data_n5_r0_t15.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n5_r1_t15.mat -O acoustic2D/data/data_n5_r1_t15.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n5_r2_t15.mat -O acoustic2D/data/data_n5_r2_t15.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n7_r0_t15.mat -O acoustic2D/data/data_n7_r0_t15.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n7_r1_t15.mat -O acoustic2D/data/data_n7_r1_t15.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n7_r2_t15.mat -O acoustic2D/data/data_n7_r2_t15.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p1.mat -O acoustic2D/data/box_spectrum_p1.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p2.mat -O acoustic2D/data/box_spectrum_p2.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p3.mat -O acoustic2D/data/box_spectrum_p3.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p4.mat -O acoustic2D/data/box_spectrum_p4.mat
wget https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p5.mat -O acoustic2D/data/box_spectrum_p5.mat
```

#### `curl` data download
```
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n5_r0_t15.mat -o acoustic2D/data/data_n5_r0_t15.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n5_r1_t15.mat -o acoustic2D/data/data_n5_r1_t15.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n5_r2_t15.mat -o acoustic2D/data/data_n5_r2_t15.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n7_r0_t15.mat -o acoustic2D/data/data_n7_r0_t15.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n7_r1_t15.mat -o acoustic2D/data/data_n7_r1_t15.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/data_n7_r2_t15.mat -o acoustic2D/data/data_n7_r2_t15.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p1.mat -o acoustic2D/data/box_spectrum_p1.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p2.mat -o acoustic2D/data/box_spectrum_p2.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p3.mat -o acoustic2D/data/box_spectrum_p3.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p4.mat -o acoustic2D/data/box_spectrum_p4.mat
curl https://media.githubusercontent.com/media/bfam/GDComplexGeometries/master/acoustic2D/data/box_spectrum_p5.mat -o acoustic2D/data/box_spectrum_p5.mat
```

### Installing nodal-dg codes
As noted above, the nodal-dg codes from
[HesthavenWarburton2008](http://dx.doi.org/10.1007/978-0-387-72067-8) are used
to handle the simplicial domain. The codes in this repository assume that these
have have been installed in the base directory of this repository; if this is
not the case then you should set the variable `NODAL_DG_ROOT` in the file
`src/Globals2D_gdgd.m`

To install the nodal-dg code you can use the git command:

```
git clone https://github.com/tcew/nodal-dg
cd nodal-dg
git checkout 50c09d1 -b gddg
patch -p1 < ../nodal-dg.patch
```

The third command will ensure that you are on the commit that is known to work
with the codes and the fourth command patches the nodal-dg codes so that they
work with the gddg codes

# References:
```
@BOOK{HesthavenWarburton2008,
  title = {Nodal Discontinuous {G}alerkin Methods: {A}lgorithms, Analysis, and
           Applications},
  publisher = {Springer},
  year = {2008},
  author = {Hesthaven, Jan S. and Warburton, Tim},
  volume = {54},
  series = {Texts in Applied Mathematics},
  doi = {10.1007/978-0-387-72067-8}
}
```
