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
