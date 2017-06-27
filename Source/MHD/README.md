#Rouitnes for Magnetohydrodynamics 
All of the subroutines in this directory (except fort_advance_mhd) 
are called from subroutine fort_advance_mhd, 
which is called directly from Nyx::just_the_mhd in Nyx_mhd.cpp

In Nyx_advection_mhd_3d.f90:
  fort_advance_mhd
  advance_mhd_tile
  umeth3d
  ctoprim
  consup
  cmpflx
  hlld
  divu

In normalize_species_3d.f90:
  normalize_species_fluxes
  normalize_new_species

FORT_ADVANCE_MHD  --- (tiling)    ---> ADVANCE_MHD_TILE

ADVANCE_MHD_TILE  ---   ctoprim   ---  uflaten  (in flatten_3d.f90)

->                ---   umeth3d   ---  ppm*     (in mhd_ppm_3d.f90)
->                                ---  trace*   (in mhd_trace_ppm_3d.f90)
->                                ---  uslope   (in mhd_slope_3d.f90)
->                                ---  trans*   (in mhd_trans_3d.f90)
->                                ---  cmpflx  ---  riemann --- analriem

->                ---   divu 

->                ---   consup    ---  normalize_species_fluxes 

->                ---   enforce_minimum_density

->                ---   add_grav_source  (in add_grav_source_3d.f90)

->                ---   enforce_nonnegative_species (NOTE: FORT_ENFORCE_NONNEGATIVE_SPECIES is also called from the C++;

->                ---    FORT_ENFORCE_NONNEGATIVE_SPECIES then calls enforce_nonnegative_species;

->                ---     both of these routines are in Nyx/Source/Src_3d/enforce_nonnegative_species_3d.f90)

->                ---   normalize_new_species 
