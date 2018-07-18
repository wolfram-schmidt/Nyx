# AMREX_HOME defines the directory in which we will find all the BoxLib code
# AMREX_HOME = /project/projectdirs/nyx/src/amrex
 AMREX_HOME ?= /home/sireeves/amrex
# AMREX_HOME ?= /home/steven/amrex


# TOP defines the directory in which we will find Source, Exec, etc
TOP = /home/sireeves/Nyx
#TOP = ../../..

# compilation options
COMP    = gcc

USE_MPI = TRUE 
USE_OMP = FALSE 

PRECISION = DOUBLE
DEBUG     = TRUE
#DEBUG     = FALSE

# physics
DIM      = 3
USE_MHD = TRUE
USE_GRAV = FALSE

Bpack := ./Make.package
Blocs := .

include $(TOP)/Exec/Make.Nyx
