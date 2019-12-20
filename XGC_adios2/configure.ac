# -*- Autoconf -*-
# Process this file with autoconf 2.61 or newer to produce configure script.
AC_PREREQ([2.61])

# Set package name and version number and give contact e-mail address
AC_INIT([XGC-1],[3.0],[sku@pppl.gov])

# Set location of auxiliary autoconf files
AC_CONFIG_AUX_DIR([../acfiles])

# Specify config header file
AC_CONFIG_HEADERS([config.h])

# Set prefix to current directory by default
AC_PREFIX_DEFAULT([`pwd`])

# CPES default settings 
AX_CPES_DEFAULTS

# Initialize automake version 1.10 or newer
AM_INIT_AUTOMAKE([1.10])

# Option for maintainer mode
AM_MAINTAINER_MODE

# XGC-1 configure message
AC_MSG_NOTICE([Configuring $PACKAGE_STRING for $build])

# Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_SUBST(INSTALL,"$INSTALL -p")
AC_PROG_RANLIB

# Checks for Fortran 90 compiler.
AX_CPES_FC_OPTIONS

# Checks for C compiler (used to build binpack and camtimers libraries).
AX_CPES_CC_OPTIONS
# Need AM_PROG_CC_C_O macro for building binpack with per-target flags.
AM_PROG_CC_C_O

# PETSc option
AX_CPES_WITH_PETSC

# Compile target options
AX_CPES_ENABLE_OPT
AC_LANG_PUSH([Fortran])
AX_CPES_ENABLE_MPI
AC_LANG_POP
AX_CPES_ENABLE_64BIT
AX_CPES_ENABLE_OMP

# Optional packages
AX_CPES_WITH_BINPACK
AX_CPES_WITH_CAMTIMERS

# Third-party packages
AX_CPES_WITH_NETCDF
AX_CPES_WITH_PSPLINE
AX_CPES_WITH_IMSL

# Fix CPP defines option for IBM Fortran compiler
AC_CONFIG_COMMANDS_PRE([
  case "$FC" in
    *xlf*) # convert this DEFS flag into separate FCFLAGS and CFLAGS lists
      FCFLAGS="-WF,$DEFS $FCFLAGS"
      CFLAGS="$DEFS $CFLAGS"
      DEFS=""
    ;;
  esac])

# List of Makefiles to produce
AC_CONFIG_FILES([Makefile])

# Output files
AC_OUTPUT

# Final message
echo 
echo Now change directory to $ax_cpes_defaults_bldtarget and run GNU make to build XGC-1

