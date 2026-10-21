# gnet2gtc

`gnet2gtc` computes GTC markers from a GNET 5D distribution. It reads the
GNET HDF5 file, interpolates the distribution onto random points in the
GTC phase space, keeps points by rejection sampling, and writes the markers to
an HDF5 file. 

## Requirements

- GCC with C++20 support (tested 12.3)
- MPI with `mpicc`/`mpicxx` (tested MPICH 4.1.1)
- CUDA with `nvcc` (tested 12.1.1)
- CMake >= 3.20 (tested 3.26)
- HDF5 C library (tested 1.14.1)
- git

Load these first, for example with modules (names differ per system):

```bash
module load gcc/12.3.0 mpich/4.1.1 cuda/12.1.1 cmake/3.26.3 hdf5/1.14.1
```

## Build

The scripts in `build/` build all other dependencies from source and then
pcms with `gnet2gtc`, installing everything into one directory.

1. Get the scripts:

   ```bash
   git clone -b gnet-gtc-coupling https://github.com/SCOREC/pcms.git
   cd pcms/tools/gnet/build
   ```

2. Edit `env.sh`: set `GPU_ARCH` for your GPU (see the table below) and,
   if needed, `GNET_ROOT` (default `$HOME/gnet2gtc`).

3. Build:

   ```bash
   source env.sh
   ./build-deps.sh
   ./build-pcms.sh
   ```

The executable is installed at `$GNET_ROOT/install/bin/gnet2gtc`. See
[Usage](#usage) to run it.

Directory layout under `$GNET_ROOT`:

| Directory  | Contents                     |
|------------|------------------------------|
| `src/`     | git checkouts                |
| `build/`   | build trees                  |
| `install/` | installed libraries and `bin/gnet2gtc` |

### GPU architecture

| `GPU_ARCH` | GPUs                         |
|------------|------------------------------|
| `VOLTA70`  | V100                         |
| `AMPERE80` | A100                         |
| `AMPERE86` | RTX 30xx, A40, A10           |
| `ADA89`    | RTX 40xx, L4, L40            |
| `HOPPER90` | H100, H200                   |

Kokkos is built for one architecture. If you change `GPU_ARCH`, delete
`$GNET_ROOT/build` and `$GNET_ROOT/install` and rebuild.

## Usage

In a new shell, load the same modules as for the build, then:

```bash
source /path/to/pcms/tools/gnet/build/env.sh   # puts gnet2gtc on PATH
gnet2gtc --input gnet.h5 --output gtc.h5
```

Or run it by its full path (default `GNET_ROOT` shown):

```bash
~/gnet2gtc/install/bin/gnet2gtc --input gnet.h5 --output gtc.h5
```

With the defaults it writes 30 million markers.
Use `--samples` for fewer or more.

| Option                  | Default    | Description                                  |
|-------------------------|------------|----------------------------------------------|
| `--input <file>`        | required   | GNET HDF5 file                               |
| `--output <file>`       | required   | GTC output HDF5 file (overwritten)           |
| `--samples <n>`         | 30000000   | number of GTC markers (max 400000000)        |
| `--batch <n>`           | 1000000    | candidate points per sampling round          |
| `--seed <n>`            | 14000      | random seed                                  |
| `--degree <n>`          | 1          | MLS polynomial degree                        |
| `--lambda <x>`          | 1e-4       | MLS regularization                           |
| `--threshold-scale <x>` | 1e10       | threshold = scale / max(`volume3d`)          |
| `--verbose`             | off        | print some source points                     |
| `--help`                |            | show the options                             |


## Input

A GNET HDF5 file with:

| Dataset                                   | Type   | Description                     |
|-------------------------------------------|--------|---------------------------------|
| `/distribution5D/dimensions/mxrho`        | int    | radial bins                     |
| `/distribution5D/dimensions/mxzeta`       | int    | toroidal bins                   |
| `/distribution5D/dimensions/mxenergy`     | int    | energy bins                     |
| `/distribution5D/dimensions/mxpitch`      | int    | pitch bins                      |
| `/distribution5D/range/engmax_eV`         | double | maximum energy (eV)             |
| `/distribution5D/range/zetamax_rad`       | double | maximum toroidal angle (rad)    |
| `/distribution5D/data/f5d`                | double | 3D: (spatial cell, energy, pitch) |
| `/distribution5D/data/volume3d`           | double | 1D: volume of each spatial cell |

## Output

One dataset, `/gtc_n6d`: `N x 6` doubles, one row per marker
(`N` = `--samples`).

| Column | Quantity | Range        |
|--------|----------|--------------|
| 0      | rho      | [0, 1]       |
| 1      | theta    | [0, 2π]      |
| 2      | zeta     | [0, 2π/5]    |
| 3      | energy   | [2000, 200000] eV |
| 4      | pitch    | [-1, 1]      |
| 5      | value    | interpolated GNET distribution at the marker |

## Tested versions

`build-deps.sh` checks out these versions:

| Package        | Version               |
|----------------|-----------------------|
| Kokkos         | `4764a9a79` (develop, 4.5.99) |
| Kokkos Kernels | `4.6.01`              |
| Omega_h        | `3cd28bcd` (SCOREC)   |
| ADIOS2         | `v2.10.2`             |
| redev          | `1e4c876`             |
| perfstubs      | `2676497`             |
| Catch2         | `v3.8.0` (tests only) |
