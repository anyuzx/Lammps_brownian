This is the custom `LAMMPS` fix for [overdamped langevin dynamics](https://en.wikipedia.org/wiki/Brownian_dynamics)(Brownian dynamics). `LAMMPS` has `fix langevin` for langevin dynamics simulation. However the algorithm LAMMPS use to integrate langevin dynamics equation is Velocity-Verlet. This is not suiable for high friction(overdamped) simulation. Here are several files I modified from original `fix_langevin.cpp` and some other files to give `LAMMPS` a fix for overdamped langevin dynamics simulation. This fix is only suitable for high friction case since velocity is overdamped. The detailed physics for this fix can be viewed in `derivation.pdf`. For low friction, please use original `fix langevin`. As for what value is approriate for high and low friction, one should do experiments themsevles.

## HOW TO USE

### Build LAMMPS

1. Put `fix_\*.cpp` and their header files into `src` folder of `LAMMPS`. Put `fix_\*_omp.cpp` and its header file in `src/USER-OMP` folder.
2. If you want to use `omp` acceraltion for this fix, do `package yes-user-omp` to install `omp` package
3. Build `LAMMPS`.

### Use it in input file

```
fix ID group-ID bd temperature damp seed
```

* ID, group-ID
* bd = style name of this fix command
* temperature = desired temperature of run (temperature units)
* damp = damping parameter (time units)
* seed = random number seed to use for white noise (positive integer)

### Examples

```
fix 1 all bd 1.00 0.01 324231
```

## IMPORTANT NOTES

* 1/damp should be much larger than one. This is to ensure the parameters are indeed for overdamped/high friction/low intertia regime. For instance, damp = 0.01 is sufficiently large for brownian dynamics simulation. damp = 100 is too small to use this fix.

## Other Options

* `fix bd/baoab`: BAOAB algorithm for Brownian Dynamics simulation.
* `fix bd/srk`: SRK algorithm for Brownian Dynamics simulation.
* `fix bd/omp`: OMP version of Brownian Dynamics simulation.
