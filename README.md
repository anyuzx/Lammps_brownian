This is the custom `LAMMPS` fix for [overdamped langevin dynamics](https://en.wikipedia.org/wiki/Brownian_dynamics)(Brownian dynamics). `LAMMPS` has `fix langevin` for langevin dynamics simulation. However the algorithm LAMMPS use to integrate langevin dynamics equation is Velocity-Verlet. This is not suiable for high friction(overdamped) simulation. Here are several files I modified from original `fix_langevin.cpp` and some other files to give `LAMMPS` a fix for overdamped langevin dynamics simulation. This fix is only suitable for high friction case since velocity is overdamped. The detailed physics for this fix can be viewed in `derivation.pdf`. For low friction, please use original `fix langevin`. As for what value is approriate for high and low friction, one should do experiments themsevles.

## HOW TO USE

### Build LAMMPS

1. Put `fix_langevin_overdamp.cpp`, `fix_nve_overdamp.cpp` and their header files into `src` folder of `LAMMPS`. Put `fix_nve_overdamp_omp.cpp` and its header file in `src/USER-OMP` folder.
2. If you want to use `omp` acceraltion for this fix, do `package yes-user-omp` to install `omp` package
3. Build `LAMMPS`.

### Use it in input file

```
fix ID group-ID langevin/overdamp Tstart Tstop damp seed
```

* ID, group-ID
* lanevin/overdamp = style name of this fix command
* Tstart, Tstop = desired temperature at start/end of run (temperature units)
* damp = damping parameter (time units)
* seed = random number seed to use for white noise (positive integer)

Like `fix langevin`, `fix langevin/overdamp` also doesn't do integration. To use `fix langevin/overdamp` one also should include `fix nve/overdamp` in their input file before `fix langevin/overdamp`. The command for `fix nve/overdamp` is 

```
fix ID group-ID nve/overdamp damp
```

* damp = damping parameter (time units)

This damping parameter must be the same as the damping parameter specified in `fix langevin/overdamp`.

### Examples

```
fix 1 all nve/overdamp 0.01
fix 2 all langevin/overdamp 1.00 1.00 0.01 3301384
```

## IMPORTANT NOTES

* damping parameters in `fix nve/overdamp` and `fix langevin/overdamp` must have the same value
* 1/damp should be much larger than one. This is to ensure the parameters are indeed for overdamped/high friction/low intertia regime. For instance, damp = 0.01 is sufficiently large for brownian dynamics simulation. damp = 100 is too small to use this fix.
