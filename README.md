:green_circle: **Latest version of LAMMPS now contains official fix command for overdamped BD**. See [here](https://docs.lammps.org/fix_brownian.html).

---

This is the custom `LAMMPS` fix for [overdamped Langevin dynamics](https://en.wikipedia.org/wiki/Brownian_dynamics)(Brownian dynamics). `LAMMPS` has `fix langevin` for Langevin dynamics simulation. However, the algorithm LAMMPS use to integrate Langevin dynamics equation is Velocity-Verlet. This is not suitable for high friction(overdamped) simulation. Here are several files I modified from original `fix_langevin.cpp` and some other files to give `LAMMPS` a fix for overdamped Langevin dynamics simulation. This fix is only suitable for high friction case since velocity is overdamped. For more information and derivation of equation of motions, you can read this blog post as well https://www.guangshi.io/posts/simulating-brownian/

## How to use

### Build LAMMPS

1. Put `fix_bd.cpp`/`fix_bd_baoab.cpp` and their header files into `src` folder of `LAMMPS`.
2. Build `LAMMPS`.

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

## Remarks

* For low friction, please use original `fix langevin`. As for what value is appropriate for high and low friction, one should do experiments themselves.

* damp should be much smaller than one. This is to ensure the parameters are indeed for overdamped/high friction/low inertia regime. For instance (with LJ unit for simplicity), damp = 0.01 is sufficiently small for Brownian dynamics simulation. damp = 100 is too large to use this fix. Instead, just use `fix langevin` for a large value of damp constant.

* :heavy_exclamation_mark: **Only fix bd and fix bd/baoab works. Do not use OMP versions.**

## Other Options

* `fix bd/baoab`: BAOAB algorithm for Brownian Dynamics simulation (see https://arxiv.org/abs/1203.5428).
* ~~`fix bd/srk`: SRK algorithm for Brownian Dynamics simulation.~~ (**discarded**)
* ~~`fix bd/omp`/`fix bd/srk/omp`/`fix bd/baoab/omp`: OMP version of Brownian Dynamics simulation.~~ (**discarded**)

#### Reference

* Leimkuhler, Benedict, and Charles Matthews. "Rational construction of stochastic numerical methods for molecular sampling." Applied Mathematics Research eXpress 2013.1 (2012): 34-56.
