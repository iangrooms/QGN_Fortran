# QGN: An N-Layer QG code in Fortran
Periodic domain size and other physical parameters are set in parameters.f90. Units are MKS. The number of time steps is set in QGN_Driver.

Background stratification, relative layer thicknesses, and mean flow profile are set in `set_background_profiles.ipynb`. You can alternatively just provide them as binary files.

There are four time-integration options.
Three explicit schemes are Adams-Bashforth 3, an adaptive RK3/2 pair (Bogacki & Shampine 1989), and RK4. A fourth-order IMEX Runge-Kutta scheme is available (ARK4(3)6L[2]SA -- Kennedy & Carpenter 2003) that treats PV diffusion implicitly and either uses a fixed hyperviscous coefficient (`ARK43.f90`) or uses a QG-Leith scaled hyperviscous coefficient (`ARK43_QGLeith.f90`). You pick the time integrator using an "include" statement in QGN_Module's TimeStep subroutine. AB3 and RK use a fixed time step. BSRK32 and ARK43 (can) use a PI.3.4 adaptive time step (Soderlind 2002).

The current version compiles and runs without errors on Cheyenne (NCAR) using intel + MKL with the following compile and environment commands:
1. `export OMP_NUM_THREADS=whatever`
2. `export OMP_STACKSIZE=whatever`
3. `ulimit -s unlimited`
4. `ifort -o QGN.x mkl_dfti.f90 QGN_Module.f90 QGN_Driver.f90 -assume byterecl -mkl=parallel -qopenmp -parallel -O3 -xHost -mcmodel medium -shared-intel`

At NCAR, after loading intel and mkl modules, you need to copy or link `/glade/u/apps/opt/intel/2019u5/mkl/include/mkl_dfti.f90` into the current directory. (location may change if you're using a different version of mkl.)

The code was first described in [Grooms (2016)](https://doi.org/10.1016/j.ocemod.2016.09.005). SD Bachman ported the code to [Chapel](https://github.com/sdbachman/Chapel_QGN). Rachel Robey contributed to the jupyter notebook that generates the background and added improved restarts.
