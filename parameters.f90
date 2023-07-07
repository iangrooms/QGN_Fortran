!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameter statements:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! double precision, pi
    integer, parameter :: dp = C_DOUBLE !, dpc = C_DOUBLE_COMPLEX
    real(dp), parameter :: pi = 4._dp*atan(1.0_dp)

! Physical & grid parameters:
    ! Domain, meters
    real(dp), parameter :: Lx = 1.3824E6_dp, Ly = 1.3824E6_dp, Htot = 4800._dp
    ! Grid: Horizontal part must be divisible by 3.
    integer, parameter :: nx = (3*256)/2, ny = (3*256)/2, nz = 32 
    ! Coriolis coefficients
    real(dp), parameter :: f0 = 1.E-4_dp, beta = 2.E-11_dp
    ! Gravity, reference density
    real(dp), parameter :: g = 9.81_dp, rho0 = 1028.8_dp
    ! Viscosities & Ekman friction: r = r0*Htot/H(nz) fixes barotropic drag;r0=fdE/(2Htot)
    real(dp), parameter :: A2 = 0._dp, A8 = 1.E23_dp, r0 = 1.E-7_dp
    ! Quadratic drag parameter. C_d = c_d/h_BL where h_BL is the bottom
    ! boundary layer thickness. This is the drag felt by the barotropic mode,
    ! not by the bottom layer. C_d L_d is a measure of drag strength.
    real(dp), parameter :: C_d = 0.00001_dp

! Time step error tolerance; max-norm vorticity
    real(dp), parameter :: TOL = 1.E-10_dp

! Is PV diffusion treated explicitly or implicitly
    logical, parameter :: is_explicit = .FALSE.

! CFL for the fastest Rossby wave = 1 / (2 * frequency)
    real(dp), parameter :: dt_max =  pi / (beta * Lx)
