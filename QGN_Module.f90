module QGN_Module
use ISO_C_BINDING
use iso_fortran_env, only: int64
use :: MKL_DFTI 
implicit none ! Applies to all contained subroutines
private ! All module vars and contained subroutines are default private
public dp, nx, ny, nz, pi, Initialize, TimeStep, WriteQP, WriteQPmodal, Cleanup, &
       WriteWB, WriteKE, GetTE

include "parameters.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module Variables:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop variables; dangerous to make these module variables... but convenient
integer :: i, j, k, kk
! output identifier(s)
integer :: outQ=30, outP=31, outPsi=32, outWR=33, outWL=34, outWN=35
integer :: outQm=35, outPm=36, outJm=37, outKE=38, outB=39, outWB=40
! Location of horizontal grid
real(dp) :: x(nx,ny), y(nx,ny)
! Vertical layer depths: H(nz) is the bottom layer; S=f^2/N^2(z)
real(dp) :: H(nz), S(0:nz)
! Zonal mean velocity profile, associated meridional PV and buoyancy gradients
real(dp) :: uBar(nz), qyBar(nz), byBar(nz-1) 
! Vertical modes, EVals = -kd^2
real(dp) :: Modes(nz,nz), EVals(nz)
! W Modes, W Evals
real(dp) :: WModes(nz-1,nz-1), WEVals(nz-1)
! MKL_DFTI plans
type(DFTI_DESCRIPTOR), pointer :: g2s_q, s2g_q, g2s_b, s2g_b, g2s_1, s2g_1
! Generic holders for FFT inputs
real(dp) :: grid(nx,ny,nz), gridE(nx*ny*nz)
Equivalence (grid, gridE)
real(dp) :: grid_b(nx,ny,nz+1), gridE_b(nx*ny*(nz+1))
Equivalence (grid_b, gridE_b)
real(dp) :: grid_1(nx,ny), gridE_1(nx*ny)
Equivalence (grid_1, gridE_1)
! Generic holders for FFT outputs
complex(dp) :: spec(nx/2 +1,ny,nz), specE((nx/2 +1)*ny*nz)
Equivalence (spec, specE)
complex(dp) :: spec_b(nx/2 +1,ny,nz+1), specE_b((nx/2 +1)*ny*(nz+1))
Equivalence (spec_b, specE_b)
complex(dp) :: spec_1(nx/2 +1,ny), specE_1((nx/2 +1)*ny)
Equivalence (spec_1, specE_1)
! Status for MKL DFTI calls
integer :: Status
! Wavenumbers
real(dp), dimension(nx/2 +1,ny) :: kx, ky, k2
! Some module variables; saves having to declare them in multiple subroutines
real(dp), dimension(nx,ny,nz) :: psi_phys, u_phys, v_phys, q_phys, jaco_phys
complex(dp), dimension(nx/2+1,ny,nz) :: psi_hat, u_hat, v_hat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module Functions & Subroutines:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize:
! Allocate aligned arrays and create FFTW plans
! Initialize threads
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Initialize(q_hat,N0)
    complex(dp), intent(out) :: q_hat(nx/2+1,ny,nz)
    integer,      intent(in) :: N0
    real(dp) :: qp(nx,ny,nz) 
    real(dp) :: dx = Lx / real(nx,dp), dy = Ly / real(ny, dp)
    real(dp) :: UBarBC(1:nz), UBarOC(1:nz)
    ! L matrix, LeftEVs
    real(dp) :: L(nz,nz),ModesL(nz,nz),L_copy(nz,nz)
    ! Matrices for W inversion
    real(dp) :: SD2S(nz-1,nz-1),WModesL(nz-1,nz-1)
    ! Imaginary part of eigenvalues
    real(dp) :: EValsI(nz)
    ! LAPACK workspace size
    integer, parameter :: LWORK = 34*nz
    ! LAPACK workspace
    real(dp) :: WORK(LWORK)
    ! Output for LAPACK
    integer :: INFO
    ! Some MKL vars 
    integer :: Dims(2)
    integer :: real_strides(3), cplx_strides(3)
    integer :: zeroInd, maxInd ! index of zero and 1st baroclinic eigenvalues
    ! Used to generate ubar
    real(dp) :: tmp(nz)
    ! Used to initialize q0
    real(dp) :: zc, qscale, q0_profile(nz)
    character(len=16) :: file_name

    print *,'---------------------------------------------------------------'
    print *,' dx = ', dx,' dy = ', dy

! Horizontal locations of grid points
    do j=1,ny,1
        y(:,j) = -(Ly+dy)/2._dp + dy*real(j,dp) 
    end do
    do i=1,nx,1
        x(i,:) = -(Lx+dx)/2._dp + dx*real(i,dp)
    end do

! Stratification and vertical grid
    open(unit=99,file='H.dat',access='STREAM',status='OLD')
    read(99) H 
    close(99)
    H = Htot*H/sum(H)
    open(unit=99,file='S.dat',access='STREAM',status='OLD')
    read(99) S
    close(99)

    print *,'---------------------------------------------------------------'
    print *,' Layer depths  '
    print *, H(:)
    print *,'---------------------------------------------------------------'
    print *,' f^2/N^2  '
    print *, S(:)
    print *,'---------------------------------------------------------------'

! Initialize vertical modes
    ! Set elements of L ([Del+L]psi=q)
    L = 0._dp
    L(1,2) = 2._dp*S(1)/(H(1)*(H(1)+H(2)))
    L(1,1) =-L(1,2)
  if( nz > 2 ) then
    do k=2,nz-1
        L(k,k-1) = 2._dp*S(k-1)/(H(k)*(H(k)+H(k-1)))
        L(k,k+1) = 2._dp*S(k)/(H(k)*(H(k)+H(k+1)))
        L(k,k)   =-(L(k,k-1)+L(k,k+1))
    end do
  end if
    L(nz,nz-1) = 2._dp*S(nz-1)/(H(nz)*(H(nz)+H(nz-1)))
    L(nz,nz)   =-L(nz,nz-1)

    open(unit=99,file='L.dat',access='STREAM',status='REPLACE')
    write(99) L
    close(99)

    ! Eigenvalue decomposition; DGEEV over-writes L_copy
    L_copy = L
    call DGEEV('N','V',nz,L_copy,nz,EVals,EValsI,ModesL,nz,Modes,nz,WORK,LWORK,INFO)

    ! Now normalize so that depth average of Modes(:,k)**2 is 1 and sign at top
    ! is positive
    do k=1,nz
        Modes(:,k) = sign(1._dp,Modes(1,k)) * Modes(:,k) / sqrt(sum(H*Modes(:,k)**2)/Htot)
    end do

    open(unit=77,file='modes.dat',access='STREAM',status='REPLACE')
    write(77) Modes
    close(77)

    print *,'---------------------------------------------------------------'
    print *,' If the following line is zero then vertical modes initialized '
    print *, INFO
    print *,'---------------------------------------------------------------'
    print *,' The eigenvalues are '
    print *, EVals(:)
    print *,'---------------------------------------------------------------'
    print *,' The deformation radii are '
    print *, 1._dp/SQRT(-EVals(:))
    print *,'---------------------------------------------------------------'

! Initialize W Modes
    ! Elements of SD2S
    ! Omega eqn is (-k^2N^2 + f0^2 D2)w = wRHS; rescale to
    !        (N/f0)(-k^2I + (f0/N) D2 (f0/N)) (Nw/f0) = wRHS/f0**2
    ! First construct D2
    SD2S(:,:) = 0._dp
    SD2S(1,2) = 2._dp/(H(2)*(H(1)+H(2)))
    SD2S(1,1) =-2._dp/(H(1)*H(2))
    do k=2,nz-2
        SD2S(k,k-1) = 2._dp/(H(k  )*(H(k)+H(k+1)))
        SD2S(k,k+1) = 2._dp/(H(k+1)*(H(k)+H(k+1)))
        SD2S(k,k)   =-2._dp/(H(k)*H(k+1))
    end do
    SD2S(nz-1,nz-2) = 2._dp/(H(nz-1)*(H(nz)+H(nz-1)))
    SD2S(nz-1,nz-1) =-2._dp/(H(nz-1)*H(nz))
    ! Now scale rows by S
    do k=1,nz-1
        SD2S(k,:) = S(k)*SD2S(k,:)
    end do
    ! Now scale columns by S
    do k=1,nz-1
        SD2S(:,k) = S(k)*SD2S(:,k)
    end do

    open(unit=99,file='SD2S.dat',access='STREAM',status='REPLACE')
    write(99) SD2S
    close(99)

    ! Get eigenvalue decomposition 
    call DGEEV('N','V',nz-1,SD2S,nz-1,WEVals,EValsI(1:nz-1),WModesL,nz-1,WModes,nz-1,WORK,LWORK,INFO)

    ! Note that D2 is self-adjoint with respect to the inner product
    ! defined by the diagonal scaling with H(k)+H(k+1). Diagonal
    ! scaling by S does not change this. So the eigenvectors of 
    ! SD2S should be orthogonal with respect to the same inner
    ! product.
    ! Now normalize 
    do k=1,nz-1
        tmp(1) = 0._dp
        do kk=1,nz-1
            tmp(1) = tmp(1) + (H(kk)+H(kk+1))*WModes(kk,k)**2
        end do
        WModes(:,k) = WModes(:,k) / sqrt(tmp(1))
    end do
    open(unit=99,file='WModes.dat',access='STREAM',status='REPLACE')
    write(99) WModes
    close(99)

    print *,'---------------------------------------------------------------'
    print *,' If the following line is zero then W  modes initialized '
    print *, INFO
    print *,'---------------------------------------------------------------'
    print *,' The W deformation radii are '
    print *, 1._dp/SQRT(-WEVals(:))
    print *,'---------------------------------------------------------------'


! Mean zonal velocity, PV and buoyancy gradient
    open(unit=99,file='UBar.dat',access='STREAM',status='OLD')
    read(99) uBar
    close(99)
    qyBar = matmul(L,-uBar)
    open(unit=99,file='Qy.dat',access='STREAM',status='REPLACE')
    write(99) qyBar
    close(99)
    do k=1,nz-1
        byBar(k) = -(2._dp * f0 / (H(k) + H(k+1)) ) * (uBar(k)-uBar(k-1))
    end do
    open(unit=99,file='By.dat',access='STREAM',status='REPLACE')
    write(99) byBar
    close(99)

! Initialize wavenumbers
    do i=1,nx/2+1
        kx(i,:) = (2._dp*pi/Lx) * real(i-1,dp)
    end do

    do j=1,ny/2+1
        ky(:,j) = (2._dp*pi/Ly) * real(j-1,dp)
    end do
    do j=ny/2+2,ny
        ky(:,j) = (2._dp*pi/Ly) * real(j-ny-1,dp)
    end do

    k2 = kx**2+ky**2
! Initialize q0 from file or from an initial condition
   if( N0 > 1 ) then
       write (file_name,'(A2,I0.9,A4)') "q.",N0,'.dat'
       print *,' Opening q0 file  '
       print *,  file_name
       open(unit=10,file=file_name,access='STREAM',status='OLD')
       read(unit=10) qp
       close(unit=10)
   else
! Set initial condition for freely decaying case: small m, large k
       qscale = 1.e-5
       open(unit=99,file='q0_profile.dat',access='STREAM',status='OLD')
       read(99) q0_profile
       close(99)
       zc = 0._dp
       do k=1,nz
           zc = zc - 0.5*H(k)/Htot
           do i=1,nx
               do j=1,ny
                   !qp(i,j,k) = qscale*(1._dp+tanh(28*(1._dp-zc/Htot - 0.95))) - 5.415876413396763e-06
                   !qp(i,j,k) = qscale*(exp(-(zc/0.2_dp)**2 / (0.5_dp - zc/0.16_dp))- 0.3_dp)
                   qp(i,j,k) = qscale * q0_profile(k) * cos(2._dp*pi*x(i,j)/Lx + pi/4._dp*sin(2._dp*pi*y(i,j)/Ly)**2)
               end do
           end do
           zc = zc - 0.5*H(k)/Htot
       end do
   end if


! Set up MKL DFTI plans
    include "MKL_DFTI_Setup.f90"

! Transform the initial condition
    print *,' Transforming the initial condition '
    grid = qp
    Status = DftiComputeForward(g2s_q, gridE, specE)
    print *, Status
    q_hat = spec
    call DeAlias(q_hat)

    print *,'---------------------------------------------------------------'
    print *,' Initialization complete '
    print *,'---------------------------------------------------------------'

end subroutine Initialize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DeAlias: zeros out top 1/3 of coefficients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DeAlias(field)
    complex(dp), intent(inout)  :: field(nx/2+1,ny,nz)
    do i=(nx/3+1),(nx/2+1),1
        field(i,:,:) = cmplx(0._dp,0._dp)
    end do
    do j=(ny/3+1),(2*(ny/3)+1),1
        field(:,j,:) = cmplx(0._dp,0._dp)
    end do
    field(1,1,:) = 0._dp

end subroutine DeAlias


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Jacobian
! Computes 2/3-rule dealiased jacobian, returns 
! Fourier coefficients thereof in jaco_spec
! Also computes spectral coefficients of v for beta
! and of psi for viscosity
! Assumes dealiased input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Jacobian(q_hat,jaco_hat)
    complex(dp), intent(in)  :: q_hat(nx/2+1,ny,nz)
    complex(dp), dimension(nx/2+1,ny,nz), intent(out) :: jaco_hat
    complex(dp), dimension(nx/2+1,ny,nz) :: uq_hat, vq_hat

    ! Get psi_hat, u_hat, v_hat
    call GetPsi(q_hat)
    do k=1,nz
        u_hat(:,:,k) =-cmplx(0._dp,1._dp)*ky*psi_hat(:,:,k)
        v_hat(:,:,k) = cmplx(0._dp,1._dp)*kx*psi_hat(:,:,k)
    end do
    ! Get u, v, q
    spec = q_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    q_phys = grid/real(nx*ny,dp)
    spec = u_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    u_phys = grid/real(nx*ny,dp)
    spec = v_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    v_phys = grid/real(nx*ny,dp)
    ! Get uq_hat, vq_hat
    grid = u_phys*q_phys
    Status = DftiComputeForward(g2s_q, gridE, specE)
    uq_hat = spec
    grid = v_phys*q_phys
    Status = DftiComputeForward(g2s_q, gridE, specE)
    vq_hat = spec
    ! Compute jacobian_spec
    do k=1,nz
        jaco_hat(:,:,k) = cmplx(0._dp,1._dp)*kx*uq_hat(:,:,k) &
                         +cmplx(0._dp,1._dp)*ky*vq_hat(:,:,k)
    end do

end subroutine Jacobian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! JacobianFG
! Computes 2/3-rule dealiased jacobian, returns 
! Fourier coefficients thereof in jaco_spec
! Same as Jacobian, except it computes J[F,G]
! for arbitrary F and G instead of J[psi,q]
! Assumes dealiased input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine JacobianFG(f_hat,g_hat,jaco_hat)
    complex(dp), dimension(nx/2+1,ny,nz), intent(in)  :: f_hat, g_hat
    complex(dp), dimension(nx/2+1,ny,nz), intent(out) :: jaco_hat
    complex(dp), dimension(nx/2+1,ny,nz) :: uf_hat, vf_hat
    complex(dp), dimension(nx/2+1,ny,nz) :: ug_hat, vg_hat
    real(dp), dimension(nx,ny,nz) :: g_phys, uf_phys, vf_phys

    ! Get uf_hat, vf_hat
    do k=1,nz
        uf_hat(:,:,k) =-cmplx(0._dp,1._dp)*ky*f_hat(:,:,k)
        vf_hat(:,:,k) = cmplx(0._dp,1._dp)*kx*f_hat(:,:,k)
    end do
    ! Get u, v, q
    spec = g_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    g_phys = grid/real(nx*ny,dp)
    spec = uf_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    uf_phys = grid/real(nx*ny,dp)
    spec = vf_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    vf_phys = grid/real(nx*ny,dp)
    ! Get uq_hat, vq_hat
    grid = uf_phys*g_phys
    Status = DftiComputeForward(g2s_q, gridE, specE)
    ug_hat = spec
    grid = vf_phys*g_phys
    Status = DftiComputeForward(g2s_q, gridE, specE)
    vg_hat = spec
    ! Compute jacobian_spec
    do k=1,nz
        jaco_hat(:,:,k) = cmplx(0._dp,1._dp)*kx*ug_hat(:,:,k) &
                         +cmplx(0._dp,1._dp)*ky*vg_hat(:,:,k)
    end do
    ! DeAlias output
    call DeAlias(jaco_hat)

end subroutine JacobianFG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetRHS: computes FFT of 
!-J[psi,q]-uBar*qx-(beta+qyBar)*v-nu*del[q] - Ekman
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetRHS(q_hat, RHS)
    complex(dp), intent(in)  :: q_hat(nx/2+1,ny,nz)
    complex(dp), intent(out) :: RHS(nx/2+1,ny,nz)
    complex(dp), dimension(nx/2+1,ny,nz) :: jaco_hat
    complex(dp), dimension(nx/2+1,ny) :: drag_hat
    real(dp), dimension(nx,ny) :: URMS
    real(dp), dimension(nx/2+1,ny) :: k8

    ! Advection
    call Jacobian(q_hat,jaco_hat)
    RHS = -jaco_hat
    ! Mean advection, beta
    do k=1,nz
        RHS(:,:,k) = RHS(:,:,k) - uBar(k)*cmplx(0._dp,1._dp)*kx*q_hat(:,:,k) &
                    - (beta + qyBar(k))*v_hat(:,:,k)
    end do
    ! (hyper)viscosity
    if (is_explicit) then
        k8 = k2**4
        do k=1,nz
            RHS(:,:,k) = RHS(:,:,k) - (A2*k2 + A8*k8)*q_hat(:,:,k)
        end do
    end if
    ! Ekman
    RHS(:,:,nz) = RHS(:,:,nz) + (r0*Htot/H(nz)) * k2 * psi_hat(:,:,nz)
    ! Quadratic drag
    drag_hat = (0._dp,0._dp)
    URMS = sqrt(u_phys(:,:,nz)**2+v_phys(:,:,nz)**2)
    grid_1 = URMS*u_phys(:,:,nz)
    Status = DftiComputeForward(g2s_1, gridE_1, specE_1)
    drag_hat = cmplx(0._dp,1._dp)*ky*spec_1
    grid_1 = URMS*v_phys(:,:,nz)
    Status = DftiComputeForward(g2s_1, gridE_1, specE_1)
    drag_hat = drag_hat - cmplx(0._dp,1._dp)*kx*spec_1
    RHS(:,:,nz) = RHS(:,:,nz) + (C_d*Htot/H(nz))*drag_hat
    ! Dealias
    call DeAlias(RHS)

end subroutine GetRHS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TimeStep: as it sounds.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TimeStep(q_hat,t,dt)
    complex(dp), intent(inout) :: q_hat(nx/2+1,ny,nz)
    real(dp),    intent(inout) :: t, dt

    include "ARK43_QGLeith.f90"
    call DeAlias(q_hat)

end subroutine TimeStep


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetPsi: Gets psi_hat from q_hat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetPsi(q_hat)
    complex(dp), dimension(nx/2+1,ny,nz), intent(in) :: q_hat
    complex(dp), dimension(nx/2+1,ny,nz) :: q_hat_mode, psi_hat_mode

    q_hat_mode(:,:,:) = (0._dp,0._dp)
    ! Get q_hat_mode and psi_hat_mode
!$OMP PARALLEL DO 
    do k=1,nz
        do kk=1,nz
            q_hat_mode(:,:,k) = q_hat_mode(:,:,k) + H(kk)*Modes(kk,k)*q_hat(:,:,kk)
        end do
        q_hat_mode(:,:,k) = q_hat_mode(:,:,k)/Htot
        psi_hat_mode(:,:,k) = q_hat_mode(:,:,k)/(-k2+EVals(k))
        psi_hat_mode(1,1,k) = (0._dp,0._dp)
    end do
!$OMP END PARALLEL DO
    ! Get psi_hat
    psi_hat = (0._dp,0._dp)
!$OMP PARALLEL DO 
    do k=1,nz
        do kk=1,nz
           psi_hat(:,:,k) = psi_hat(:,:,k) + psi_hat_mode(:,:,kk)*Modes(k,kk)
        end do
    end do
!$OMP END PARALLEL DO

end subroutine GetPsi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetTE: Returns total energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetTE(q_hat,TE)
    complex(dp), intent(in) :: q_hat(nx/2+1,ny,nz)
    real(dp),   intent(out) :: TE

    call GetPsi(q_hat)
    spec = psi_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    psi_phys = grid/real(nx*ny,dp)
    spec = q_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    q_phys = grid/real(nx*ny,dp)
    TE = 0._dp
    do k=1,nz
        TE = TE - H(k) * SUM(psi_phys(:,:,k) * q_phys(:,:,k))
    end do
    TE = TE / (Htot * real(2*nx*ny, dp))

end subroutine GetTE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteKE: Writes KE(z) to a direct-access record
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteKE(q_hat,RecNo)
    complex(dp), intent(in) :: q_hat(nx/2+1,ny,nz)
    integer,     intent(in) :: RecNo

    real(dp) :: KE(nz)
    !integer(kind=int64), save :: my_pos = 1

    call GetPsi(q_hat)
    do k=1,nz
        u_hat(:,:,k) =-cmplx(0._dp,1._dp)*ky*psi_hat(:,:,k)
        v_hat(:,:,k) = cmplx(0._dp,1._dp)*kx*psi_hat(:,:,k)
    end do
    ! Get u, v
    spec = u_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    u_phys = grid/real(nx*ny,dp)
    spec = v_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    v_phys = grid/real(nx*ny,dp)

    do k=1,nz
        KE(k) = SUM(u_phys(:,:,k)**2 + v_phys(:,:,k)**2)
    end do
    KE(:) = KE(:) / real(2*nx*ny,dp)

    ! Write KE 
    if( RecNo == 1 ) then
        open(unit=outKE,file="KE.dat",access='STREAM',status='REPLACE')
    else
        open(unit=outKE,file="KE.dat",access='STREAM',status='OLD',&
             position='append')
    end if
    write(outKE) KE
    close(outKE)

end subroutine WriteKE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteQP: Writes q_phys and psi_phys to
! direct-access records
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteQP(q_hat,nt)
    complex(dp), intent(in) :: q_hat(nx/2+1,ny,nz)
    integer,     intent(in) :: nt
    character(len=50) :: file_name

    ! Get psi_hat
    call GetPsi(q_hat)
    ! Get q_phys and psi_phys
    spec = q_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    q_phys = grid/real(nx*ny,dp)
    spec = psi_hat
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    psi_phys = grid/real(nx*ny,dp)
    ! Write fields

    write (file_name,'(A2,I0.9,A4)') "q.",nt,'.dat'
    open(unit=outQ,file=file_name,access='STREAM',status='UNKNOWN')
    ! if ( RecNo == 1) then
    !     REWIND outQ
    ! end if
    write(outQ) q_phys
    close(outQ)
    write (file_name,'(A2,I0.9,A4)') "p.",nt,'.dat'
    open(unit=outP,file=file_name,access='STREAM',status='UNKNOWN')
    ! if ( RecNo == 1) then
    !     REWIND outP
    ! end if
    write(outP) psi_phys
    close(outP)

end subroutine WriteQP

subroutine WriteQPmodal(q_hat,nt)
    complex(dp), intent(in) :: q_hat(nx/2+1,ny,nz)
    integer,     intent(in) :: nt
    character(len=50) :: file_name

    complex(dp), dimension(nx/2+1,ny,nz) :: q_hat_mode, psi_hat_mode
    complex(dp), dimension(nx/2+1,ny,nz) :: jaco_hat, jaco_hat_mode

    !call   Jacobian(q_hat,jaco_hat)

    q_hat_mode(:,:,:) = (0._dp,0._dp)
    jaco_hat_mode(:,:,:) = (0._dp,0._dp)
    ! Get q_hat_mode and psi_hat_mode
!$OMP PARALLEL DO 
    do k=1,nz
        do kk=1,nz
            q_hat_mode(:,:,k) = q_hat_mode(:,:,k) + H(kk)*Modes(kk,k)*q_hat(:,:,kk)
            !jaco_hat_mode(:,:,k) = jaco_hat_mode(:,:,k) + H(kk)*Modes(kk,k)*jaco_hat(:,:,kk)
        end do
        q_hat_mode(:,:,k) = q_hat_mode(:,:,k)/Htot
        !jaco_hat_mode(:,:,k) = jaco_hat_mode(:,:,k)/Htot
        psi_hat_mode(:,:,k) = q_hat_mode(:,:,k)/(-k2+EVals(k))
        psi_hat_mode(1,1,k) = (0._dp,0._dp)
    end do
!$OMP END PARALLEL DO

    ! Get q_phys and psi_phys
    spec = q_hat_mode
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    q_phys = grid/real(nx*ny,dp)
    spec = psi_hat_mode
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    psi_phys = grid/real(nx*ny,dp)
    !spec = jaco_hat_mode
    !Status = DftiComputeBackward(s2g_q, specE, gridE)
    !jaco_phys = grid/real(nx*ny,dp)
    ! Write fields
    write (file_name,'(A7,I0.9,A4)') "q_mode.",nt,'.dat'
    open(unit=outQm,file=file_name,access='STREAM',status='UNKNOWN')
    ! if ( RecNo == 1) then
    !     REWIND outQm
    ! end if
    write(outQm) q_phys
    close(outQm)
    write (file_name,'(A7,I0.9,A4)') "p_mode.",nt,'.dat'
    open(unit=outPm,file=file_name,access='STREAM',status='UNKNOWN')
    ! if ( RecNo == 1) then
    !     REWIND outPm
    ! end if
    write(outPm) psi_phys
    close(outPm)
    !write (file_name,'(A7,I0.9,A4)') "J_mode.",nt,'.dat'
    !open(unit=outJm,file=file_name,access='STREAM',status='UNKNOWN')
    !if ( RecNo == 1) then
    !    REWIND outJm
    !end if
    !write(outJm) jaco_phys
    !close(outJm)

end subroutine WriteQPmodal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetW: Solves [N^2nabla^2 + f0^2D^2]w = wRHS
! boundary values on w should already have
! been incorporated into wRHS, and will need
! to be put into w outside this routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetW(wRHS,w_phys)
    complex(dp), dimension(nx/2+1,ny,nz-1), intent(in) :: wRHS
    real(dp), dimension(nx,ny,0:nz),       intent(out) :: w_phys
    complex(dp), dimension(nx/2+1,ny,nz-1) :: rhs_hat_mode, w_hat_mode
    complex(dp), dimension(nx/2+1,ny,0:nz) :: w_hat

    rhs_hat_mode(:,:,:) = (0._dp,0._dp)
    ! Get rhs_hat_mode and w_hat_mode
!$OMP PARALLEL DO 
    do k=1,nz-1
        do kk=1,nz-1
            rhs_hat_mode(:,:,k) = rhs_hat_mode(:,:,k) + (H(kk)+H(kk+1))*WModes(kk,k)*S(kk)*wRHS(:,:,kk)
        end do
        rhs_hat_mode(:,:,k) = rhs_hat_mode(:,:,k)/f0**2
        w_hat_mode(:,:,k) = -rhs_hat_mode(:,:,k)/(k2-WEVals(k))
        w_hat_mode(1,1,k) = (0._dp,0._dp)
    end do
!$OMP END PARALLEL DO
    ! Get w_hat
    w_hat(:,:,:) = (0._dp,0._dp)
!$OMP PARALLEL DO 
    do k=1,nz-1
        do kk=1,nz-1
           w_hat(:,:,k) = w_hat(:,:,k) + S(k)*w_hat_mode(:,:,kk)*WModes(k,kk)
        end do
    end do
!$OMP END PARALLEL DO
    spec_b = w_hat
    Status = DftiComputeBackward(s2g_b, specE_b, gridE_b)
    w_phys = grid_b/real(nx*ny,dp)

end subroutine GetW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetB: Finds b at the top and bottom surfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetB(q_hat, b_hat_top, b_hat_bot)
    complex(dp), dimension(nx/2+1,ny,nz),         intent(in) :: q_hat  
    complex(dp), dimension(nx/2+1,ny), optional, intent(out) :: b_hat_top, b_hat_bot  
    real(dp), save :: c1, c2, c3
    logical, save  :: FirstCall = .TRUE.

    if (FirstCall) then
        c1 = 2._dp*(3._dp*H(1) + 2._dp*H(2) + H(3)) / &
             ( (H(1) + H(2)) * (H(1) + H(2) + H(3)) )
        c3 = 2._dp*(2._dp*H(1) + H(2)) / &
             ( (H(2) + H(3)) * (H(1) + H(2) + H(3)) )
        c2 = -(c1 + c3)
        FirstCall = .FALSE.
    end if

    call GetPsi(q_hat)
    if (present(b_hat_top)) b_hat_top(:,:) = f0*(c1*psi_hat(:,:,1) + c2*psi_hat(:,:,2) + c3*psi_hat(:,:,3))
    if (present(b_hat_bot)) b_hat_bot(:,:) = f0*(-c1*psi_hat(:,:,nz) - c2*psi_hat(:,:,nz-1) - c3*psi_hat(:,:,nz-2))

end subroutine GetB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteWB: Writes w_phys and b_phys to
! direct-access records
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteWB(q_hat, nt)
    complex(dp), intent(in) :: q_hat(nx/2+1,ny,nz)
    integer,     intent(in) :: nt 
    complex(dp), dimension(nx/2+1,ny,nz)   :: temp, jaco_hat
    complex(dp), dimension(nx/2+1,ny,nz-1) :: wRHS
    complex(dp), dimension(nx/2+1,ny,0:nz) :: b_hat
    real(dp), dimension(nx,ny,0:nz)  :: w_phys, b_phys
    real(dp), dimension(nx,ny)       :: URMS
    character(len=50) :: file_name

    ! Get  surface buoyancy; also updates psi_hat
    call GetB(q_hat,b_hat_top=b_hat(:,:,0),b_hat_bot=b_hat(:,:,nz))
    ! Get b_hat = f0*dpsi/dz 
    do k=1,nz-1
        b_hat(:,:,k)  = (2._dp*f0/(H(k)+H(k+1)))*(psi_hat(:,:,k) - psi_hat(:,:,k+1))
    end do
    spec_b = b_hat
    Status = DftiComputeBackward(s2g_b, specE_b, gridE_b)
    b_phys = grid_b/real(nx*ny,dp)
    ! Write b_phys
    write (file_name,'(A2,I0.9,A4)') "b.",nt,'.dat'
    open(unit=outB,file=file_name,access='STREAM',status='UNKNOWN')
    write(outB) b_phys
    close(outB)

    ! Get wRHS for bottom drag
    ! wRHS for drag should be -2 f0^2 w_n / H(n)*(H(n-1)+H(n))
    ! I think the Ekman layer has fw/H(n) = r*omega + c_d*curl(|u|u) 
    ! So wRHS in the bottom layer should be
    ! -(2f0/(H(n)+H(n-1)) * ( r0*Htot/H(n) * omega(n) + C_d*Htot/H(n) * curl(|u|u) )
    wRHS(:,:,:) = (0._dp,0._dp)
    ! Linear drag
    if( r0 > 0._dp ) wRHS(:,:,nz-1) = 2._dp*(r0*k2*Htot*f0 / (H(nz)*(H(nz)+H(nz-1))))*psi_hat(:,:,nz)
    ! quadratic drag
    if( C_d > 0._dp ) then
        u_hat(:,:,nz) =-cmplx(0._dp,1._dp)*ky*psi_hat(:,:,nz)
        v_hat(:,:,nz) = cmplx(0._dp,1._dp)*kx*psi_hat(:,:,nz)
        spec_1 = u_hat(:,:,nz)
        Status = DftiComputeBackward(s2g_1, specE_1, gridE_1)
        u_phys(:,:,nz) = grid_1/real(nx*ny,dp)
        spec_1 = v_hat(:,:,nz)
        Status = DftiComputeBackward(s2g_1, specE_1, gridE_1)
        v_phys(:,:,nz) = grid_1/real(nx*ny,dp)
        temp(:,:,nz) = (0._dp,0._dp)
        URMS(:,:) = sqrt(u_phys(:,:,nz)**2+v_phys(:,:,nz)**2)
        grid_1 = URMS*u_phys(:,:,nz)
        Status = DftiComputeForward(g2s_1, gridE_1, specE_1)
        temp(:,:,nz) = cmplx(0._dp,1._dp)*ky*spec_1
        grid_1 = URMS*v_phys(:,:,nz)
        Status = DftiComputeForward(g2s_1, gridE_1, specE_1)
        temp(:,:,nz) = temp(:,:,nz) - cmplx(0._dp,1._dp)*kx*spec_1
        wRHS(:,:,nz-1) = wRHS(:,:,nz-1) - 2._dp*(f0*C_d*Htot/(H(nz)*(H(nz)+H(nz-1))))*temp(:,:,nz)
    end if
    if( (r0 > 0._dp) .or. (C_d > 0._dp) ) then
        call GetW(wRHS,w_phys)
        spec_1 = (Htot/f0)*(-r0*k2*psi_hat(:,:,nz) + C_d*temp(:,:,nz))
        Status = DftiComputeBackward(s2g_1, specE_1, gridE_1)
        w_phys(:,:,nz) = grid_1/real(nx*ny,dp)
        ! Write w for bottom drag
        write (file_name,'(A7,I0.9,A4)') "w_drag.",nt,'.dat'
        open(unit=outWR,file=file_name,access='STREAM',status='UNKNOWN')
        write(outWR) w_phys
        close(outWR)
        if( nt == 1 ) then
            open(unit=outWB,file="WB_drag.dat",access='STREAM',status='REPLACE')
        else
            open(unit=outWB,file="WB_drag.dat",access='STREAM',status='OLD',&
                 position='append')
        end if
        write(outWB) SUM(SUM(w_phys*b_phys,1),1)/real(nx*ny,dp)
        close(outWB)
    end if
    

    ! Get wRHS for the linear terms
    if( sum(abs(uBar)) > 0._dp ) then
        ! First f0 d/dz [uBar dot grad omega + beta v]
        do k=1,nz
            temp(:,:,k) = (0._dp,1._dp)*kx*(-k2)*psi_hat(:,:,k)
            temp(:,:,k) = uBar(k)*temp(:,:,k) + beta * v_hat(:,:,k)
        end do
        do k=1,nz-1
            wRHS(:,:,k) = (2._dp*f0/(H(k)+H(k+1))) * (temp(:,:,k)-temp(:,:,k+1))
        end do
        ! Now -Del[uBar dot grad b + u dot grad bBar]
        ! Hard coded for equispaced
        do k=1,nz-1
            wRHS(:,:,k) = wRHS(:,:,k) &
                & + k2*(0.5_dp*(uBar(k)+uBar(k+1))*(0._dp,1._dp)*kx*b_hat(:,:,k) &
                & +     0.5_dp*(v_hat(:,:,k)+v_hat(:,:,k+1))*byBar(k) )
        end do
        call GetW(wRHS,w_phys)
        ! Write w for linear terms
        write (file_name,'(A9,I0.9,A4)') "w_linear.",nt,'.dat'
        open(unit=outWL,file=file_name,access='STREAM',status='UNKNOWN')
        write(outWL) w_phys
        close(outWL)
        if( nt == 1 ) then
            open(unit=outWB,file="WB_linear.dat",access='STREAM',status='REPLACE')
        else
            open(unit=outWB,file="WB_linear.dat",access='STREAM',status='OLD',&
                 position='append')
        end if
        write(outWB) SUM(SUM(w_phys*b_phys,1),1)/real(nx*ny,dp)
        close(outWB)
    end if

    ! Get wRHS for the nonlinear terms
    ! First f0 d/dz [ J[psi,omega] ]
    do k=1,nz
        temp(:,:,k) = -k2*psi_hat(:,:,k)
    end do
    call JacobianFG(psi_hat,temp,jaco_hat)
    do k=1,nz-1
        wRHS(:,:,k) = (2._dp*f0/(H(k)+H(k+1)))*(jaco_hat(:,:,k) - jaco_hat(:,:,k+1))
    end do
    ! Now -Del[ J[psi,b] ]
    do k=1,nz-1
        temp(:,:,k) = (H(k+1)*psi_hat(:,:,k)+H(k)*psi_hat(:,:,k+1))/(H(k)+H(k+1))
    end do
    call JacobianFG(temp,b_hat,jaco_hat)
    do k=1,nz-1
        wRHS(:,:,k) = wRHS(:,:,k) + k2*jaco_hat(:,:,k)
    end do
    call GetW(wRHS,w_phys)
    ! Write w for nonlinear terms
    write (file_name,'(A2,I0.9,A4)') "w.",nt,'.dat'
    open(unit=outWN,file=file_name,access='STREAM',status='UNKNOWN')
    write(outWN) w_phys
    close(outWN)
    if( nt == 1 ) then
        open(unit=outWB,file="WB.dat",access='STREAM',status='REPLACE')
    else
        open(unit=outWB,file="WB.dat",access='STREAM',status='OLD',&
             position='append')
    end if
    write(outWB) SUM(SUM(w_phys*b_phys,1),1)/real(nx*ny,dp)
    close(outWB)

end subroutine WriteWB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cleanup:
! De-allocate aligned arrays and destroy FFTW plans
! cleanup threads
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Cleanup

print *, '------------------------'
print *, ' Cleaning up '
    Status = DftiFreeDescriptor(g2s_q)
    print *, Status
    Status = DftiFreeDescriptor(s2g_q)
    print *, Status
print *, '------------------------'

end subroutine Cleanup

end module QGN_Module
