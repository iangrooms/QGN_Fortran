module QGN_Module
use ISO_C_BINDING
use iso_fortran_env, only: int64
use :: MKL_DFTI 
implicit none ! Applies to all contained subroutines
private ! All module vars and contained subroutines are default private
public dp, nx, ny, nz, pi, Initialize, TimeStep, WriteQP, WriteQPmodal, Cleanup

include "parameters.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module Variables:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop variables; dangerous to make these module variables... but convenient
integer :: i, j, k, kk
! output identifier(s)
integer :: outQ=30, outP=31, outPsi=32
integer :: outQm=35, outPm=36, outJm=37
! Location of horizontal grid
real(dp) :: x(nx,ny), y(nx,ny)
! Vertical layer depths: H(nz) is the bottom layer; S=f^2/N^2(z)
real(dp) :: H(nz), S(0:nz)
! Zonal mean velocity profile, associated meridional PV gradient
real(dp) :: uBar(nz), qyBar(nz)
! Vertical modes, EVals = -kd^2
real(dp) :: Modes(nz,nz), EVals(nz)
! MKL_DFTI plans
type(DFTI_DESCRIPTOR), pointer :: g2s_q, s2g_q, g2s_b, s2g_b
! Generic holders for FFT inputs
real(dp) :: grid(nx,ny,nz), gridE(nx*ny*nz)
Equivalence (grid, gridE)
! Generic holders for FFT outputs
complex(dp) :: spec(nx/2 +1,ny,nz), specE((nx/2 +1)*ny*nz)
Equivalence (spec, specE)
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
subroutine Initialize(qp,q_hat)
    real(dp), intent(in)     :: qp(nx,ny,nz)
    complex(dp), intent(out) :: q_hat(nx/2+1,ny,nz)
    real(dp) :: dx = Lx / real(nx,dp), dy = Ly / real(ny, dp)
    ! L matrix, LeftEVs
    real(dp) :: L(nz,nz),ModesL(nz,nz),L_copy(nz,nz)
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
    ! include "Stratification.f90"
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

! Mean zonal velocity, PV and buoyancy gradient
    open(unit=99,file='UBar.dat',access='STREAM',status='OLD')
    read(99) uBar
    close(99)
    qyBar = matmul(L,-uBar)
    print *,'---------------------------------------------------------------'
    print *,' ubar  '
    print *, uBar(:)
    print *,'---------------------------------------------------------------'
    print *,' qybar  '
    print *, qyBar(:)
    print *,'---------------------------------------------------------------'
    open(unit=99,file='Qy.dat',access='STREAM',status='REPLACE')
    write(99) qyBar
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

! Set up MKL DFTI plans
    include "MKL_DFTI_Setup.f90"

! Transform the initial condition
    print *,' Transforming the initial condition '
    grid = qp
    Status = DftiComputeForward(g2s_q, gridE, specE)
    print *, Status
    q_hat = spec
    call DeAlias(q_hat)

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetRHS: computes FFT of 
!-J[psi,q]-uBar*qx-(beta+qyBar)*v+diffusion+drag
! when is_explicit is false, does not compute
! diffusion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetRHS(q_hat, RHS)
    complex(dp), intent(in)  :: q_hat(nx/2+1,ny,nz)
    complex(dp), intent(out) :: RHS(nx/2+1,ny,nz)
    complex(dp), dimension(nx/2+1,ny,nz) :: jaco_hat
    complex(dp), dimension(nx/2+1,ny) :: drag_hat
    real(dp), dimension(nx,ny,nz) :: URMS
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
    URMS = sqrt(u_phys**2+v_phys**2)
    grid = URMS*u_phys
    Status = DftiComputeForward(g2s_q, gridE, specE)
    drag_hat = cmplx(0._dp,1._dp)*ky*spec(:,:,nz)
    grid = URMS*v_phys
    Status = DftiComputeForward(g2s_q, gridE, specE)
    drag_hat = drag_hat - cmplx(0._dp,1._dp)*kx*spec(:,:,nz)
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
! WriteQP: Writes q_phys and psi_phys to
! direct-access records
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteQP(q_hat,RecNo,nt)
    complex(dp), intent(in) :: q_hat(nx/2+1,ny,nz)
    integer,     intent(in) :: RecNo, nt
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
    if ( RecNo == 1) then
        REWIND outQ
    end if
    write(outQ) q_phys
    close(outQ)
    write (file_name,'(A2,I0.9,A4)') "p.",nt,'.dat'
    open(unit=outP,file=file_name,access='STREAM',status='UNKNOWN')
    if ( RecNo == 1) then
        REWIND outP
    end if
    write(outP) psi_phys
    close(outP)

end subroutine WriteQP

subroutine WriteQPmodal(q_hat,RecNo,nt)
    complex(dp), intent(in) :: q_hat(nx/2+1,ny,nz)
    integer,     intent(in) :: RecNo, nt
    character(len=50) :: file_name

    complex(dp), dimension(nx/2+1,ny,nz) :: q_hat_mode, psi_hat_mode
    complex(dp), dimension(nx/2+1,ny,nz) :: jaco_hat, jaco_hat_mode

    call   Jacobian(q_hat,jaco_hat)

    q_hat_mode(:,:,:) = (0._dp,0._dp)
    jaco_hat_mode(:,:,:) = (0._dp,0._dp)
    ! Get q_hat_mode and psi_hat_mode
!$OMP PARALLEL DO 
    do k=1,nz
        do kk=1,nz
            q_hat_mode(:,:,k) = q_hat_mode(:,:,k) + H(kk)*Modes(kk,k)*q_hat(:,:,kk)
            jaco_hat_mode(:,:,k) = jaco_hat_mode(:,:,k) + H(kk)*Modes(kk,k)*jaco_hat(:,:,kk)
        end do
        q_hat_mode(:,:,k) = q_hat_mode(:,:,k)/Htot
        jaco_hat_mode(:,:,k) = jaco_hat_mode(:,:,k)/Htot
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
    spec = jaco_hat_mode
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    jaco_phys = grid/real(nx*ny,dp)
    ! Write fields
    write (file_name,'(A7,I0.9,A4)') "q_mode.",nt,'.dat'
    open(unit=outQm,file=file_name,access='STREAM',status='UNKNOWN')
    if ( RecNo == 1) then
        REWIND outQm
    end if
    write(outQm) q_phys
    close(outQm)
    write (file_name,'(A7,I0.9,A4)') "p_mode.",nt,'.dat'
    open(unit=outPm,file=file_name,access='STREAM',status='UNKNOWN')
    if ( RecNo == 1) then
        REWIND outPm
    end if
    write(outPm) psi_phys
    close(outPm)
    write (file_name,'(A7,I0.9,A4)') "J_mode.",nt,'.dat'
    open(unit=outJm,file=file_name,access='STREAM',status='UNKNOWN')
    if ( RecNo == 1) then
        REWIND outJm
    end if
    write(outJm) jaco_phys
    close(outJm)

end subroutine WriteQPmodal


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
