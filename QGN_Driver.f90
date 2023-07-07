program main
use QGN_Module
implicit none
real(dp) :: q0(nx,ny,nz) ! Initial condition
complex(dp) :: q_hat(nx/2+1,ny,nz)
real(dp) :: t = 0._dp, dt = 750._dp
integer :: i, j, N0 = 1, Nt = 10, diag_freq = 10000, chkpt_freq = 1, out_freq = 1 
character(len=32) :: arg
character(len=16) :: file_name

do i = 1, command_argument_count()
    call get_command_argument(i,arg)
    if (len_trim(arg) == 0) exit 
    read (arg,*) N0 
end do
print *,'---------------------------------------------------------------'
print *,' Starting iteration  '
print *, N0
print *,'---------------------------------------------------------------'

if( N0 > 1 ) then
   ! Initialize from a saved file
   write (file_name,'(A2,I0.9,A4)') "q.",N0,'.dat'
   print *,'---------------------------------------------------------------'
   print *,' Opening q0 file  '
   print *,  file_name
   print *,'---------------------------------------------------------------'
   open(unit=10,file=file_name,access='STREAM',status='OLD')
   read(unit=10) q0
   close(unit=10)

else
   ! Define the initial condition below. 
   ! call RANDOM_NUMBER(q0)
   ! do i=1,nz
   !     q0(:,:,i) = 1.E-9_dp*(q0(:,:,i) - sum(q0(:,:,i))/real(nx*ny,dp))
   ! end do
   do i=1,nx
       do j=1,ny
           q0(i,j,:) = 1.E-5_dp*sin(48._dp*pi*real(i-1,dp)/real(nx,dp))
           q0(i,j,:) = q0(i,j,:) + 2.E-5_dp*cos(32._dp*pi*real(j-1,dp)/real(ny,dp))
       end do
   end do

end if

call Initialize(q0,q_hat)
call WriteQP(q_hat,1,0) ! Writes q and psi to q.dat and p.dat.
do i=N0,N0+Nt
    call TimeStep(q_hat,t,dt)
    if( mod(i,out_freq)==0 ) print *, 'Time since inception = ',real(t/86400._dp),'days; time step size = ',real(dt),'s'
    if( mod(i,diag_freq)==0 ) then
        print *, ' Diagnostics, RecNo = ', i/diag_freq
        ! If you want any diagnostics beyond QP and QP Modal, you have to write your own subroutine and call it here.
    end if
    if( mod(i,chkpt_freq)==0 ) then
        call WriteQP(q_hat,1,i) ! Writes q and psi to q.dat and p.dat.
                              ! Good for restarts.
        call WriteQPmodal(q_hat,1,i) ! Writes q and psi to q.dat and p.dat.
    end if
end do
call Cleanup

end program main
