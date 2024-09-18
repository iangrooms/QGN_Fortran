program main
use QGN_Module
implicit none
complex(dp) :: q_hat(nx/2+1,ny,nz)
real(dp) :: t = 0._dp, dt = 86400._dp, diag_freq = 5._dp*86400._dp, chkpt_freq = 10._dp*86400._dp
integer :: i, j, N0 = 1, Nt = 14400, out_freq = 1
character(len=32) :: arg
character(len=16) :: file_name

real(dp) :: TE
integer :: n_diag, n_chkpt
integer :: last_diag = 0
integer :: last_chkpt = 10

do i = 1, command_argument_count()
    call get_command_argument(i,arg)
    if (len_trim(arg) == 0) exit
    read (arg,*) N0
end do
print *,'---------------------------------------------------------------'
print *,' Starting iteration  '
print *, N0

call Initialize(q_hat,N0)

if( N0 == 1 ) then
   call WriteQP(q_hat,0) ! Writes q and psi
   call WriteKE(q_hat,1) ! Writes KE.dat
end if

print *,'---------------------------------------------------------------'

do i=N0,N0+Nt-1
    call TimeStep(q_hat,t,dt)

    if( mod(i,out_freq)==0 ) then
        call GetTE(q_hat, TE)
        print *, i-N0+1,',',real(t/86400._dp),',',real(dt),',',real(TE)
    end if

    n_diag = INT(t/diag_freq)
    if( n_diag > last_diag ) then
        !print *, 'Diagnostics at t ', real(t/86400._dp), 'days'
        call WriteKE(q_hat, 0) ! Writes KE.dat
        call WriteWB(q_hat, i)
        last_diag = n_diag
    end if
    n_chkpt = INT(t/chkpt_freq)
    if( n_chkpt > last_chkpt ) then
        call WriteQP(q_hat,i) ! Writes q and psi
        !call WriteQPmodal(q_hat,i)
        !call WriteWB(q_hat, i)
        last_chkpt = n_chkpt
    end if

    dt = MAX(1._dp, MIN(dt, (n_diag+1)*diag_freq - t, (n_chkpt+1)*chkpt_freq - t))
end do

n_chkpt = INT(t/chkpt_freq)
if( n_chkpt > last_chkpt ) then
    call WriteQP(q_hat,i) ! Writes q and psi
    !call WriteQPmodal(q_hat,i)
    call WriteWB(q_hat, i)
end if

call Cleanup

end program main
