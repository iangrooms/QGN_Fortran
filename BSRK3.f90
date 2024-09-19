    ! RK stages for Bogacki-Shampine
    complex(dp), dimension(nx/2+1,ny,nz), save :: RHS0
    complex(dp), dimension(nx/2+1,ny,nz) :: RHS1, RHS2, RHS3, err_spec
    complex(dp)  :: q_tmp(nx/2+1,ny,nz)
    real(dp) :: err_phys(nx,ny,nz)
    real(dp), save :: r0, r1 ! current and previous errors for adaptive time step
    logical, save  :: FirstCall = .TRUE.
    logical :: reject

    reject = .TRUE.

    if (FirstCall) then
        ! First RK stage, t=0
        call GetRHS(q_hat,RHS0)
        r0 = TOL
        FirstCall = .FALSE.
    end if

do while (reject)
    ! Second RK stage t=dt/2
    q_tmp = q_hat + 0.5_dp*dt*RHS0
    call GetRHS(q_tmp,RHS1)
    ! Third RK stage, t=3dt/4
    q_tmp = q_hat + 0.75_dp*dt*RHS1
    call GetRHS(q_tmp,RHS2)
    ! Fourth RK stage, FSAL, t=dt
    q_tmp = q_hat + (2._dp*dt/9._dp)*RHS0 + (1._dp*dt/3._dp)*RHS1 + (4._dp*dt/9._dp)*RHS2
    call GetRHS(q_tmp,RHS3)
    ! Error control
    err_spec = ( (5._dp/72._dp)*RHS0 - (1._dp/12._dp)*RHS1 &
               - (1._dp/9._dp)*RHS2 + (1._dp/8._dp)*RHS3 )
    spec = err_spec
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    err_phys = grid/real(nx*ny,dp)
    r1 = dt*MAXVAL(ABS(err_phys))
    if( r1>TOL ) then
        dt = 0.5_dp*dt
        reject = .TRUE.
    else
        ! RK, compute update
        q_hat = q_tmp
        t = t + dt
        ! FSAL
        RHS0 = RHS3
        ! Stepsize adjustment
        dt = MIN(dt*((0.2_dp*TOL/r1)**0.1_dp)*((r0/r1)**(0.4_dp/3._dp)),dt_max)
        r0 = r1
        reject = .FALSE.
    end if
end do

