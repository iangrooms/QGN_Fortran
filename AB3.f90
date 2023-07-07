complex(dp), dimension(nx/2+1,ny,nz), save :: RHS0, RHS1, RHS2
logical, save :: FirstCall = .TRUE.
logical, save :: SecondCall = .FALSE.

if (FirstCall) then
    ! Forward Euler
    call GetRHS(q_hat,RHS0)
    q_hat = q_hat + dt*RHS0
    FirstCall = .FALSE.
    SecondCall = .TRUE.
else if (SecondCall) then
    ! AB2
    call GetRHS(q_hat,RHS1)
    q_hat = q_hat + dt*(1.5_dp*RHS1 - 0.5_dp*RHS0)
    SecondCall = .FALSE.
else
    ! AB3
    call GetRHS(q_hat,RHS2)
    q_hat = q_hat + (dt/12._dp)*(23._dp*RHS2-16._dp*RHS1+5._dp*RHS0)
    RHS0 = RHS1
    RHS1 = RHS2
end if
t = t + dt
