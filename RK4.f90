! RK stages for Bogacki-Shampine
complex(dp), dimension(nx/2+1,ny,nz) :: RHS0, RHS1, RHS2, RHS3
complex(dp)  :: q_tmp(nx/2+1,ny,nz)

! First RK stage, t=0
call GetRHS(q_hat,RHS0)
! Second RK stage t=dt/2
q_tmp = q_hat + (0.5_dp*dt)*RHS0
call GetRHS(q_tmp,RHS1)
! Third RK stage, t=dt/2
q_tmp = q_hat + (0.5_dp*dt)*RHS1
call GetRHS(q_tmp,RHS2)
! Fourth RK stage, t=dt
q_tmp = q_hat + dt*RHS2 
call GetRHS(q_tmp,RHS3)
! RK, compute update
q_hat = q_hat + (dt/6._dp)*(RHS0 + 2._dp*RHS1+2._dp*RHS2+RHS3)
t = t + dt 
