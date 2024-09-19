    ! RK stages for ARK4(3)6L[2]SA 
    complex(dp), dimension(nx/2+1,ny,nz) :: N1, N2, N3, N4, N5, N6
    complex(dp), dimension(nx/2+1,ny,nz) :: L1, L2, L3, L4, L5, L6
    complex(dp), dimension(nx/2+1,ny,nz) :: q_tmp, Mq
    real(dp), dimension(nx/2+1,ny), save :: k8
    real(dp), save :: err0, err1 ! current and previous errors for adaptive time step
    real(dp), dimension(6,6), save :: ae, ai ! explicit & implicit RK coefficients
    real(dp), dimension(6), save :: b, be ! be is error coefficients
    logical, save  :: FirstCall = .TRUE.
    logical :: reject

    reject = .TRUE.

    if (FirstCall) then
        k8 = k2**4
        err0 = TOL
        ae = 0._dp
        ae(2,1) = 0.5_dp
        ae(3,1) = 13861._dp/62500._dp
        ae(3,2) = 6889._dp/62500._dp
        ae(4,1) = -116923316275._dp/2393684061468._dp
        ae(4,2) = -2731218467317._dp/15368042101831._dp
        ae(4,3) = 9408046702089._dp/11113171139209._dp
        ae(5,1) = -451086348788._dp/2902428689909._dp
        ae(5,2) = -2682348792572._dp/7519795681897._dp
        ae(5,3) = 12662868775082._dp/11960479115383._dp
        ae(5,4) = 3355817975965._dp/11060851509271._dp
        ae(6,1) = 647845179188._dp/3216320057751._dp
        ae(6,2) = 73281519250._dp/8382639484533._dp
        ae(6,3) = 552539513391._dp/3454668386233._dp
        ae(6,4) = 3354512671639._dp/8306763924573._dp
        ae(6,5) = 4040._dp/17871._dp
        b(1) = 82889._dp/524892._dp
        b(2) = 0._dp
        b(3) = 15625._dp/83664._dp
        b(4) = 69875._dp/102672._dp
        b(5) = -2260._dp/8211._dp
        b(6) = 0.25_dp
        ai = 0._dp
        ai(2,1) = 0.25_dp
        ai(2,2) = 0.25_dp ! All diagonal elements except ai(1,1) are 1/4
        ai(3,1) = 8611._dp/62500._dp
        ai(3,2) = -1743._dp/31250._dp
        ai(4,1) = 5012029._dp/34652500._dp
        ai(4,2) = -654441._dp/2922500._dp
        ai(4,3) = 174375._dp/388108._dp
        ai(5,1) = 15267082809._dp/155376265600._dp
        ai(5,2) = -71443401._dp/120774400._dp
        ai(5,3) = 730878875._dp/902184768._dp
        ai(5,4) = 2285395._dp/8070912._dp
        ai(6,:) = b
        be(1) = 31666707._dp/9881966720._dp
        be(2) = 0._dp
        be(3) = -256875._dp/105007616._dp
        be(4) = -2768025._dp/128864768._dp
        be(5) = 169839._dp/3864644._dp
        be(6) = -5247._dp/225920._dp
        FirstCall = .FALSE.
    end if

    ! First RK stage, t=0
    call GetRHS(q_hat, N1)
    do k=1,nz
        L1(:,:,k) = -(A2*k2 + A8*k8)*q_hat(:,:,k)
    end do
    
do while (reject)
    do k=1,nz
        Mq(:,:,k) = 1._dp/(1._dp + 0.25_dp*dt*(A2*k2+A8*k8))
    end do
    ! Second RK stage
    q_tmp = Mq*(q_hat + dt*(ae(2,1)*N1+ai(2,1)*L1))
    call GetRHS(q_tmp,N2)
    do k=1,nz
        L2(:,:,k) = -(A2*k2 + A8*k8)*q_tmp(:,:,k)
    end do
    ! Third RK stage
    q_tmp = Mq*(q_hat + dt*(ae(3,1)*N1+ae(3,2)*N2 &
                           +ai(3,1)*L1+ai(3,2)*L2))
    call GetRHS(q_tmp,N3)
    do k=1,nz
        L3(:,:,k) = -(A2*k2 + A8*k8)*q_tmp(:,:,k)
    end do
    ! Fourth RK stage
    q_tmp = Mq*(q_hat + dt*(ae(4,1)*N1+ae(4,2)*N2+ae(4,3)*N3 &
                          &+ai(4,1)*L1+ai(4,2)*L2+ai(4,3)*L3 ))
    call GetRHS(q_tmp,N4)
    do k=1,nz
        L4(:,:,k) = -(A2*k2 + A8*k8)*q_tmp(:,:,k)
    end do
    ! Fifth RK stage
    q_tmp = Mq*(q_hat+dt*(ae(5,1)*N1+ae(5,2)*N2+ae(5,3)*N3+ae(5,4)*N4 &
                         +ai(5,1)*L1+ai(5,2)*L2+ai(5,3)*L3+ai(5,4)*L4 ))
    call GetRHS(q_tmp,N5)
    do k=1,nz
        L5(:,:,k) = -(A2*k2 + A8*k8)*q_tmp(:,:,k)
    end do
    ! Sixth RK stage
    q_tmp = Mq*(q_hat + dt*(ae(6,1)*N1+ae(6,2)*N2+ae(6,3)*N3 &
                           +ae(6,4)*N4+ae(6,5)*N5+ai(6,1)*L1 &
                           +ai(6,2)*L2+ai(6,3)*L3+ai(6,4)*L4+ai(6,5)*L5 ))
    call GetRHS(q_tmp,N6)
    do k=1,nz
        L6(:,:,k) = -(A2*k2 + A8*k8)*q_tmp(:,:,k)
    end do
    ! Error control, upper layer only
    spec = be(1)*(N1+L1)+be(3)*(N3+L3)  &
          +be(4)*(N4+L4)+be(5)*(N5+L5)+be(6)*(N6+L6)
    Status = DftiComputeBackward(s2g_q, specE, gridE)
    err1 = dt*MAXVAL(ABS(grid))/real(nx*ny,dp)
    if( err1>TOL ) then
        dt = 0.5_dp*dt
        reject = .TRUE.
    else
        ! RK, compute update
        q_hat = q_hat + dt*(b(1)*(N1+L1)+b(3)*(N3+L3) &
                           +b(4)*(N4+L4)+b(5)*(N5+L5) &
                           +b(6)*(N6+L6))
        q_hat(1,1,:) = cmplx(0._dp,0._dp)
        t = t + dt
        ! Stepsize adjustment PI.3.4, divide by 4 for 4th order method with 3rd embedded
        dt = MIN(dt*((0.2_dp*TOL/err1)**0.075_dp)*((err0/err1)**0.1_dp),dt_max)
        err0 = err1
        reject = .FALSE.
    end if
end do
