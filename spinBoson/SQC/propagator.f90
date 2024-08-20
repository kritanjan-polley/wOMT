module propagator
    use params
    use functions
    use omp_lib
    implicit none
    ! external :: cgemv

    integer:: i, j, k, l, ii
    integer :: i1, k1, rint, int1
    real(kind=wp1), dimension(nsys) :: local_n
    real(kind=wp1), dimension(bath) :: dp, dpp
    complex(kind=wp1), dimension(nsys) :: md, mdd
    real(kind=wp1):: cmn1, q11_m_q22_dot, cp, cq
    !$omp threadprivate(int1, i1, k1, rint, local_n, dp, dpp)
    !$omp threadprivate(md, mdd, cmn1, q11_m_q22_dot, cp, cq)

    public 
    contains

    subroutine propagation1d(pnew, qnew, mnew, time, hsys, n1val)
    implicit none
    integer,intent(in) ::time
    real(kind=wp1), intent(in) :: n1val
    complex(kind=wp1), dimension(nsys, time+2) :: mnew
    real(kind=wp1), dimension(bath, time+2)  :: pnew, qnew
    real(kind=wp1), dimension(nsys, nsys), intent(in) :: hsys

    cmn1 = one / (n1val + gamma)

    !Go one step backward
    !Common things first

    local_n(:) = cmn1*abs2(mnew(:,2)) - gamma
    q11_m_q22_dot = four * hsys(1,2) * aimag(mnew(2,2)*conjg(mnew(1,2))) * cmn1
    cq = dot_product(ci, qnew(:,2))
    cp = dot_product(ci, pnew(:,2))

    !! compute derivatives
    md = -ima*matmul(hsys + czero, mnew(:,2))
    md = md - ima*cq*matmul(sigmazc, mnew(:,2))

    mdd = -ima*matmul(hsys + czero, md) 
    mdd = mdd - ima*cq*matmul(sigmazc, md)
    mdd = mdd - ima*cp*matmul(sigmazc, mnew(:,2))


    dp(:) = -wi2(:)*qnew(:,2) - ci(:) * (local_n(1) - local_n(2))
    dpp(:) = -wi2(:)*pnew(:,2) - ci(:) * q11_m_q22_dot

    !Now move backward
    do i1=1,nsys
        mnew(i1,1) = mnew(i1,2) - dt*md(i1) + dt205*mdd(i1)
    end do
    qnew(:,1) = qnew(:,2) - dt*pnew(:,2) + dt205*dp(:)
    pnew(:,1) = pnew(:,2) - dt*dp(:) + dt205*dpp(:)

    !End of backward movement
    !Time for forward movement
    do int1=2,time+1
        !! Common things first, again
        do i1=1,nsys
            local_n(i1) = cmn1*abs2(mnew(i1,i)) - gamma
        end do
        local_n(:) = cmn1*abs2(mnew(:,int1)) - gamma
        q11_m_q22_dot = four * hsys(1,2) * aimag(mnew(2,int1)*conjg(mnew(1,int1))) * cmn1
        cq = dot_product(ci, qnew(:,int1))
        cp = dot_product(ci, pnew(:,int1))

        !! compute derivatives
        md = -ima*matmul(hsys + czero, mnew(:,int1))
        md = md - ima*cq*matmul(sigmazc, mnew(:,int1))
    
        mdd = -ima*matmul(hsys + czero, md) 
        mdd = mdd - ima*cq*matmul(sigmazc, md)
        mdd = mdd - ima*cp*matmul(sigmazc, mnew(:,int1))

        dp(:) = -wi2(:)*qnew(:,int1) - ci(:) * (local_n(1) - local_n(2))
        dpp(:) = -wi2(:)*pnew(:,int1) - ci(:) * q11_m_q22_dot

        ! Move  forward m's, p's and q's
        do i1=1,nsys
            mnew(i1,int1+1) = two*mnew(i1,int1) - mnew(i1,int1-1) + dt2*mdd(i1)
        end do
        qnew(:,int1+1) = two*qnew(:,int1) - qnew(:,int1-1) + dt2*dp(:)
        pnew(:,int1+1) = two*pnew(:,int1) - pnew(:,int1-1) + dt2*dpp(:)

    end do
    end subroutine propagation1d

    function calc_state(n) result(state)
        implicit none
        !! n is the population element
        real(kind=wp1), dimension(nsys), intent(in) :: n
        integer :: state, i
        logical, dimension(nsys) :: whichState

        do i=1,nsys
            whichState(i) = all(abs(n-diag_arr(i)) <= gamma, 1)
        end do

        do i=1,nsys
            if (whichState(i)) then
                state = i
            end if
        end do

        if (all(whichState .eqv. .false.)) then
            state = 0
        end if
    end function calc_state

    function calc_coh_state(n) result(state)
        implicit none
        !! n is the population element
        real(kind=wp1), dimension(nsys), intent(in) :: n
        integer :: mi, mj, mk
        integer, dimension(2) :: state

        state = 0
        do mi=1,nsys
            if ( n(mi) < half + gamma .and. n(mi) >= half - gamma) then
            state(1) = mi
            do mj=1,nsys
                if (mi /= mj) then
                    if ( n(mj) > half + gamma .or. n(mj) <= half - gamma ) cycle
                    do mk = 1,nsys
                        if (mk /= mi .and. mk /= mj) then
                            if ( n(mk) > gamma .or. n(mk) <= -gamma ) go to 99
                        end if
                    end do
                    state(2) = mj
                    99 continue
                end if
            end do
            end if
        end do

        if (state(1) == 0 .or. state(2) == 0) then
            state = 0
        end if

    end function calc_coh_state

end module propagator
