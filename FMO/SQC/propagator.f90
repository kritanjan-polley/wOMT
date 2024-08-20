module propagator
    use params
    use functions
    use omp_lib
    implicit none
    ! external :: cgemv

    integer:: i, j, k, l
    integer :: i1, k1
    real(kind=wp1), dimension(nsys) :: local_n, ndot
    real(kind=wp1), dimension(nsys, bath) :: dp, dpp
    complex(kind=wp1), dimension(nsys) :: md, mdd
    real(kind=wp1):: cmn1
    !$omp threadprivate(i1, k1, local_n, dp, dpp, ndot)
    !$omp threadprivate(md, mdd, cmn1)

    public 
    private :: c_diagonal 
    contains

    subroutine propagation1d(pnew, qnew, mnew, time, hsys, n1val)
    implicit none
    integer,intent(in) ::time
    real(kind=wp1), intent(in) :: n1val
    complex(kind=wp1), dimension(nsys, time+2) :: mnew
    real(kind=wp1), dimension(nsys, bath, time+2)  :: pnew, qnew
    real(kind=wp1), dimension(nsys, nsys), intent(in) :: hsys

    cmn1 = one / (n1val + gamma)

    !Go one step backward
    !Common things first
    do i1=1,nsys
        local_n(i1) = cmn1*abs2(mnew(i1,2)) - gamma
    end do
    ndot = rzero
    do i1=1,nsys
        do k1=1,nsys
            if (i1 .ne. k1) then
                ndot(i1) = ndot(i1) + two*cmn1*hsys(i1,k1)*aimag(conjg(mnew(i1,2))*mnew(k1,2))
            end if
        end do
    end do

    !! compute derivatives
    md = -ima*matmul(hsys + czero, mnew(:,2))
    md = md - ima*matmul(c_diagonal(qnew(:,:,2)) + czero, mnew(:,2))

    mdd = -ima*matmul(hsys + czero, md) 
    mdd = mdd - ima*matmul(c_diagonal(qnew(:,:,2)) + czero, md )
    mdd = mdd - ima*matmul(c_diagonal(pnew(:,:,2)) + czero, mnew(:,2))

    do i1=1,nsys
        dp(i1,:) = -wi2(:)*qnew(i1,:,2) - local_n(i1)*ci(:)
        dpp(i1,:) = -wi2(:)*pnew(i1,:,2) - ndot(i1)*ci(:)
    end do

    !Now move backward
    do i1=1,nsys
        mnew(i1,1) = mnew(i1,2) - dt*md(i1) + dt205*mdd(i1)
        qnew(i1,:,1) = qnew(i1,:,2) - dt*pnew(i1,:,2) + dt205*dp(i1,:)
        pnew(i1,:,1) = pnew(i1,:,2) - dt*dp(i1,:) + dt205*dpp(i1,:)
    end do
    !End of backward movement
    !Time for forward movement
    do i=2,time+1
        !! Common things first, again
        do i1=1,nsys
            local_n(i1) = cmn1*abs2(mnew(i1,i)) - gamma
        end do

        ndot = rzero
        do i1=1,nsys
            do k1=1,nsys
                if (i1 .ne. k1) then
                    ndot(i1) = ndot(i1) + two*cmn1*hsys(i1,k1)*aimag(conjg(mnew(i1,i))*mnew(k1,i))
                end if
            end do
        end do

        !! compute derivatives
        md = -ima*matmul(hsys + czero, mnew(:,i))
        md = md - ima*matmul(c_diagonal(qnew(:,:,i)) + czero, mnew(:,i))
    
        mdd = -ima*matmul(hsys+czero, md) 
        mdd = mdd - ima*matmul(c_diagonal(qnew(:,:,i)) + czero, md)
        mdd = mdd - ima*matmul(c_diagonal(pnew(:,:,i)) + czero, mnew(:,i))

        do i1=1,nsys
            dp(i1,:) = -wi2(:)*qnew(i1,:,i) - local_n(i1)*ci(:)
            dpp(i1,:) = -wi2(:)*pnew(i1,:,i) - ndot(i1)*ci(:)
        end do

        !! Move  forward m's, p's and q's
        do i1=1,nsys
            mnew(i1,i+1) = two*mnew(i1,i) - mnew(i1,i-1) + dt2*mdd(i1)
            qnew(i1,:,i+1) = two*qnew(i1,:,i) - qnew(i1,:,i-1) + dt2*dp(i1,:)
            pnew(i1,:,i+1) = two*pnew(i1,:,i) - pnew(i1,:,i-1) + dt2*dpp(i1,:)
        end do

    end do
    end subroutine propagation1d

    function c_diagonal(pmat) result(output)
        implicit none
        real(kind=wp1), dimension(nsys, bath), intent(in) :: pmat
        real(kind=wp1), dimension(nsys, nsys) :: output
        integer :: ii

        output = rzero
        do ii=1,nsys
            output(ii,ii) = dot_product(ci(:),pmat(ii,:))
        end do
    end function c_diagonal

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
