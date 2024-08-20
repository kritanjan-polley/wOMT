program main
    use omp_lib
    use params
    use functions
    use propagator
    implicit none

    !Variable directory
    real(kind=wp1), dimension(:,:,:), allocatable :: pin, qin
    real(kind=wp1), dimension(:,:), allocatable :: randtab, hamil, pop, popf
    complex(kind=wp1), dimension(:,:), allocatable :: randtab2, mmat
    real(kind=wp1), dimension(:), allocatable :: shorttab1, shorttab2, nin_loop, ninitial
    complex(kind=wp1), dimension(:,:,:), allocatable :: coh, cohf
    real(kind=wp1), dimension(:,:,:), allocatable :: qbath, pbath
    !time and save txt file
    real(kind=wp1) :: start, finish, n1_in_loop
    integer ::  ierror, m, fixint
    character(len=100) :: nStr, lamStr, tempStr, tauStr, filename, file_end_st

    !Time counting starts right from here
    call get_walltime(start)

    allocate(randtab(nsys+1, lmax), stat=ierror)
    allocate(randtab2(nsys+1, lmax), stat=ierror)
    allocate(hamil(nsys, nsys))
    allocate(mmat(nsys, tlen+2))

    allocate(qbath(nsys, bath, tlen+2))
    allocate(pbath(nsys, bath, tlen+2))

    allocate(pop(nsys, tlen+1), coh(nsys, nsys, tlen+1))
    allocate(popf(nsys, tlen+1), cohf(nsys, nsys, tlen+1))
    allocate(shorttab1(bath), shorttab2(bath))
    allocate(pin(bath,lmax,nsys), qin(bath,lmax,nsys))
    allocate(nin_loop(nsys+1), ninitial(nsys+1))

    !Initial angles
    call random_number(randtab)
    shorttab1=sqrt(wi/(two*tanh(half*beta*wi))) !! momenta
    shorttab2=sqrt(one/(two*tanh(half*beta*wi)*wi)) !! coordinates

    do j=1,lmax
        do i=1,bath
            do k=1,nsys
                pin(i,j,k) = rand_normal(shorttab1(i))
                qin(i,j,k) = rand_normal(shorttab2(i))
            end do
        end do
    end do

    deallocate(shorttab1, shorttab2)
    randtab2 = im2 * randtab
    deallocate(randtab)

    filename = 'hamil.txt'
    if (access(filename, 'r') == 0) then
        print*, 'Reading Hamiltonian from txt file'
        open(11, file=filename)
        do i = 1,nsys
            read(11, * ) hamil(i,:)
        end do
        read(11,*) file_end_st
        close(11)
        if (file_end_st .ne. 'x') then
            print*, 'Something wrong with reading file'
            stop
        end if
    else
        print*, 'Using user-defined (2 by 2) Hamiltonian'
        hamil = rzero
        hamil(1,1) = 50.0_wp1
        hamil(2,2) = -50.0_wp1
        hamil(1,2) = 100.0_wp1
        hamil(2,1) = 100.0_wp1
    end if
    hamil = hamil/omegac
    call check_hermitian(hamil, "Hamiltonian")
    call print2d(hamil*omegac)

    popf = rzero; cohf = czero; m = 0

    print*, "It's about to start"
    ! call omp_set_num_threads(6)
    print*,"With number of available processors:", omp_get_num_procs()
    print*,"No of bath modes=", bath
    print*, 'zero point energy parameter ', gamma
    print*,"Time step dt=",dt,"(In reduced units)"
    print*, "No. of initial conditions:",lmax
    print*,"Total time of propagation:", nint(tlen*dt), "(In reduced units)"
    print*,"No of time slices in a single propagation:",tlen
    print*,"Bath parameters: lambda=", lambda, '(reduced)'
    print*, "omegac=", omegac,"(cm inverse)"
    if (bath .eq. 1) then
        print*,"The civalue is", ci
        print*,"The wivalue is",  wi
    end if

    ninitial = rzero
    ninitial(2) = one

    !$omp parallel do  private(j,i,k) &
    !$omp private(mmat, fixint) &
    !$omp private(qbath, pbath)  &
    !$omp private(pop, coh, n1_in_loop, nin_loop) &
    !$omp shared(pin, qin, randtab2, ninitial) &
    !$omp reduction(+:m) &
    !$omp reduction(+:popf) &
    !$omp reduction(+:cohf)
    outer : do j=1,lmax
        pop = rzero; coh = czero

        do k=1,nsys
            qbath(k,:,2) = qin(:,j,k)
            pbath(k,:,2) = pin(:,j,k)
        end do

        n1_in_loop = ninitial(1)
        do k=1,nsys
            mmat(k,2) = sqrt((n1_in_loop+gamma)*(ninitial(k+1)+gamma)) * &
                        exp(randtab2(1,j)-randtab2(k+1,j))
        end do

        call propagation1d(pbath, qbath, mmat, tlen, hamil, n1_in_loop)

        do k = 1,nsys
            pop(k,:) = real(mmat(k,2:)*conjg(mmat(k,2:)))/(n1_in_loop+gamma) - gamma
        end do
        do i=1,nsys
            do k=i+1,nsys
                coh(i,k,:) = mmat(i,2:)*conjg(mmat(k,2:))/(n1_in_loop+gamma)
                coh(k,i,:) = conjg(coh(i,k,:))
            end do
        end do

        if ( maxval(sum(pop, dim=1)) > one + 0.01_wp1)  then
            print*, 'bad trajectory', j, ' total population is ', maxval(sum(pop, dim=1))
            call flush()
            m = m + 1
            cycle outer
        end if

        if ( any(abs(pop)>filter) ) then
            print*, 'bad trajectory', j
            call flush()
            m = m + 1
            cycle outer
        end if

        popf = popf + pop
        cohf = cohf + coh

        if (mod(j,(lmax/20)).eq.0) then
            call get_walltime(finish)
            print*,"done lmax=",j, 'after', finish-start, 'seconds'
            call flush()
        end if

    end do outer !end of the lmax loop
    !$omp end parallel do

    popf = popf/real(lmax-m, kind=wp1)
    cohf = cohf/real(lmax-m, kind=wp1)

    print*, real(m,kind=wp1)*ilmax*100.0_wp1,"% bad initial conditions"

    deallocate(randtab2, nin_loop, ninitial)
    deallocate(mmat)
    deallocate(qbath, pbath)
    deallocate(pin, qin, pop, coh)

    nStr = str(nsys)
    lamStr = str(nint(lambda*omegac))
    tempStr = str(nint(one/beta*omegac/kelvintocminv))
    tauStr = str(nint(cm_inv_to_fs/omegac))

    filename = "OMT_"//trim(nStr)//"LS_lambda_"//trim(lamStr)//"_T_"// & 
                trim(tempStr)//"_tau_"//trim(tauStr)//".txt"
    !! filename = 'sqc_population.txt'
    print*, "Writing data file ..."
    open(16, file=filename, form='formatted', status='replace', iostat=ierror)
    do i=1,tlen+1
        write(16,'(f15.6$)')(i-1)*dt
        do k=1,nsys
            write(16,'(f15.6$)') popf(k,i)
        end do
        write(16,'(f15.6$)')cohf(1,2,i)
        if (nsys > 2) then
            write(16,'(f15.6$)')cohf(2,3,i)
        end if
        write(16,*)
    end do
    close(16)

    deallocate(popf, cohf, hamil)

    call get_walltime(finish)
    write(*,*) "wall time     : ",(finish-start),"seconds"

end program main
