program main
    use params
    use functions
    use force
    use modMultiMash
    implicit none

    !! ref 1 : https://arxiv.org/pdf/2305.08835.pdf

    integer :: mi, mj, err_main, im, jm, track_bad_trj, int1, int2
    real(kind=wp1) :: start, finish
    real(kind=wp1), allocatable :: q(:), p(:), qe(:), pe(:), Et(:), pN(:), population(:)
    real(kind=wp1), allocatable :: qt(:,:), pt(:,:), qet(:,:), pet(:,:), store_pop(:,:)
    real(kind=wp1), allocatable :: pop_step(:,:)
    complex(kind=wp1), allocatable :: coherence(:,:), coherence_step(:,:,:), store_coherence(:,:,:)
    character(len=100) :: nStr, lamStr, tempStr, tauStr, filename, command_arg

    allocate(population(nsys), coherence(nsys,nsys), stat=err_main)
    allocate(store_pop(nsteps + 1, nsys), pop_step(nsteps + 1, nsys), stat=err_main)
    allocate(coherence_step(nsteps + 1, nsys, nsys), store_coherence(nsteps+1, nsys, nsys), stat=err_main)

    call get_walltime(start)

    call get_command_argument(1,command_arg)

    print*, hN, alphaN

    store_pop = rzero
    store_coherence = rzero
    track_bad_trj = 0

    do im = 1, num_traj
        call initial_sampling_el
        call initial_sampling_bath
        call runtraj(q, p, qe, pe, qt, pt, qet, pet, Et, dt, err_main, nsteps, pop_step, coherence_step)

        if (err_main == 0) then
            store_pop = store_pop + pop_step
            store_coherence = store_coherence + coherence_step
        end if
        if (err_main /= 0) then
            track_bad_trj = track_bad_trj + 1
        end if

        if (mod(im,(num_traj/20)).eq.0) then
                print*,"Done with n_traj=",im
        end if

    end do

    store_pop = store_pop/real(num_traj - track_bad_trj, kind=wp1)
    store_coherence = store_coherence/real(num_traj - track_bad_trj, kind=wp1)
    
    print*,track_bad_trj ," out of ",num_traj ," initial conditions are bad"
    print*,(real(track_bad_trj,kind=wp1)/num_traj*100.0_wp1),"% bad initial conditions"

    nStr = str(nsys)
    lamStr = str(nint(lambda*omegac))
    tempStr = str(nint(one/beta*omegac/kelvintocminv))
    tauStr = str(nint(cm_inv_to_fs/omegac))

    filename = "mash_"//trim(nStr)//"LS_lambda_"//trim(lamStr)//"_T_"// & 
                trim(tempStr)//"_tau_"//trim(tauStr)//"_"//trim(method)// &
                trim(command_arg)//".txt"
    open(102,file=filename,status="replace")
    do im = 1, nsteps + 1
        write(102,'(f8.3$)') (im-1)*dt
        do jm = 1,nsys
            write(102,'(f12.7$)') store_pop(im, jm)
        end do
        write(102,'(f12.7$)') real(store_coherence(im, 1, 2))
        write(102,'(f12.7$)') aimag(store_coherence(im, 1, 2))
        if (nsys > 2) then
            write(102,'(f12.7$)') real(store_coherence(im, 2, 3))
            write(102,'(f12.7$)') aimag(store_coherence(im, 2, 3))
        end if
        write(102,*)
    end do
    close(102)

    deallocate(pN, stat = err_main)
    call get_walltime(finish)
    write(*,'(A, f16.6, A)') "wall time  : ",(finish-start)," seconds"

    contains

    function rho_to_c(rho_val)
        implicit none
        real(kind=wp1) :: rho_to_c
        real(kind=wp1), intent(in) :: rho_val

        rho_to_c = (rho_val - one/nsys)/alphaN + one/nsys
    end function rho_to_c

    subroutine initial_sampling_el
        implicit none
        real(kind=wp1), allocatable :: randn(:), qe_temp(:), pe_temp(:), pN(:)
        real(kind=wp1), allocatable :: temp_check(:), rangle(:)

        allocate(pe(nsys), qe(nsys), stat=err_main)

        allocate(randn(2*nsys), pN(nsys), rangle(nsys), stat=err_main)
        allocate(qe_temp(nsys), pe_temp(nsys), temp_check(nsys), stat=err_main)

        if (method =='RM') then
            !! described on page 7 of Ref 1, a bit confusing
            do 
                do int1 = 1, 2*nsys
                    randn(int1) = rand_normal(one)
                end do
                randn = randn/sqrt(sum(randn*randn))
                qe_temp = randn(1:nsys)
                pe_temp = randn(nsys+1:2*nsys)
                temp_check = qe_temp**2 + pe_temp**2
                int1 = 1
                do int2=1,nsys
                    if (temp_check(int2) > temp_check(int1)) then
                        int1 = int2
                    end if
                end do
                if (int1 == initial_state) then
                    pe = pe_temp
                    qe = qe_temp
                    exit
                end if
            end do
        elseif (method =='AA') then
            !! using Eq 29 in Ref 1
            ! pN = (/ (one, int1=1,nsys) /)/nsys*(one-one/alphaN)
            ! pN(initial_state) = one/nsys + one/alphaN*(one-one/nsys)

            do int1=1,nsys
                pN(int1) = rho_to_c(rzero)
            end do
            pN(initial_state) = rho_to_c(half)
            pN(second_state) = rho_to_c(half)

            call random_number(rangle)
            rangle = rangle*pi2

            qe = cos(rangle)*sqrt(pN)
            pe = sin(rangle)*sqrt(pN)
        end if

        ! print*, one/nsys + alphaN*(qe(1)**2 + pe(1)**2 - one/nsys)

        deallocate(qe_temp, pe_temp, temp_check, randn, rangle, pN)
        
    end subroutine initial_sampling_el

    subroutine initial_sampling_bath
        ! implicit none

        allocate(q(nbath_total), p(nbath_total), stat=err_main)
        ! allocate(randsys(nsys), stat=err_main)

        allocate(qt(nsteps + 1, nbath_total), pt(nsteps + 1 , nbath_total), stat=err_main)
        allocate(qet(nsteps + 1, nsys), pet(nsteps + 1, nsys), stat=err_main)
        allocate(Et(nsteps + 1), stat=err_main)

        !! initial bath, Wigner
        do mj=1,nsys
            do mi=1,nbath
                p((mj-1)*nsys + mi) = rand_normal(shorttab1(mi))
                q((mj-1)*nsys + mi) = rand_normal(shorttab2(mi))
            end do
        end do

    end subroutine initial_sampling_bath

    subroutine runtraj(q, p, qe, pe, qt, pt, qet, pet, Et, dt_here, ierr, nsteps, temp_pop, temp_coherence)
        ! implicit none
        integer, intent(in) :: nsteps
        real(kind=wp1), intent(in) :: dt_here
        real(kind=wp1), intent(inout) :: q(nbath_total), p(nbath_total), qe(nsys), pe(nsys)
        real(kind=wp1), intent(out) :: qt(nsteps + 1, nbath_total), pt(nsteps + 1, nbath_total), &
                                       qet(nsteps + 1, nsys), pet(nsteps + 1, nsys), &
                                       Et(nsteps + 1)
        integer, intent(out) :: ierr
        real(kind=wp1), intent(inout) :: temp_pop(nsteps + 1, nsys)
        complex(kind=wp1), intent(inout) :: temp_coherence(nsteps+1, nsys, nsys)

        real(kind=wp1), allocatable :: Vad(:), U(:,:), dvdq(:)
        integer :: a, ierror, kk, jj, ll

        allocate(Vad(nsys), U(nsys, nsys), dvdq(nbath_total), stat=ierror)

        call pot(q, U)
        call eigsys(U, Vad)
        call cstate(qe, pe, U, a)
        call grad_adiabatic_diagonal(q, U, a, dvdq)

        temp_pop = rzero
        do kk = 1, nsteps
            call store_temp(q, p, qe, pe, kk, qt, pt, qet, pet)
            Et(kk) = ham_phase_space(q, p, qe, pe)
            call evolve(q, p, qe, pe, Vad, U, dvdq, a, dt_here, ierr)

            if (ierr > 0) then
                write(*,'(A,I6)') 'Too many hops(bad trajectory)', kk
                exit
            end if

            call get_density(q, qe, pe, coherence, 'd')
            do jj=1,nsys
                temp_pop(kk, jj) = real(coherence(jj,jj))
            end do

            do jj=1,nsys
                do ll=jj+1,nsys
                    temp_coherence(kk,jj,ll) = coherence(jj,ll)
                    temp_coherence(kk,ll,jj) = coherence(ll,jj)
                end do
            end do
        end do

        if (ierr.eq.0) then
            call store_temp(q, p, qe, pe, nsteps + 1, qt, pt, qet, pet)
            Et(nsteps + 1) = ham_phase_space(q, p, qe, pe)

            call get_popplation(q, qe, pe, population, 'd')
            temp_pop(nsteps + 1,:) = population

            call get_density(q, qe, pe, coherence, 'd')
            do jj=1,nsys
                do ll=jj+1,nsys
                    temp_coherence(nsteps+1,jj,ll) = coherence(jj,ll)
                    temp_coherence(nsteps+1,ll,jj) = coherence(ll,jj)
                end do
            end do
        end if

        deallocate(Vad, stat=ierror)
        deallocate(U, stat=ierror)
        deallocate(dvdq, stat=ierror)

    end subroutine runtraj

end program main
