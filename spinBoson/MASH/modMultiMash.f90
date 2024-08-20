module modMultiMash
    use params
    use functions
    use force

    implicit none

    ! real(kind=wp1) :: alphaN, hN
    real(kind=wp1) :: basis(nsys,nsys)

    contains

    !! potential surface, for a given q
    subroutine mash_potential(q, qe, pe, vtot)
        !! eq 9 in ref 1
        real(kind=wp1), intent(in) :: q(nbath_total), qe(nsys), pe(nsys)
        real(kind=wp1), intent(out) :: vtot
        integer :: a, ierror
        real(kind=wp1), allocatable :: U(:,:), Vad(:)

        allocate(U(nsys, nsys), Vad(nsys), stat=ierror)

        call pot(q, U)
        call eigsys(U, Vad)
        call cstate(qe, pe, U, a)
        vtot = Vad(a)

        deallocate(U, stat=ierror)
        deallocate(Vad, stat=ierror)

    end subroutine  mash_potential

   function ham_phase_space(q, p, qe, pe)
        !! Function to compute Hamiltonian at a phase-space point
        real(kind=wp1) :: ham_phase_space
        real(kind=wp1), intent(in) :: q(nbath_total), p(nbath_total), qe(nsys), pe(nsys)
        real(kind=wp1) :: vtot

        call mash_potential(q, qe, pe, vtot)
        ham_phase_space = vtot + half*sum(p*p) 
        !! reduced description, mass not needed
    end function ham_phase_space

    !! determine current adiabatic state
    subroutine cstate(qe, pe, U, a)
        !! Get current adiabatic state with the highest population |c_a|^2
        real(kind=wp1), intent(in) :: qe(nsys), pe(nsys), U(nsys,nsys)
        integer, intent(out) :: a
        real(kind=wp1), allocatable :: qa(:), pa(:)
        complex(kind=wp1), allocatable :: ca(:)
        integer :: ierror

        allocate(qa(nsys), pa(nsys), ca(nsys), stat=ierror)
        !! Convert to adiabatic rep.
        qa = matmul(qe, U)
        pa = matmul(pe, U)
        ca = cmplx(qa, pa, kind=wp1)

        call cstate_ad(ca, a)

        deallocate(qa, stat=ierror)
        deallocate(pa, stat=ierror)
        deallocate(ca, stat=ierror)

    end subroutine

    subroutine cstate_ad(ca, a)
        !! Get current adiabatic state with the highest population |c_a|^2
        ! implicit none
        complex(kind=wp1), intent(in) :: ca(:)
        integer, intent(out) :: a

        integer :: ierror
        real(kind=wp1), allocatable :: pop(:)
        allocate(pop(nsys))
        ! Convert to adiabatic rep.
        pop = abs(ca)**2
        a = maxloc(pop, dim=1)
        deallocate(pop, stat=ierror)
    end subroutine cstate_ad

    subroutine cstate_ad_2(ca, a, b)
        !! Get second-highest population
        complex(kind=wp1), intent(in) :: ca(:)
        integer :: a, b, ierror
        real(kind=wp1), allocatable :: pop(:)

        allocate(pop(nsys), stat=ierror)
        ! Convert to adiabatic rep.
        pop = abs(ca)**2
        pop(a) = rzero
        !! otherwise implementing sorting algorithms can be brutal
        b = maxloc(pop, dim=1)
    
        deallocate(pop, stat=ierror)
    end subroutine cstate_ad_2

    !! compute density matrix elements
    subroutine density_mash(c, density)
        complex(kind=wp1), intent(inout) :: c(:)
        complex(kind=wp1), intent(inout) :: density(:,:)
        real(kind=wp1) :: ninv
        integer :: i,j

        ninv = one/nsys
        density = czero
        !! diagonal terms, Eq 29 in Ref 1
        do i=1,nsys
            density(i,i) = ninv + alphaN*(abs(c(i))**2 - ninv)
        end do
        !! coherences, Eq 34 in Ref 1
        do i=1,nsys
            do j=i+1,nsys
                density(i,j) = alphaN*conjg(c(i))*c(j)
                density(j,i) = conjg(density(i,j))
            end do
        end do

    end subroutine density_mash

    subroutine get_density(q, qe, pe, density, rep)
        !! density in diabatic or adiabatic representation
        real(kind=wp1), intent(in) :: q(:), qe(:), pe(:)
        complex(kind=wp1), intent(out) :: density(:,:)
        character, intent(in) :: rep !! diabatic('d') or adiabatic('a')

        real(kind=wp1), allocatable :: Vad(:), U(:,:)
        complex(kind=wp1), allocatable :: c(:)
        integer :: ierror, i

        real(kind=wp1), allocatable :: identity(:,:)

        allocate(identity(nsys, nsys), stat=ierror)
        allocate(c(nsys), stat=ierror)

        identity = rzero
        forall(i=1:nsys) identity(i,i) = one 

        if (rep.eq.'a') then
            allocate(Vad(nsys), U(nsys,nsys))

            call pot(q, U)
            call eigsys(U, Vad)

            !! c = qe + i*pe
            c = cmplx(matmul(qe,U), matmul(pe,U), kind=wp1)
            call density_mash(c, density)
            
            deallocate(Vad, stat=ierror)
            deallocate(U, stat=ierror)

        else if (rep.eq.'d') then
            c = cmplx(qe, pe, kind=wp1)
            call density_mash(c, density)
        else
            stop 'wrong representation'
        end if
        deallocate(c, stat=ierror)
    end subroutine get_density

    subroutine get_popplation(q, qe, pe, pop, rep)
        real(kind=wp1), intent(in) :: q(:), qe(:), pe(:)
        real(kind=wp1), intent(inout) :: pop(:)
        character, intent(in) :: rep
        integer :: i, ierror

        complex(kind=wp1), allocatable :: density(:,:)
        allocate(density(nsys, nsys), stat=ierror)

        call get_density(q, qe, pe, density, rep)
        do i=1,nsys
            pop(i) = real(density(i,i))
        end do

        deallocate(density, stat=ierror)
    end subroutine get_popplation

    !! verlet propagation 
    subroutine move_p(p, dvdq, dt)
        !! Eq E1b or E1d in Ref 1
        real(kind=wp1), intent(inout) :: p(:)
        real(kind=wp1), intent(in) :: dvdq(:), dt

        p = p - dvdq*dt
    end subroutine move_p

    subroutine move_q(q, p, dt)
        !! Eq E1c in Ref 1
        real(kind=wp1), intent(inout) :: q(:), p(:)
        real(kind=wp1), intent(in) :: dt

        q = q + dt*p/mass !! mass is set to 1 or in reduced units
    end subroutine move_q

    subroutine move_e(qe, pe, Vad, U, dt)
        !! Eq 4 in Ref 1
        real(kind=wp1), intent(inout) :: qe(:), pe(:)
        real(kind=wp1), intent(in) :: Vad(:), U(:,:)
        real(kind=wp1), intent(in) :: dt
        integer :: ierror
   
        complex(kind=wp1), allocatable :: c(:)
        allocate(c(nsys), stat=ierror)

        !! U transforms from adiabatic to diabatic basis
        !! c = q + ip
        qe = matmul(qe, U)
        pe = matmul(pe, U)
        c = cmplx(qe, pe, kind=wp1)
        c = exp(-ima * dt * Vad) * c
        !! back to adiabatic basis
        qe = matmul(U, real(c)) 
        pe = matmul(U, aimag(c))

        deallocate(c, stat=ierror)
    end subroutine move_e

    subroutine verlet(q, p, qe, pe, Vad, U, dvdq, a, dt)
        real(kind=wp1), intent(inout) :: q(:), p(:), qe(:), pe(:), Vad(:), &
                                         U(:,:), dvdq(:)
        integer, intent(in) :: a !! state
        real(kind=wp1), intent(in) :: dt
        real(kind=wp1) :: dt2

        dt2 = dt*half
        !! E1a - E1e in Ref 1
        call move_e(qe, pe, Vad, U, dt2)
        call move_p(p, dvdq, dt2)
        call move_q(q, p, dt)
        call pot(q, U)
        call eigsys(U, Vad)
        call grad_adiabatic_diagonal(q, U, a, dvdq)
        call move_p(p, dvdq, dt2)
        call move_e(qe, pe, Vad, U, dt2)

    end subroutine verlet

    function deltaP(ca, a)
        ! implicit none
        real(kind=wp1) :: deltaP
        complex(kind=wp1) :: ca(nsys)
        integer :: a, b

        if ((a<=0) .or. (a>nsys)) then
            print*, 'something wrong with indexing'
            stop
        endif

        call cstate_ad_2(ca, a, b)
        deltaP = abs(ca(a))**2 - abs(ca(b))**2

    end function deltaP

    subroutine evolve(q, p, qe, pe, Vad, U, dvdq, a, dtbase, ierr3)
        real(kind=wp1), intent(inout) :: q(nbath_total), p(nbath_total), &
                                         qe(nsys), pe(nsys), Vad(nsys), &
                                         U(nsys,nsys), dvdq(nbath_total)
        integer, intent(inout) :: a, ierr3
        real(kind=wp1), intent(in) :: dtbase
        real(kind=wp1), allocatable :: q0(:),p0(:),qe0(:),pe0(:), &
                                       Vad0(:), U0(:,:), dvdq0(:)
        complex(kind=wp1), allocatable :: ca0(:), ca1(:)
        integer :: b, maxhop, ierror, iii, jj
        logical :: accepted
        real(kind=wp1) :: tx, tr, tl, tm, fl , fr, fm, dt_here

        allocate(q0(nbath_total), stat=ierror) 
        allocate(p0(nbath_total), stat=ierror)
        allocate(qe0(nsys), stat=ierror)
        allocate(pe0(nsys), stat=ierror)
        allocate(Vad0(nsys), U0(nsys, nsys), dvdq0(nbath_total), stat=ierror)
        allocate(ca0(nsys), ca1(nsys), stat=ierror)

        ierr3 = 0
        maxhop = 30
        dt_here = dtbase
        !! Appendix E in Ref 1
        do iii=1,maxhop ! Limit number of hops to look for in a timestep dt
            ! Store initial values
            call store_data(q, p, qe, pe, Vad, U, dvdq, q0, p0, qe0, pe0, Vad0, U0, dvdq0)
        
            ! Store initial adiabatic wavefunction
            ca0 = cmplx(matmul(qe,U), matmul(pe,U), kind=wp1)
            ! print*,'test: ', size(matmul(qe,U)), sum(ca0)
            ! Attempt full step
            call verlet(q, p, qe, pe, Vad, U, dvdq, a, dt)

            ! Calculate new active state
            ca1 = cmplx(matmul(qe,U),matmul(pe,U), kind=wp1)
            call cstate_ad(ca1, b)

            if (a.eq.b) then
                ! Stayed on state - we're done
                exit
            else
                ! States have changed - find crossing time with bisection root search
                tl = rzero
                tr = dt_here
                fl = deltaP(ca0, a)
                fr = deltaP(ca1, a)
                do jj=1,10
                    tm = (tl+tr)*half ! Mid point
                    call store_data(q0, p0, qe0, pe0, Vad0, U0, dvdq0, q, p, qe, pe, Vad, U, dvdq)
                    call verlet(q, p, qe, pe, Vad, U, dvdq, a, tm)
                    ca1 = cmplx(matmul(qe,U), matmul(pe,U), kind=wp1)
                    fm = deltaP(ca1, a)
                    if (fm.gt.0) then
                        tl = tm
                    else
                        tr = tm
                    end if
                end do
                ! print*,a,b 
                call cstate_ad_2(ca1, a, b)
                ! print*,a,b
                call hopping(q, p, ca1, a, b, Vad, U, accepted)
                if (accepted) then
                    tx = tr
                else
                    tx = tl
                end if
                call store_data(q0, p0, qe0, pe0, Vad0, U0, dvdq0, q, p, qe, pe, Vad, U, dvdq)
                call verlet(q, p, qe, pe, Vad, U, dvdq, a, tx)
                ca1 = cmplx(matmul(qe,U), matmul(pe,U), kind=wp1)
                call hopping(q, p, ca1, a, b, Vad, U, accepted)
                if (accepted) then
                    call grad_adiabatic_diagonal(q, U, b, dvdq)
                    a = b
                end if
                dt_here = dt_here - tx
            end if
            if (iii.eq.maxhop) then
                ierr3=1
            end if
        end do

        deallocate(ca0, stat=ierror)
        deallocate(ca1, stat=ierror)
        deallocate(q0, stat=ierror)
        deallocate(p0, qe0, pe0, Vad0, U0, stat=ierror)

    end subroutine evolve

    !! check hopping
    subroutine hopping(q, p, cad, n, m, Vad, U, accepted)
        ! implicit none
        real(kind=wp1), intent(inout) :: p(:)
        complex(kind=wp1), intent(inout) :: cad(:)
        real(kind=wp1), intent(in) :: q(:), Vad(:), U(:,:)
        logical, intent(out) :: accepted
        integer, intent(in) :: n, m
        integer :: ierror

        real(kind=wp1), allocatable :: d(:,:,:), dj(:), pnac(:), porth(:)
        real(kind=wp1) :: Ekin, Vdiff

        allocate(d(nbath_total, nsys, nsys),dj(nbath_total),pnac(nbath_total),porth(nbath_total), stat=ierror)

        call nonAdCoupDirection(q, cad, Vad, U, n, m, dj)
        dj = dj/sqrt(mass)
        
        ! Use mass-scaled momenta
        p = p/sqrt(mass)
        if (nbath_total.eq.1) then
            pnac = p
        else
            pnac = mydot(p,dj)/mydot(dj,dj)*dj
        end if
        porth = p - pnac
        Ekin = half*sum(pnac**2) 
        Vdiff = Vad(m)-Vad(n)
        if ((Ekin-Vdiff).gt.rzero) then
            ! Rescale momentum 
            pnac = sqrt(two*(Ekin-Vdiff)) & ! (no masses since p is mass-scaled)
                * pnac/mynorm(pnac)
            p = porth + pnac
            accepted = .true.
        else
            ! Reverse momentum along NAC vector
            pnac = -pnac
            p = porth + pnac
            accepted = .false.
        end if
        p = p*sqrt(mass)

        deallocate(d, dj, pnac, porth, stat=ierror)

    end subroutine hopping

    !! save data
    subroutine store_data(q, p, qe, pe, Vad, U, dvdq, q0, p0, qe0, pe0, Vad0, U0, dvdq0)
        ! implicit none
        real(kind=wp1), intent(inout) :: q(:), p(:), qe(:), pe(:), Vad(:), U(:,:), &
                                         dvdq(:), q0(:), p0(:), qe0(:), pe0(:), & 
                                         Vad0(:), U0(:,:), dvdq0(:)   

        q0 = q
        p0 = p
        qe0 = qe
        pe0 = pe
        Vad0 = Vad
        U0 = U
        dvdq0 = dvdq
    end subroutine store_data

    subroutine store_temp(q, p, qe, pe, it, qt, pt, qet, pet)
        ! implicit none
        real(kind=wp1), intent(in) :: q(:), p(:), qe(:), pe(:)
        integer, intent(in) :: it
        real(kind=wp1), intent(inout) :: qt(:,:), pt(:,:), qet(:,:), pet(:,:)

        qt(it,:) = q
        pt(it,:) = p
        qet(it,:) = qe
        pet(it,:) = pe

    end subroutine store_temp

end module modMultiMash