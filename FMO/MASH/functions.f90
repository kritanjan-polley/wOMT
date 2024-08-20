module functions
    use params
    implicit none
    real(kind=wp1), external :: ddot, dnrm2

    contains
    subroutine get_walltime(wctime)
        use iso_fortran_env
        implicit none
        real(kind=wp1), intent(inout) :: wctime
        integer :: r, c
        call system_clock(c, r)
        wctime = real(c,kind=wp1)/r
    end subroutine get_walltime

    function str(k)
        implicit none
        integer, intent(in) :: k
        character(len=20) :: str
        write (str, *) k
        str = adjustl(str)
    end function str

    subroutine eigsys(amat, w)
        !! evaluates eigen values and eigenvectors (first entry in dsyevd is 'V')
        !! input, real symmetric matrix :: amat
        !! w -> eigenvalues, ascending order 
        !! amat (output) -> orthonormal eigenvectors of amat
        external :: dsyevd

        real(kind=wp1), intent(inout) :: amat(:,:)
        real(kind=wp1) :: w(:), delta
        real(kind=wp1), allocatable :: work(:), amat_local(:,:)
        integer :: lda, n, liwork, lwork, eiginfo, eigint
        integer, allocatable :: iwork(:)

        lda = size(amat, 1) !! leading dimension of amat
        n = size(amat, 2) !! order of the matrix (amat)
        lwork = 1 + 6*n + 2*n**2
        liwork = 3 + 5*n

        allocate(amat_local(size(amat,1), size(amat,2)))

        if (sum(amat - transpose(amat)) > eps) then
            print*,'Un-symmetric matrix'
            stop
        end if

        if (size(w) /= n) then
            print*,'size mismatch'
            stop
        end if


        allocate(iwork(liwork) , work(lwork), stat=eigint)

        call dsyevd('V', 'U', n, amat, lda, w, work, lwork, iwork, liwork, eiginfo)

        if (eiginfo /= 0) then
            print *, 'Error with eigen decomposition:', eiginfo
            stop
        end if

        deallocate(iwork, work, stat=eigint)


        deallocate(amat_local)

    end subroutine eigsys

    function mydot(x, y)
        !! dot product between two vectors
        implicit none
        ! integer, intent(in) :: lenx, incx, incy
        real(kind=wp1), intent(in) :: x(:), y(:)
        real(kind=wp1) :: mydot

        mydot = ddot(size(x), x, 1, y, 1)
    end function mydot

    function mynorm(x)
        implicit none
        real(kind=wp1), intent(in) :: x(:)
        real(kind=wp1) :: mynorm

        mynorm = dnrm2(size(x), x, 1)
    end function mynorm

    function rand_normal(std)
        !! box-muller
        implicit none
        real(kind=wp1), intent(in) :: std
        real(kind=wp1):: rand_normal,r,theta, temp(2)

        call random_number(temp)
        r = sqrt(-two*log(temp(1)))
        theta = pi2*temp(2)
        rand_normal = std*r*cos(theta)
    end function rand_normal

end module functions