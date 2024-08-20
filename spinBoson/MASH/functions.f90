module functions
    use params
    implicit none
    real(kind=wp1), external :: ddot, dnrm2

    interface check_hermitian
        module procedure check_hermitian_real
        module procedure check_hermitian_cmplx
    end interface

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

    subroutine print2d(array)
        implicit none
        real(kind=wp1), dimension(:,:), intent(in) :: array
        integer :: len1, li
        len1 = size(array,1)
        do li=1,len1
            write (*,'(*(f12.6))') array(li,:)
        end do
    end subroutine print2d

    subroutine eigsys(amat, w)
        !! evaluates eigen values and eigenvectors (first entry in dsyevd is 'V')
        !! input, real symmetric matrix :: amat
        !! w -> eigenvalues, ascending order 
        !! amat (output) -> orthonormal eigenvectors of amat
        external :: dsyevd

        real(kind=wp1), intent(inout) :: amat(:,:)
        real(kind=wp1) :: w(:)
        real(kind=wp1), allocatable :: work(:)
        integer :: lda, n, liwork, lwork, eiginfo, eigint
        integer, allocatable :: iwork(:)

        lda = size(amat, 1) !! leading dimension of amat
        n = size(amat, 2) !! order of the matrix (amat)
        lwork = 1 + 6*n + 2*n**2
        liwork = 3 + 5*n

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

    subroutine check_hermitian_real(amat, opt_str)
        implicit none
        real(kind=wp1), intent(in) :: amat(:,:)
        character(len=*), intent(in), optional :: opt_str
        character(len=30) :: out_str

        if (present(opt_str)) then
            out_str = opt_str
        else
            out_str = ''
        end if

        if (size(amat,2) .ne. size(amat,1)) then
            print*, 'Input matrix must be a square matrix'
            stop
        end if

        if ( abs(sum(amat-transpose(amat))) > eps ) then
            print*, trim(out_str)//trim(' Not a hermitian matrix')
        else
            print*, trim(out_str)//trim(' is a Hermitian matrix')
        end if
    end subroutine check_hermitian_real

    subroutine check_hermitian_cmplx(amat, opt_str)
        implicit none
        complex(kind=wp1), intent(in) :: amat(:,:)
        character(len=*), intent(in), optional :: opt_str
        character(len=30) :: out_str

        if (present(opt_str)) then
            out_str = opt_str
        else
            out_str = ''
        end if

        if (size(amat,2) .ne. size(amat,1)) then
            print*, 'Input matrix must be a square matrix'
            stop
        end if

        if ( abs(sum(amat - conjg(transpose(amat)))) > eps ) then
            print*, trim(out_str)//trim(' Not a hermitian matrix')
        else
            print*, trim(out_str)//trim(' Is Hermitian matrix')
        end if
    end subroutine check_hermitian_cmplx

end module functions