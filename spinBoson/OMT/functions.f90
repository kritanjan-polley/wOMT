module functions
    use params
    implicit none

    interface check_hermitian
        module procedure check_hermitian_real
        module procedure check_hermitian_cmplx
    end interface

    contains 

    subroutine get_walltime(wctime)
        use iso_fortran_env
        implicit none
        integer, parameter :: dp = kind(one)
        real(dp) :: wctime
        integer :: r, c
        call system_clock(c, r)
        wctime = real(c, dp) / real(r, dp)
    end subroutine get_walltime

    function rand_normal(std)
        use omp_lib
        implicit none
        real(kind=wp1), intent(in)::std
        real(kind=wp1):: rand_normal,r,theta,temp(2)
  
        call random_number(temp)
        r = sqrt(-two*log(temp(1)))
        theta = pi2*temp(2)
        rand_normal = std*r*cos(theta)
    end function rand_normal
  
    function rand_uniform(a,b) result(res)
        implicit none
        !! random uniform between a, b
        real(kind=wp1), intent(in) :: a, b
        real(kind=wp1) :: res, temp

        call random_number(temp)
        res = temp * abs(a-b) + min(a,b)
    end function rand_uniform

    function random_int(low, high) result(rint)
        implicit none
        integer, intent(in) :: low, high
        integer :: rint
        real(kind=wp1) :: u

        call random_number(u)
        rint = low + floor(real(high-low + 1, kind=wp1)*u)
    end function random_int

    function str(k)
        implicit none
        integer, intent(in) :: k
        character(len=20) :: str
        write (str, *) k
        str = adjustl(str)
    end function str

    elemental real(kind=wp1) function abs2(x)
        implicit none
        complex(kind=wp1), intent(in) :: x
        abs2 = real(x)*real(x) + aimag(x)*aimag(x)
    end function abs2

    subroutine print2d(array)
        implicit none
        real(kind=wp1), dimension(:,:), intent(in) :: array
        integer :: len1, li
        len1 = size(array,1)
        do li=1,len1
            write (*,'(*(f12.6))') array(li,:)
        end do
    end subroutine print2d

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

    function bar(j) result(k)
        implicit none
        integer, intent(in) :: j
        integer :: k
        
        select case (j)
        case (1)
            k = 2
        case (2)
            k = 1
        case default
            print*, 'invalid input in bar function'
            stop
        end select

    end function bar

end module functions