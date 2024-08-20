module params
    !!use functions
    implicit none
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!To be changed as a variable!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: wp1 = selected_real_kind(15,37)
    integer, parameter :: wp2 = selected_real_kind(16)
    integer, parameter :: k15 = selected_int_kind(15)


    integer(kind=k15), parameter :: lmax = 1E6
    integer , parameter :: bath = 100
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real(kind=wp1), parameter :: kelvintocminv = 0.695028_wp1
    real(kind=wp1), parameter :: cm_inv_to_fs = 5305.16477_wp1

    !!!!double precision numbers
    real(kind=wp1), parameter :: rzero = 0.0_wp1
    complex(kind=wp1), parameter :: czero = complex(rzero, rzero)
    real(kind=wp1), parameter :: quart = 0.25_wp1
    real(kind=wp1), parameter :: half = 0.5_wp1
    real(kind=wp1), parameter :: tquart = 0.75_wp1
    real(kind=wp1), parameter :: one = 1.0_wp1
    real(kind=wp1), parameter :: oah = 1.5_wp1
    real(kind=wp1), parameter :: two = 2.0_wp1
    real(kind=wp1), parameter :: three = 3.0_wp1
    real(kind=wp1), parameter :: four = 4.0_wp1
    real(kind=wp1), parameter :: five = 5.0_wp1
    real(kind=wp1), parameter :: six = 6.0_wp1
    real(kind=wp1), parameter :: eight = 8.0_wp1
    real(kind=wp1), parameter :: pi = four*atan(one)
    real(kind=wp1), parameter :: pi2 = pi + pi

    !!!!!time step and total time and their derivatives
    integer, parameter :: nsys = 2, model = 1
    real(kind=wp1), parameter :: gamma = half * (sqrt(three) - one)
    real(kind=wp1), parameter :: dt = 0.0005_wp1, dt2 = dt*dt, &
                                 dt205 = dt2*half, dthalf = dt * half
    integer, parameter :: tlen = nint(20.0_wp1/dt)
    real(kind=wp1), parameter :: ilmax = one/real(lmax,kind=wp1)

    !!!!!imaginary numbers and time derivatives
    complex(kind=wp1), parameter :: ima = complex(rzero, one)
    complex(kind=wp1), parameter :: cone = complex(one, rzero)
    complex(kind=wp1), parameter :: imdt = ima*dt, im2 = ima*pi2

    !!!!!bath and system parameters
    real(kind=wp1), parameter :: omegac = 1.0_wp1
    real(kind=wp1), parameter ::  lambda = 0.125_wp1
    real(kind=wp1), parameter :: capomega = atan(20.0_wp1)
    real(kind=wp1), parameter :: beta = 0.25 !!
    real(kind=wp1), parameter :: e11 = 1.0_wp1 * omegac, j12 = 1.0_wp1 * omegac

    integer :: int11, int2
    real(kind=wp1), parameter :: s1=sqrt(four*lambda*capomega/(pi*real(bath,kind=wp1)))
    real(kind=wp1), parameter :: temp1=capomega/real(bath,kind=wp1)
    real(kind=wp1), dimension(bath), parameter :: wi=(/( tan(real(int11,kind=wp1)*temp1) ,int11=1,bath )/)
    real(kind=wp1), dimension(bath), parameter :: ci=s1*wi
    real(kind=wp1), dimension(bath), parameter :: wi2=wi*wi

    real(kind=wp1), parameter :: filter = three
    real(kind=wp1), parameter :: eps = real(1e-9,kind=wp1)
    real(kind=wp1), parameter :: epsl = real(10000.0,kind=wp1)

    real(kind=wp1), dimension(2,2), parameter :: sigmazr = reshape( [ one, rzero, &
                                                                     rzero, -one], shape(sigmazr), order=(/2,1/))
    complex(kind=wp1), dimension(2,2), parameter :: sigmazc = sigmazr + czero

end module params
