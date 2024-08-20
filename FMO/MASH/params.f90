module params
    implicit none

    !! precision
    integer,parameter::wp1=selected_real_kind(15,100)
    integer,parameter::wp2=1.0D0
    integer,parameter::k15=selected_int_kind(15)

    real(kind=wp1), parameter :: kelvintocminv = 0.695028_wp1
    real(kind=wp1), parameter :: cm_inv_to_fs = 5305.16477_wp1

    !! global parameters
    real(kind=wp1),parameter :: rzero = 0.0_wp1
    real(kind=wp1),parameter :: quart = 0.25_wp1
    real(kind=wp1),parameter :: half = 0.5_wp1
    real(kind=wp1),parameter :: tquart = 0.75_wp1
    real(kind=wp1),parameter :: one = 1.0_wp1
    real(kind=wp1),parameter :: oah = 1.5_wp1
    real(kind=wp1),parameter :: two = 2.0_wp1
    real(kind=wp1),parameter :: three = 3.0_wp1
    real(kind=wp1),parameter :: four = 4.0_wp1
    real(kind=wp1),parameter :: seven = 7.0_wp1
    real(kind=wp1),parameter :: eight = 8.0_wp1
    real(kind=wp1),parameter :: ten = 10.0_wp1
    real(kind=wp1),parameter :: fifteen = 15.0_wp1
    real(kind=wp1),parameter :: onethird = one/three
    real(kind=wp1),parameter :: pi = four*atan(one)
    real(kind=wp1),parameter :: spi = sqrt(pi)
    real(kind=wp1),parameter :: pi2 = pi+pi
    real(kind=wp1), parameter :: obs = one/6.0_wp1

    complex(kind=wp1),parameter :: czero = complex(rzero,rzero)
    complex(kind=wp1),parameter :: ima = complex(rzero,one)

    real(kind=wp1),parameter :: eps=1E-10

    integer, parameter :: initial_state = 1
    integer, parameter :: second_state = 2

    !! input parameters, time
    real(kind=wp1), parameter :: dt = 0.001_wp1
    real(kind=wp1), parameter :: total_time = 10.0_wp1
    integer, parameter :: nsteps = nint(total_time/dt)
    integer, parameter :: num_traj = nint(1e3)

    !! input parameters, system & bath
    integer, parameter :: nsys = 3
    integer, parameter :: nbath = 100
    integer, parameter :: nbath_total = nbath * nsys

    real(kind=wp1), parameter :: omegac = 53.07_wp1, lambda = 50.0_wp1/omegac
    real(kind=wp1), parameter :: capomega = atan(10.0_wp1)
    real(kind=wp1), parameter :: beta = omegac/53.51774_wp1!!208.51037_wp1!!208.51037_wp1!!

    real(kind=wp1), dimension(nsys, nsys), parameter :: Hamil_e(nsys, nsys) = reshape([&
       200.0_wp1, -87.7_wp1, 5.5_wp1, &
       -87.7_wp1, 320.0_wp1, 30.8_wp1, &
       5.5_wp1, 30.8_wp1, 0.0_wp1 ], shape(Hamil_e), order=[2,1])/omegac

    integer :: pi_1, pi_2
    real(kind=wp1), parameter :: pr_1 = sqrt(four*lambda*capomega/(pi*real(nbath,kind=wp1)))
    real(kind=wp1), parameter :: pr_2 = capomega/real(nbath,kind=wp1)
    real(kind=wp1), dimension(nbath), parameter :: wi_temp = (/( tan(real(pi_1,kind=wp1)*pr_2) ,pi_1=1,nbath )/)
    real(kind=wp1), dimension(nbath), parameter :: ci = pr_1*wi_temp
    
    real(kind=wp1), dimension(nbath_total), parameter :: wi = (/ (wi_temp, pi_1=1,nsys) /)
    real(kind=wp1), dimension(nbath_total), parameter :: mass = (/ (one, pi_1=1,nbath_total) /)
    real(kind=wp1), dimension(nbath), parameter :: mw2 = mass(1:nbath)*wi_temp*wi_temp

    real(kind=wp1), dimension(nbath), parameter :: shorttab1 = sqrt(wi_temp/(two*tanh(half*beta*wi_temp))) !p
    real(kind=wp1), dimension(nbath), parameter :: shorttab2 = sqrt(one/(two*tanh(half*beta*wi_temp)*wi_temp)) !q

    !!!
    real(kind=wp1), dimension(nsys), parameter :: nsys_array = (/ (pi_1, pi_1=1,nsys) /)
    real(kind=wp1), parameter :: hN = sum(one/nsys_array)
    real(kind=wp1), parameter :: alphaN = (nsys - one)/(hN - one)

    character(len=30), parameter :: method = 'AA'

end module params