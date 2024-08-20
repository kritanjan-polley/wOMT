module mod_spectra_heom4
!! Chen, Zheng, Shi, Tan, JCP 131, 094502 (2009)
!use, intrinsic :: iso_fortran_env
use omp_lib
implicit none

!! precision
!integer, parameter :: wp1 = real64
integer,parameter::wp1=selected_real_kind(15,37)
integer,parameter::wp2=selected_real_kind(16)

!! Numbers
real(kind=wp1),parameter::rzero=0.0_wp1
complex(kind=wp1),parameter::czero=complex(rzero,rzero)
real(kind=wp1),parameter::half=0.5_wp1
real(kind=wp1),parameter::one=1.0_wp1
real(kind=wp1),parameter::two=2.0_wp1
real(kind=wp1),parameter::three=3.0_wp1
real(kind=wp1),parameter::four=4.0_wp1
real(kind=wp1),parameter::obs=one/6.0_wp1
complex(kind=wp1),parameter :: ima = complex(rzero,one)
real(kind=wp1), parameter ::pi=four*atan(one)
real(kind=wp1), parameter :: tolerance = 1.0e-6

!! Parameters
real(kind=wp1), parameter :: hbar = 5308.8_wp1!!in cm-1 fs
complex(kind=wp1), parameter :: ih = ima/hbar
real(kind=wp1), parameter :: kb = 0.69503_wp1!! in cm-1 K
real(kind=wp1), parameter :: speed_of_light = 299792458.0_wp1 * 100.0_wp1 * 1.0E-15 !! cm/fs
real(kind=wp1), parameter :: wc_to_tau = one/speed_of_light/(pi + pi)
integer, parameter :: nlevel = 1
!! if model == 1
! integer, parameter :: model = 1
! real(kind=wp1), parameter :: wc = 50.0_wp1
! real(kind=wp1), parameter :: delta = wc
! real(kind=wp1), parameter :: e11 = wc
! real(kind=wp1), parameter :: eta = wc / 8.0_wp1 * two

! real(kind=wp1), parameter :: tau_c = wc_to_tau/wc
! real(kind=wp1), parameter :: temperature = 4.0 * hbar /tau_c / kb
! integer, parameter :: KK = 0, LL = 20
! real(kind=wp1), parameter :: tend = 20.0_wp1/delta * hbar
! real(kind=wp1), parameter :: dt = one !! fs
!!!!!!!!!!!!!!!
!! if model == 2
! integer, parameter :: model = 2
! real(kind=wp1), parameter :: wc = 50.0_wp1
! real(kind=wp1), parameter :: delta = wc
! real(kind=wp1), parameter :: e11 = wc
! real(kind=wp1), parameter :: eta = wc / 8.0_wp1 * two

! real(kind=wp1), parameter :: tau_c = wc_to_tau/wc
! real(kind=wp1), parameter :: temperature = 0.20 * hbar /tau_c / kb
! real(kind=wp1), parameter :: tend = 20.0_wp1/delta * hbar
! real(kind=wp1), parameter :: dt = one/10.0_wp1 !! fs
! integer, parameter :: KK = 3, LL = 25
!!!!!!!!!!!!!!!!
!! if model == 3
integer, parameter :: model = 3
real(kind=wp1), parameter :: wc = 50.0_wp1
real(kind=wp1), parameter :: delta =  wc
real(kind=wp1), parameter :: e11 = 1.5_wp1* wc
real(kind=wp1), parameter :: eta = wc * 0.75_wp1 * two

real(kind=wp1), parameter :: tau_c = wc_to_tau/wc !! 1/fs
real(kind=wp1), parameter :: temperature = one/0.25_wp1 * hbar /tau_c / kb !! K
integer, parameter :: KK = 2, LL = 25
real(kind=wp1), parameter :: tend = 15.0_wp1/delta * hbar !! fs
real(kind=wp1), parameter :: dt = one/10.0_wp1 !! fs
!!!!!!!!!!!!!!!!
! !! if model == 4
! integer, parameter :: model = 4
! real(kind=wp1), parameter :: wc = 50.0_wp1
! real(kind=wp1), parameter :: delta = wc
! real(kind=wp1), parameter :: e11 = 3.0_wp1 * wc
! real(kind=wp1), parameter :: eta = wc * 0.375_wp1 * two

! real(kind=wp1), parameter :: tau_c = wc_to_tau/wc
! real(kind=wp1), parameter :: temperature = one/0.25_wp1 * hbar /tau_c / kb
! real(kind=wp1), parameter :: tend = 40.0_wp1/delta * hbar
! real(kind=wp1), parameter :: dt = one/10.0_wp1 !! fs
! integer, parameter :: KK = 2, LL = 25


real(kind=wp1), parameter ::  omega_c = one/tau_c !!in 1/fs
!! HEOM variables
integer :: nn_tot
integer,allocatable :: nn(:,:,:)
integer,allocatable :: map_nplus(:,:,:),map_nneg(:,:,:),map_sum(:)
integer,allocatable :: zero(:)
complex(kind=wp1),allocatable :: rho(:,:,:),rho_t0(:,:)
complex(kind=wp1),allocatable :: mu_p_tot(:,:,:),mu_p_t0(:,:),mu_vib_p(:,:)
complex(kind=wp1),allocatable :: mu_m_tot(:,:,:),mu_m_t0(:,:),mu_vib_m(:,:)
complex(kind=wp1),allocatable :: c_matsubara(:),omg_matsubara(:)
real(kind=wp1) :: sum_c_over_nu

!! Response function
complex(kind=wp1),allocatable :: dip_mom_corr(:)

!! Response function
complex(kind=wp1),allocatable :: rup(:,:), rdn(:,:), nup(:,:), ndn(:,:)

!! System specific
integer :: nlevel2
complex(kind=wp1),allocatable :: Hamil_e(:,:),dip_moment(:,:),dip_vib(:,:),dip_tot(:,:,:)
complex(kind=wp1),allocatable :: dip_full(:,:)
real(kind=wp1),allocatable :: Q_op(:,:,:)

!! Evolution
integer :: nsteps, t1steps, t2steps, t3steps
real(kind=wp1) :: curr_time
logical, parameter:: printfile=.true.!!.true.!!.true.
logical, parameter:: initialCheck=.true.!!.false.
logical, parameter:: bath =.true.
character(len=100) :: nStr, lamStr, tempStr, tauStr, filename, file_end_st


!! Misc
real(kind=wp1) :: var1, mu12, mu13, y12, j23
logical, parameter:: flag_spectra=.true.!!.false.!!calculate population/spectra
complex(kind=wp1),allocatable :: mat_mu_all(:,:,:)
integer :: ierror
contains
!----------------------------------------------------------
!----------------------------------------------------------

subroutine setup
  implicit none

  if (model >4 .or. model <= 0) then
    print*, 'Invalid model identifier'
    stop
  end if

  print*, "Number of Hierarchy levels : ", LL
  print*, "Number of Matsubara terms : ", KK
  print*, "Total propagation time ", tend, "fs"
  print*, "Time step ", dt, "fs"

  !----------------------------------------------------------
  nlevel2 = 2


  call compute_nn_tot
  write(6,*) "number of auxiliary matrices:",nn_tot
  nsteps = nint(tend/dt) + 1
  print*,"number of levels is", nlevel
  print*,"length of the Hamilonian is", nlevel2
  print*,"time step dt is:", dt

  allocate(dip_mom_corr(nsteps),dip_moment(nlevel2,nlevel2),dip_vib(nlevel2,nlevel2),dip_tot(nlevel2,nlevel2,nn_tot))
  allocate(nn(nlevel,0:KK,nn_tot),map_nplus(nlevel,0:KK,nn_tot),map_nneg(nlevel,0:KK,nn_tot))
  allocate(map_sum(0:LL),zero(nn_tot))
  allocate(rho(nlevel2,nlevel2,nn_tot),rho_t0(nlevel2,nlevel2))
  allocate(mu_p_tot(nlevel2,nlevel2,nn_tot),mu_p_t0(nlevel2,nlevel2),mu_vib_p(nlevel2,nlevel2))
  allocate(mu_m_tot(nlevel2,nlevel2,nn_tot),mu_m_t0(nlevel2,nlevel2),mu_vib_m(nlevel2,nlevel2))
  allocate(c_matsubara(0:KK),omg_matsubara(0:KK))
  allocate(Hamil_e(nlevel2,nlevel2))
  allocate(Q_op(nlevel2,nlevel2,nlevel),mat_mu_all(nlevel2,nlevel2,nsteps))

  call compute_nn
  call compute_map

end subroutine setup
!----------------------------------------------------------

subroutine main
  implicit none

  call setup_parameters
  call init_cond
  call computePop
  ! call linearSpectra

end subroutine main
!----------------------------------------------------------

subroutine setup_parameters
  implicit none
  integer :: k
  character(len=100) :: filename
  real(kind=wp1), dimension(:,:), allocatable :: hamil_r

  allocate(hamil_r(nlevel2, nlevel2))
  omg_matsubara(0)=omega_c
  c_matsubara(0)=eta*omega_c*half* (one/tan(hbar*omega_c/(two*kb*temperature))-ima)
  do k=1,KK
    omg_matsubara(k)=two*k*pi*kb*temperature/hbar
    c_matsubara(k)=two*eta*omega_c*kb*temperature/hbar * omg_matsubara(k)/(omg_matsubara(k)**2-omega_c**2)
  end do

  sum_c_over_nu=real(sum(c_matsubara(:)/omg_matsubara(:)))/hbar
  var1=eta*kb*temperature/(omega_c*hbar*hbar)

  mu12=one
  mu13=-0.2_wp1
  dip_moment=rzero

  !!Add electronic elements
  dip_moment(1,2) = mu12
  dip_moment(2,1) = mu12
  ! dip_moment(1,3) = mu13
  ! dip_moment(3,1) = mu13
 
  filename = 'hamil.txt'
  if (access(filename, 'r') == 0) then
      print*, 'Reading Hamiltonian from txt file'
      open(11, file=filename)
      do k = 1,nlevel2
          read(11, * ) hamil_r(k,:)
      end do
      read(11,*) file_end_st
      close(11)
      if (file_end_st .ne. 'x') then
          write(*,*) "problem in reading input file hamil.txt"
          stop
      end if
  else
      print*, 'Using user-defined (2 by 2) Hamiltonian'
      hamil_r = rzero
      hamil_r(1,1) = e11
      hamil_r(2,2) = -e11
      hamil_r(1,2) = delta
      hamil_r(2,1) = delta
  end if

  call print2d(hamil_r)
  Hamil_e = hamil_r + czero
  deallocate(hamil_r)

  Q_op=rzero
  Q_op(1,1,1) = one
  Q_op(2,2,1) = -one
  ! do k=1,nlevel2
  !   Q_op(k,k,k) = one
  ! end do

  zero=0

write(6,*) "Parameters Set ..."

end subroutine setup_parameters
!-----------------------------------------------------------------

subroutine init_cond
  implicit none
  !integer :: i

  rho=rzero
  if(flag_spectra.eqv..true.) then
    rho(1,1,1)=one
    dip_tot=rzero
    dip_tot(:,:,1)=matmul(dip_moment,rho(:,:,1))
    rho_t0=rho(:,:,1)
  else if(flag_spectra.eqv..false.) then
    rho(2,2,1)=one
    rho_t0=rho(:,:,1)
  end if

  curr_time = rzero
  dip_mom_corr = rzero


write(6,*) "Intitial Conditions set ... "

end subroutine init_cond
!-----------------------------------------------------------------

subroutine computePop
  implicit none
  integer :: i, j
  complex(kind=wp1),allocatable :: re_t1(:,:,:,:), dummy1(:,:,:)

  allocate(dummy1(nlevel2,nlevel2,nn_tot))
  allocate(re_t1(nlevel2,nlevel2,nn_tot,nsteps))

  dummy1 = czero
  dummy1(1,1,1) = one !! rho(:,:,1)
  do i=1,nsteps
    re_t1(:,:,:,i) = dummy1
    if(mod(i,1)==1)call filter
    call rk4(dummy1, dt)
    if (mod(i,(nsteps/20)).eq.0) then
      write(*,'(A,F16.8,A)')"done with = ",i*dt , " fs"
      call flush()
    end if
    if (real(re_t1(1,1,1,i)) > two) then
      print*, 'check propagation scheme, density matrix is diverging'
      stop
    end if
  end do

  print*,"Now writing population file"
  nStr = str(model)

  ! filename = "SB_HEOM_"//trim(nStr)//"LS_lambda_"//trim(lamStr)//"_T_"// &
  !               trim(tempStr)//"_tau_"//trim(tauStr)//".txt"
  filename = "SB_model"//trim(nStr)//"_heom.txt"
  open(104, file = filename, status = "replace")
  do i=1,nsteps
    write(104,'(f15.7$)')real(i-1)*dt*delta/hbar
    do j=1,nlevel2
      write(104,'(f15.7$)')real(re_t1(j,j,1,i)) 
    end do
    write(104,'(f15.7$)')real(re_t1(1,2,1,i))
    write(104,'(f15.7$)')aimag(re_t1(1,2,1,i))
    if (nlevel2 > 2) then
      write(104,'(f15.6$)')real(re_t1(2,3,1,i))
      write(104,'(f15.6$)')aimag(re_t1(2,3,1,i))
    end if 
    write(104,*)
  end do
  close(104)

  deallocate(dummy1, re_t1)

end subroutine computePop
!-----------------------------------------------------------------

function densityProp(in1,in2,fi1,fi2,time)
  implicit none
  integer :: i
  integer, intent(in)::in1,in2,fi1,fi2
  integer, intent(in)::time
  complex(kind=wp1),dimension(time) :: densityProp

  rho=czero
  rho(in1,in2,1)=one

  print*,'Propagating with initial condition in ', in1, in2, 'to', fi1, fi2
  do i=1,time,1
    densityProp(i) = rho(fi1,fi2,1)
    if(mod(i,10)==1)call filter
    call rk4(rho,dt)
    if (mod(i,(time/20)).eq.0) then
      if ((fi1 == fi2) .and. (in1 == in2)) then
        print*,'Population : done with step = ', i, 'out of ', time
      end if
      if ((fi1 /= fi2) .or. (in1 /= in2)) then
        print*,'Coherence : done with step = ', i, 'out of ', time
      end if
    end if
  end do
end function densityProp
!----------------------------------------------------------------- 

subroutine linearSpectra
  implicit none
  integer :: i
  complex(kind=wp1),allocatable :: re_t1(:,:,:,:), dummy1(:,:,:), res(:)

  allocate(dummy1(nlevel2,nlevel2,nn_tot))
  allocate(re_t1(nlevel2,nlevel2,nn_tot,t1steps))
  allocate(res(t1steps))

  mu_p_t0=tril(dip_moment)
  mu_m_t0=triu(dip_moment)

  !!!!!!!!!!   1.rephasing  up-down contribution
  dummy1=czero
  dummy1(:,:,1) = matmul(mu_p_t0, rho_t0)
  !t1 propagation
  do i=1,nsteps
    re_t1(:,:,:,i)=dummy1
    if(mod(i,1)==1)call filter
    call rk4(dummy1,dt)
  end do

  do i=1,nsteps
    res(i) = ctrace( matmul(mu_m_t0, re_t1(:,:,1,i)) )
  end do

  !!if(printfile.eqv..true.) then
  print*,"Now writing 1D spectra file"
  open(102, file = "heom_linres.out", status = "replace")
  do i=1,nsteps
    write(102,'(f15.7$)')real(i-1)*dt/hbar * delta
    write(102,'(f15.7$)')real(res(i))
    write(102,'(f15.7$)')aimag(res(i))
    write(102,*)
  end do
  close(102)
  !!end if

  deallocate(re_t1, dummy1, res)

end subroutine linearSpectra
!-----------------------------------------------------------------

subroutine rk4(rho,dt)
  implicit none
  complex(kind=wp1),intent(inout)::rho(nlevel2,nlevel2,nn_tot)
  real(kind=wp1),intent(in)::dt
  complex(kind=wp1),dimension(nlevel2,nlevel2,nn_tot) :: k1,k2,k3,k4


  call compute_deriv(rho,k1)
  call compute_deriv(rho+half*dt*k1,k2)
  call compute_deriv(rho+half*dt*k2,k3)
  call compute_deriv(rho+dt*k3,k4)

  rho = rho + dt*obs*(k1+two*(k2+k3)+k4)

end subroutine rk4
!-----------------------------------------------------------------

subroutine compute_deriv(rho,drho_dt)
  implicit none
  complex(kind=wp1),intent(in)::rho(nlevel2,nlevel2,nn_tot)
  complex(kind=wp1),intent(out)::drho_dt(nlevel2,nlevel2,nn_tot)
  complex(kind=wp1) :: tmp(nlevel2,nlevel2)
  integer :: n

  !$omp parallel default(shared) private(tmp)
  !$omp do schedule(static)
  do n=1,nn_tot
    call compute_deriv_n(n,rho,tmp)
    drho_dt(:,:,n)=tmp
  end do
  !$omp end do nowait
  !$omp end parallel

end subroutine compute_deriv
!-----------------------------------------------------------------

subroutine compute_deriv_n(n,rho,drho_dt)
  implicit none
  !! Eq. 15, JCP 131, 094502 (2009)
  integer,intent(in) :: n
  complex(kind=wp1),intent(in)::rho(nlevel2,nlevel2,nn_tot)
  complex(kind=wp1),intent(out)::drho_dt(nlevel2,nlevel2)
  complex(kind=wp1) :: tmp(nlevel2,nlevel2),mat_tmp(nlevel2,nlevel2)
  integer :: m,k,nplus,nminus
  integer :: nvec(nlevel,0:KK)

  m=n
  call index(m,nvec,0)
  !!!!part 1.0
  tmp=-ih*commute(Hamil_e,rho(:,:,n))

  if(zero(n)==0) then !! matrix at n is not filtered out

    !!!!part 1.5
    do m=1,nlevel
      do k=0,KK
        tmp=tmp-nvec(m,k)*omg_matsubara(k)*rho(:,:,n)
      end do
    end do

    !!!!part 2.0
    do m=1,nlevel
      mat_tmp=Q_op(:,:,m) !! matrix |m><m|
      tmp=tmp-(var1-sum_c_over_nu)*commute(mat_tmp,commute(mat_tmp,rho(:,:,n)))
    end do
  end if

  !!!!part 3.0
  do m=1,nlevel
    mat_tmp=Q_op(:,:,m) !! matrix |m><m|
    do k=0,KK
      nplus=map_nplus(m,k,n)

      if(nplus>0.and.nplus<=nn_tot) then
        if(zero(nplus)==0) tmp=tmp - ima*commute(mat_tmp,rho(:,:,nplus))* sqrt((nvec(m,k)+one)*abs(c_matsubara(k)))
      end if
    end do
  end do

  !!!!part 4.0
  do m=1,nlevel
    mat_tmp=Q_op(:,:,m)!! matrix |m><m|
    do k=0,KK
      nminus=map_nneg(m,k,n)
      if(nminus>0.and.nminus<=nn_tot) then
        if(zero(nminus)==0)then
          tmp=tmp-ih*(c_matsubara(k)*matmul(mat_tmp,rho(:,:,nminus))&
                -conjg(c_matsubara(k))*matmul(rho(:,:,nminus),mat_tmp) )*sqrt(nvec(m,k)/abs(c_matsubara(k)))
        end if
      end if
    end do
  end do

  drho_dt=tmp

end subroutine compute_deriv_n
!-----------------------------------------------------------------

subroutine filter
  !! Shi, Chen, Nan, Xu, Yan, JCP 130, 084105 (2009)
  implicit none
  integer :: n

  zero=0
  do n=1,nn_tot
    if(maxval(abs(rho(:,:,n)))<tolerance) then
      rho(:,:,n)=rzero
      zero(n)=one
    end if
  end do

end subroutine filter
!-----------------------------------------------------------------

subroutine index(n,nvec,iflag)
  implicit none
  integer,intent(inout)::n
  integer,intent(inout)::nvec(nlevel,0:KK)
  integer,intent(in)::iflag
  integer :: m, n_beg, n_end, l_sum

  if(iflag==0) then  !n is input, nvec is output
    nvec=nn(:,:,n)
  endif

  if(iflag==1) then  !nvec is input, n is output
    n=0
    l_sum=sum(nvec)
    if(l_sum<=LL) then
      n_beg=map_sum(l_sum)
      if(l_sum==LL) n_end=nn_tot
      if(l_sum<LL) n_end=map_sum(l_sum+1)-1
      do m=1,nn_tot
        if(all(nvec==nn(:,:,m))) then
          n=m
          exit
        end if
      end do
    end if
  end if

end subroutine index

!-----------------------------------------------------------------
subroutine compute_nn
  implicit none
  integer :: n, n_tot, n_beg, n_end

  n_tot=0
  n_beg=0;n_end=0
  do n=0,LL
    call compute_nn_sum_L(n,n_beg,n_end)
    map_sum(n)=n_beg
    write(6,*) n,n_beg,n_end
  end do
  write(6,*) n_end,nn_tot

end subroutine compute_nn
!-----------------------------------------------------------------

subroutine compute_nn_sum_L(L,n_beg,n_end)
  implicit none
  integer,intent(in)::L
  integer,intent(inout)::n_beg,n_end!! at input, state/end point of entries for sum L-1, at output, for entries for sum L
  integer :: i,j,m,k,tot_L
  integer :: flag_new,nvec(nlevel,0:KK)

  if(L==0) then
    nn(:,:,1)=0
    n_beg=1;n_end=1
  else
    tot_L=n_end
    do i=n_beg,n_end
      do m=1,nlevel
        do k=0,KK

          nvec=nn(:,:,i);nvec(m,k)=nvec(m,k)+1
          flag_new=0
          if(tot_L>n_end) then
            do j=n_end+1,tot_L
              if(all(nvec==nn(:,:,j))) then
                flag_new=1;exit
              end if
            end do
          end if

          if(flag_new==0) then
            tot_L=tot_L+1
            nn(:,:,tot_L)=nn(:,:,i)
            nn(m,k,tot_L)=nn(m,k,i)+1
          end if
        end do
      end do
    end do
    n_beg=n_end+1
    n_end=tot_L
  end if

end subroutine compute_nn_sum_L
!-----------------------------------------------------------------

subroutine compute_map
  implicit none
  integer :: n,m,k
  integer :: nvec(nlevel,0:KK),nvec_plus(nlevel,0:KK),nvec_neg(nlevel,0:KK)
  integer :: nplus,nneg

  map_nneg=0
  map_nplus=0
  do n=1,nn_tot
    m=n
    call index(m,nvec,0)
    do m=1,nlevel
      do k=0,KK
        nvec_plus=nvec;nvec_plus(m,k)=nvec_plus(m,k)+1
        call index(nplus,nvec_plus,1)
        map_nplus(m,k,n)=nplus

        if(nvec(m,k)>0) then
          nvec_neg=nvec;nvec_neg(m,k)=nvec_neg(m,k)-1
          call index(nneg,nvec_neg,1)
          map_nneg(m,k,n)=nneg
        end if
      end do
    end do
  end do

write(6,*) "map completed ..."

end subroutine compute_map
!-----------------------------------------------------------------

subroutine compute_nn_tot
  implicit none
  integer :: i
  real(kind=wp1) :: tmp

  tmp=one
  do i=1,nlevel*(KK+1)
    tmp=tmp*(LL+i)/real(i,kind=wp1)
  end do
  nn_tot=nint(tmp)

end subroutine compute_nn_tot
!-----------------------------------------------------------------

function commute(mat1,mat2) result(mat3)
  implicit none
  complex(kind=wp1),intent(in) :: mat1(:,:),mat2(:,:)
  complex(kind=wp1),allocatable :: mat3(:,:)
  integer :: i1

  i1=size(mat1,1)
  allocate(mat3(i1,i1))

  mat3 = matmul(mat1,mat2) - matmul(mat2,mat1)

end function commute
!-----------------------------------------------------------------

function ctrace(dum_arr) result(arr)
  implicit none
  complex(kind=wp1),intent(in) :: dum_arr(:,:)
  complex(kind=wp1) :: arr
  integer :: i1, i

  i1=size(dum_arr,1)

  arr=czero
  do i=1,i1
    arr=arr+dum_arr(i,i)
  end do

end function ctrace
!-----------------------------------------------------------------

function triu(matrix)
  implicit none
  complex(kind=wp1), dimension(:,:) :: matrix
  complex(kind=wp1), dimension(:,:),allocatable::triu
  integer::i1,j,k

  i1=size(matrix,1)

  allocate(triu(i1,i1))
  triu=rzero

  do j=1,i1
      do k=1,i1
          if(j.le.k) then
              triu(j,k) = matrix(j,k)
          else
              cycle
          end if
      end do
  end do
end function triu
!-----------------------------------------------------------------

function tril(matrix)
  implicit none
  complex(kind=wp1), dimension(:,:) :: matrix
  complex(kind=wp1), dimension(:,:),allocatable::tril
  integer::i1,j,k

  i1=size(matrix,1)

  allocate(tril(i1,i1))
  tril=rzero

  do j=1,i1
  do k=1,i1
      if(j.ge.k) then
          tril(j,k) = matrix(j,k)
      else
          cycle
      end if
  end do
  end do
end function tril
!-----------------------------------------------------------------

subroutine get_walltime(wctime)
  use iso_fortran_env
  implicit none
  real(kind=wp1) :: wctime
  integer :: r, c
  call system_clock(c, r)
    wctime = real(c,kind=wp1) / r
end subroutine get_walltime
!-----------------------------------------------------------------

character(len=20) function str(k)
  implicit none
  !"Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
!-----------------------------------------------------------------

subroutine print2d(array)
  implicit none
  real(kind=wp1), dimension(:,:), intent(in) :: array
  integer :: len1, li
  len1 = size(array,1)
  do li=1,len1
      write (*,'(*(f12.6))') array(li,:)
  end do
end subroutine print2d
!-----------------------------------------------------------------

end Module mod_spectra_heom4
