program spectra_heom
  use mod_spectra_heom4
  use omp_lib
  implicit none
  real*8::start,finish
  integer::numproc
  

  numproc = omp_get_num_procs()
  print*,"with processors=",numproc
  call get_walltime(start)

  call setup
  call main
  
  call get_walltime(finish)
  write(*,*) "wall time     : ",(finish-start),"seconds"

end program spectra_heom
