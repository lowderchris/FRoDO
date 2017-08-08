!****************************************************************
!program main
!----------------------------------------------------------------
! Main program: compute vector potential for glovag code
! To change number of threads, before execution:
! > setenv OMP_NUM_THREADS 4
!----------------------------------------------------------------

subroutine compA(snap, filePath)

  use shared
  use input_output
  use vectorpot
  use omp_lib
  
  implicit none

  character*(*) :: filePath
  character*(*) :: snap
  integer :: nr,nth,nph,nth1
  real(dp), dimension(:,:,:), allocatable :: br,bth,bph,ar,ath,aph
  real(dp), dimension(:), allocatable :: r,th,ph,th1
  integer, dimension(:), allocatable :: y
  real(dp) :: dEq
  integer :: threads, thread1
 
  !$OMP PARALLEL
  threads = omp_get_num_threads()
  thread1 = omp_get_thread_num()
  if (thread1.eq.0) print*,'Number of OpenMp threads:',threads
  !$OMP END PARALLEL

  print*,'--',snap,'--'

  ! Read coordinates:
  call readDims(filePath//'b_'//snap//'.nc',nr,nth,nph)
  allocate(r(nr),th(nth),ph(nph),y(nth))
  call readCoords(filePath//'b_'//snap//'.nc',r,th,ph,y)
  
  ! Make uniform array in y:
  dEq = ph(2)-ph(1)
  nth1 = nint(-log(tan(0.5_dp*th(nth)))/dEq) -  nint(-log(tan(0.5_dp*th(1)))/dEq) + 1
  allocate(th1(nth1))
  call uniformTheta(th,ph,th1)
  
  ! Allocate arrays (uniform spacing in y):
  allocate(ar(nr,nth1,nph),ath(nr,nth1,nph),aph(nr,nth1,nph))
  allocate(br(nr,nth1,nph),bth(nr,nth1,nph),bph(nr,nth1,nph))
  
  ! Read in B:
  call readB(filePath//'b_'//snap//'.nc',br,bth,bph,y)
  
  ! Compute vector potential:
  !call computeAGreens(r, th1, ph, br, bth, bph, ar, ath, aph)
  call computeAgreens(r, th1, ph, br, bth, bph, ar, ath, aph)

  ! Add vector potential to output file:
  call writeA2file(filePath//'a_'//snap//'.nc',r,th,ph,ar,ath,aph,y)
  deallocate(r,th,ph,y)
  deallocate(th1)
  deallocate(ar,ath,aph)
  deallocate(br,bth,bph)
  !end do
  
end subroutine compA

!end program main
