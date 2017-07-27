!****************************************************************
module vectorpot
! Computing the vector potential.
!****************************************************************

  use shared
  
  implicit none

contains


 !****************************************************************
  subroutine computeAGreens(r, th, ph, br, bth, bph, ar, ath, aph)
    !
    ! Compute the line integral vector potential, with the z=0 terms
    ! found using the Greens function solution.
    !
    ! WARNING: This method gives a continuous A but leads to errors in B if
    ! computed by curling with centered differences.
    !
    real(dp), intent(in), dimension(:) :: r, th, ph
    real(dp), intent(in), dimension(:,:,:) :: br, bth, bph
    real(dp), intent(inout), dimension(:,:,:) :: ar, ath, aph
    real(dp), dimension(:,:,:), allocatable :: hx, hz
    real(dp), dimension(:,:), allocatable :: br1,f,th2,ph2,psi
    real(dp), dimension(:), allocatable :: thr, phr
    double precision :: dEq
    integer :: nr, nth, nph, i, j, yMin
    real(dp), parameter :: SXTNPI = 64.0_dp*atan(1.0_dp)
    
    nr = size(r)
    nth = size(th)
    nph = size(ph)

    ! Create scale factors
    dEq = ph(2) - ph(1)
    allocate(hx(nr,nth,nph), hz(nr,nth,nph))
    do j=1,nth
       do i=1,nr
          hx(i,j,:) = dEq*r(i)*sin(th(j))
          ! Note: hy is the same as hx so save storage
       end do
    end do
    do i=1,nr
       hz(i,:,:) = dEq*r(i)
    end do
    
    ! Compute A by trapezium rule, taking advantage of equal x,y,z coords
    ar = 0.0_dp
    ath = 0.0_dp
    aph = 0.0_dp
    do i=2,nr
       aph(i,:,:) = aph(i-1,:,:) - &
            0.5_dp*(hx(i-1,:,:)*hz(i-1,:,:)*bth(i-1,:,:) &
            + hx(i,:,:)*hz(i,:,:)*bth(i,:,:))
       ! note: minus sign in ax comes from by = -bth
       ath(i,:,:) = ath(i-1,:,:) + &
            0.5_dp*(hx(i-1,:,:)*hz(i-1,:,:)*bph(i-1,:,:) &
            + hx(i,:,:)*hz(i,:,:)*bph(i,:,:))
    end do

    ! Surface terms      A = grad(psi) x e_r
    ! -------------
    allocate(th2(1:nth,1:nph), ph2(1:nth,1:nph), br1(1:nth,1:nph),f(1:nth,1:nph))
    br1 = hx(1,:,:)*hx(1,:,:)*br(1,:,:)
    do i=1,nth
       th2(i,:) = th(i)
    end do
    do i=1,nph
       ph2(:,i) = ph(i)
    end do
    allocate(thr(1:nth+1), phr(1:nph+1))
    yMin = nint(-log(tan(0.5_dp*maxval(th)))/dEq)
    do i=1,nth+1
       thr(i) = 2.0_dp*atan(exp(-dEq*(real(i-1+yMin)-0.5_dp)))
    end do
    do i=1,nph+1
       phr(i) = (real(i-1)-0.5_dp)*dEq
    end do
    !
    ! (1) Compute psi on ph-ribs
    allocate(psi(1:nth,1:nph+1))
    !$omp parallel private(j,f)
    !$omp do
    do i=1,nth
       do j=1,nph+1
          f = 1.0_dp -cos(th2)*cos(th(i)) -sin(th2)*sin(th(i))*cos(ph2 - phr(j))
          f = log(f)*br1
          psi(i,j) = -sum(f(1:nth-1,1:nph-1) + f(1:nth-1,2:nph) &
               + f(2:nth,1:nph-1) + f(2:nth,2:nph))/SXTNPI
       end do
    end do
    !$omp end do
    !$omp end parallel

    ! (2) Compute Ath at grid points
    !$omp parallel do
    do i=1,nr
       ath(i,1:nth,1:nph) = ath(i,1:nth,1:nph) + &
            psi(1:nth,2:nph+1) - psi(1:nth,1:nph)
    end do
    !$omp end parallel do
    deallocate(psi)
    !
    ! (3) Compute psi on th-ribs
    allocate(psi(1:nth+1,1:nph))
    !$omp parallel private(j,f)
    !$omp do
    do i=1,nth+1
       do j=1,nph
          f= 1.0_dp -cos(th2)*cos(thr(i)) -sin(th2)*sin(thr(i))*cos(ph2 - ph(j))
          f = log(f)*br1
          psi(i,j) = -sum(f(1:nth-1,1:nph-1) + f(1:nth-1,2:nph) &
               + f(2:nth,1:nph-1) + f(2:nth,2:nph))/SXTNPI
       end do
    end do
    !$omp end do
    !$omp end parallel
    !
    ! (4) Compute Aph at grid points
    !$omp parallel do 
    do i=1,nr
       aph(i,1:nth,1:nph) = aph(i,1:nth,1:nph) + &
            psi(2:nth+1,1:nph) - psi(1:nth,1:nph)
    end do
    !$omp end parallel do
    deallocate(psi)

    aph = aph/hx
    ath = ath/hx
                    
    deallocate(hx,hz,f,th2,ph2,thr,phr,br1)

  end subroutine computeAGreens

    
end module vectorpot
