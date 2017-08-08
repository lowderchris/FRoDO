! Fast field-line tracer.
!
! -- A. Yeates 2/10/14
!
! To compile into a python module, use:
! > f2py -c --fcompiler=gnu95 -m tracer tracer.f90
!   /opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin
!   gfortran?
!   install via homebrew...
! To compile with support for openMP:
! > f2py -m tracer --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c tracer.f90
! Then it can be called in python by "import tracer".
! To adjust the number of thread cores:
! setenv OMP_NUM_THREADS 4
!
! Note: on Macbook Air I had to run
! > xcode-select --install
! before this would compile.

subroutine fieldline(br, bth, bph, r, th, ph, r0, th0, ph0, &
     nmax, maxError, minB, nr, nth, nph, nfl0, rl, thl, phl)
  implicit none
  
  integer, intent(in) :: nr
  integer, intent(in) :: nth
  integer, intent(in) :: nph
  integer, intent(in) :: nfl0
  double precision, intent(in), dimension(nr,nth,nph) :: br
  double precision, intent(in), dimension(nr,nth,nph) :: bth
  double precision, intent(in), dimension(nr,nth,nph) :: bph
  double precision, intent(in), dimension(nr) :: r
  double precision, intent(in), dimension(nth) :: th
  double precision, intent(in), dimension(nph) :: ph
  double precision, intent(in), dimension(nfl0) :: r0
  double precision, intent(in), dimension(nfl0) :: th0
  double precision, intent(in), dimension(nfl0) :: ph0
  integer, intent(in) :: nmax
  double precision, intent(in) :: maxError
  double precision, intent(in) :: minB
  double precision, intent(out), dimension(nfl0,nmax) :: rl
  double precision, intent(out), dimension(nfl0,nmax) :: thl
  double precision, intent(out), dimension(nfl0,nmax) :: phl
  double precision :: rMin, rMax, thMin, thMax
  double precision, parameter :: TWOPI = 8.*datan(1.0d0)
  double precision :: maxds, ds, r1, th1, ph1, r2, th2, ph2
  double precision :: k1r, k1th, k1ph, k2r, k2th, k2ph, error
  double precision :: dr1, dr2, dth1, dth2, dph1, dph2, ds_dt
  integer :: nxt, dirn, i
  double precision :: ddirn
  !integer :: tid, OMP_GET_THREAD_NUM, nthr, OMP_GET_NUM_THREADS
  !real :: tstart, tfinish

  ! Extents of domain:
  rMin = minval(r)
  rMax = maxval(r)
  thMin = minval(th)
  thMax = maxval(th)

  ! Maximum allowed step-size:
  maxds = 0.05d0
  
  ! Initialise variables:
  rl(:,:) = -1.0d0
  
  ! Note: removed the private variables tid, nthr, tstart, and tfinish from
  ! below

  ! Main loop:
  !$OMP PARALLEL DO &
  !$OMP& private(ds,r1,th1,ph1,r2,th2,ph2,&
  !$OMP& k1r,k1th,k1ph,k2r,k2th,k2ph,error,&
  !$OMP& dr1,dr2,dth1,dth2,dph1,dph2,ds_dt,nxt,dirn,ddirn)
  do i=1,nfl0
     !tid = OMP_GET_THREAD_NUM()
     !nthr = OMP_GET_NUM_THREADS()
     do dirn=-1,1,2
        ddirn = dble(dirn)
        ds=maxds
        nxt = nmax/2

        ! Add startpoint to output array:
        rl(i,nxt) = r0(i)
        thl(i,nxt) = th0(i)
        phl(i,nxt) = ph0(i)

        ! Interpolate k1:
        call qkinterp(br, bth, bph, r, th, ph, &
             rl(i,nxt), thl(i,nxt), phl(i,nxt), &
             k1r, k1th, k1ph, nr, nth, nph)      
        ds_dt = dsqrt(k1r**2 + k1th**2 + k1ph**2)
        if (dabs(ds_dt) < minB) exit   ! stop if null is reached 
        if (nxt > nmax) exit     ! stop if beyond nmax
        if (nxt < 1) exit        ! stop if beyond nax
        ds_dt = ds_dt * ddirn
        k1r = k1r/ds_dt
        k1th = k1th/ds_dt/rl(i,nxt)
        k1ph = k1ph/ds_dt/rl(i,nxt)/dsin(thl(i,nxt))

        do
           ! Compute midpoint:
           r2 = rl(i,nxt) + ds*k1r
           th2 = thl(i,nxt) + ds*k1th
           ph2 = mod(phl(i,nxt) + ds*k1ph + TWOPI, TWOPI)

           ! If outside boundary, do Euler step to boundary and stop:
           if ((r2.lt.rMin).or.(r2.gt.rMax).or. &
                (th2.lt.thMin).or.(th2.gt.thMax)) then
              if ((r2.lt.rMin).or.(r2.gt.rMax)) then
                 if (k1r.lt.0.0d0) then
                    ds = (rMin - rl(i,nxt))/k1r
                 else
                    ds = (rMax - rl(i,nxt))/k1r
                 end if
              else
                 ds = 1.0d6
              end if
              if ((th2.lt.thMin).or.(th2.gt.thMax)) then
                 if (k1th.lt.0.0d0) then
                    ds = min(ds,(thMin - thl(i,nxt))/k1th)
                 else
                    ds = min(ds,(thMax - thl(i,nxt))/k1th)
                 end if
              end if
              nxt = nxt + dirn
              if (nxt > nmax) exit     ! stop if beyond nmax
              if (nxt < 1) exit        ! stop if beyond nax
              rl(i,nxt) = rl(i,nxt-dirn) + ds*k1r
              thl(i,nxt) = thl(i,nxt-dirn) + ds*k1th
              phl(i,nxt) = mod(phl(i,nxt-dirn) + ds*k1ph + TWOPI, TWOPI)
              exit
           end if

           ! Interpolate k2:
           call qkinterp(br, bth, bph, r, th, ph, r2, th2, ph2, &
                k2r, k2th, k2ph, nr, nth, nph)
           
           ds_dt = dsqrt(k2r**2 + k2th**2 + k2ph**2)
           if (dabs(ds_dt) < minB) exit   ! stop if null is reached
           if (nxt > nmax) exit     ! stop if beyond nmax
           if (nxt < 1) exit        ! stop if beyond nax
           ds_dt = ds_dt * ddirn
 
           ! Compute first and second-order update:
           k2r = k2r/ds_dt
           k2th = k2th/ds_dt/r2
           k2ph = k2ph/ds_dt/r2/dsin(th2)
           dr1 = ds*k1r
           dth1 = ds*k1th
           dph1 = ds*k1ph
           dr2 = 0.5d0*ds*(k1r + k2r)
           dth2 = 0.5d0*ds*(k1th + k2th)
           dph2 = 0.5d0*ds*(k1ph + k2ph)

           ! Estimate error from difference:
           error = (dr1-dr2)**2 + rl(i,nxt)*(dth2-dth1)**2 &
                + rl(i,nxt)*dsin(thl(i,nxt))*(dph2-dph1)**2

           ! Modify step size depending on error:
           if (error.eq.0.0d0) then
              ds = maxds
           else
              ds = min(maxds, 0.85d0*dabs(ds)*(maxerror/error)**0.25d0)
           end if

           ! Update if error is small enough:
           if (error.le.maxerror) then
              r1 = rl(i,nxt) + dr2
              th1 = thl(i,nxt) + dth2
              ph1 = mod(phl(i,nxt) + dph2 + TWOPI, TWOPI)
              ! Return midpoint if full step leaves domain:
              if ((r1.lt.rMin).or.(r1.gt.rMax).or. &
                   (th1.lt.thMin).or.(th1.gt.thMax)) then
                 r1 = r2
                 th1 = th2
                 ph1 = ph2
              end if
              nxt = nxt + dirn
              if (nxt > nmax) exit     ! stop if beyond nmax
              if (nxt < 1) exit        ! stop if beyond nax
              rl(i,nxt) = r1
              thl(i,nxt) = th1
              phl(i,nxt) = ph1

              ! Interpolate k1 at next point:
              call qkinterp(br, bth, bph, r, th, ph, r1, th1, ph1, &
                   k1r, k1th, k1ph, nr, nth, nph)      
              ds_dt = dsqrt(k1r**2 + k1th**2 + k1ph**2)
              if (dabs(ds_dt) < minB) exit   ! stop if null is reached 
              if (nxt > nmax) exit     ! stop if beyond nmax
              if (nxt < 1) exit        ! stop if beyond nax
              ds_dt = ds_dt * ddirn
              k1r = k1r/ds_dt
              k1th = k1th/ds_dt/r1
              k1ph = k1ph/ds_dt/r1/dsin(th1)
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine fieldline


subroutine fieldline_cubic(br, bth, bph, r, th, ph, r0, th0, ph0, &
     nmax, maxError, minB, nr, nth, nph, nfl0, rl, thl, phl)
  implicit none
  
  integer, intent(in) :: nr
  integer, intent(in) :: nth
  integer, intent(in) :: nph
  integer, intent(in) :: nfl0
  double precision, intent(in), dimension(nr,nth,nph) :: br
  double precision, intent(in), dimension(nr,nth,nph) :: bth
  double precision, intent(in), dimension(nr,nth,nph) :: bph
  double precision, intent(in), dimension(nr) :: r
  double precision, intent(in), dimension(nth) :: th
  double precision, intent(in), dimension(nph) :: ph
  double precision, intent(in), dimension(nfl0) :: r0
  double precision, intent(in), dimension(nfl0) :: th0
  double precision, intent(in), dimension(nfl0) :: ph0
  integer, intent(in) :: nmax
  double precision, intent(in) :: maxError
  double precision, intent(in) :: minB
  double precision, intent(out), dimension(nfl0,nmax) :: rl
  double precision, intent(out), dimension(nfl0,nmax) :: thl
  double precision, intent(out), dimension(nfl0,nmax) :: phl
  double precision :: rMin, rMax, thMin, thMax
  double precision, parameter :: TWOPI = 8.*datan(1.0d0)
  double precision :: maxds, ds, r1, th1, ph1, r2, th2, ph2
  double precision :: k1r, k1th, k1ph, k2r, k2th, k2ph, error
  double precision :: dr1, dr2, dth1, dth2, dph1, dph2, ds_dt
  integer :: nxt, dirn, i

  double precision, dimension(64,64) :: m
  double precision, dimension(8,nr,nth,nph,3) :: db

  ! Prepare tricubic interpolation:
  call tricubicPrepare(br, bth, bph, db, m, nr, nth, nph)

  ! Extents of domain:
  rMin = minval(r)
  rMax = maxval(r)
  thMin = minval(th)
  thMax = maxval(th)

  ! Maximum allowed step-size:
  maxds = 0.05d0
  
  ! Initialise variables:
  rl(:,:) = -1.0d0

  ! Main loop:
  do i=1,nfl0     
     do dirn=-1,1,2
        ! Reverse direction of field for backward and forward tracing:
        db = -db

        ds=maxds
        nxt = nmax/2

        ! Add startpoint to output array:
        rl(i,nxt) = r0(i)
        thl(i,nxt) = th0(i)
        phl(i,nxt) = ph0(i)

        ! Interpolate k1:
        call tricubicInterp(m, db, r, th, ph, &
             rl(i,nxt), thl(i,nxt), phl(i,nxt), &
             k1r, k1th, k1ph, nr, nth, nph)      
        ds_dt = dsqrt(k1r**2 + k1th**2 + k1ph**2)
        if (ds_dt < minB) exit   ! stop if null is reached 
        k1r = k1r/ds_dt
        k1th = k1th/ds_dt/rl(i,nxt)
        k1ph = k1ph/ds_dt/rl(i,nxt)/dsin(thl(i,nxt))

        do
           ! Compute midpoint:
           r2 = rl(i,nxt) + ds*k1r
           th2 = thl(i,nxt) + ds*k1th
           ph2 = mod(phl(i,nxt) + ds*k1ph + TWOPI, TWOPI)

           ! If outside boundary, do Euler step to boundary and stop:
           if ((r2.lt.rMin).or.(r2.gt.rMax).or. &
                (th2.lt.thMin).or.(th2.gt.thMax)) then
              if ((r2.lt.rMin).or.(r2.gt.rMax)) then
                 if (k1r.lt.0.0d0) then
                    ds = (rMin - rl(i,nxt))/k1r
                 else
                    ds = (rMax - rl(i,nxt))/k1r
                 end if
              else
                 ds = 1.0d6
              end if
              if ((th2.lt.thMin).or.(th2.gt.thMax)) then
                 if (k1th.lt.0.0d0) then
                    ds = min(ds,(thMin - thl(i,nxt))/k1th)
                 else
                    ds = min(ds,(thMax - thl(i,nxt))/k1th)
                 end if
              end if
              nxt = nxt + dirn
              rl(i,nxt) = rl(i,nxt-dirn) + ds*k1r
              thl(i,nxt) = thl(i,nxt-dirn) + ds*k1th
              phl(i,nxt) = mod(phl(i,nxt-dirn) + ds*k1ph + TWOPI, TWOPI)
              exit
           end if

           ! Interpolate k2:
           call tricubicInterp(m, db, r, th, ph, r2, th2, ph2, &
                k2r, k2th, k2ph, nr, nth, nph)
           
           ds_dt = dsqrt(k2r**2 + k2th**2 + k2ph**2)
           if (ds_dt < minB) exit   ! stop if null is reached
           
           ! Compute first and second-order update:
           k2r = k2r/ds_dt
           k2th = k2th/ds_dt/r2
           k2ph = k2ph/ds_dt/r2/dsin(th2)
           dr1 = ds*k1r
           dth1 = ds*k1th
           dph1 = ds*k1ph
           dr2 = 0.5d0*ds*(k1r + k2r)
           dth2 = 0.5d0*ds*(k1th + k2th)
           dph2 = 0.5d0*ds*(k1ph + k2ph)

           ! Estimate error from difference:
           error = (dr1-dr2)**2 + rl(i,nxt)*(dth2-dth1)**2 &
                + rl(i,nxt)*dsin(thl(i,nxt))*(dph2-dph1)**2

           ! Modify step size depending on error:
           if (error.eq.0.0d0) then
              ds = maxds
           else
              ds = min(maxds, 0.85d0*dabs(ds)*(maxerror/error)**0.25d0)
           end if

           ! Update if error is small enough:
           if (error.le.maxerror) then
              r1 = rl(i,nxt) + dr2
              th1 = thl(i,nxt) + dth2
              ph1 = mod(phl(i,nxt) + dph2 + TWOPI, TWOPI)
              ! Return midpoint if full step leaves domain:
              if ((r1.lt.rMin).or.(r1.gt.rMax).or. &
                   (th1.lt.thMin).or.(th1.gt.thMax)) then
                 r1 = r2
                 th1 = th2
                 ph1 = ph2
              end if
              nxt = nxt + dirn
              rl(i,nxt) = r1
              thl(i,nxt) = th1
              phl(i,nxt) = ph1

              ! Interpolate k1 at next point:
              call tricubicInterp(m, db, r, th, ph, r1, th1, ph1, &
                   k1r, k1th, k1ph, nr, nth, nph)      
              ds_dt = dsqrt(k1r**2 + k1th**2 + k1ph**2)
              if (ds_dt < minB) exit   ! stop if null is reached 
              k1r = k1r/ds_dt
              k1th = k1th/ds_dt/r1
              k1ph = k1ph/ds_dt/r1/dsin(th1)
           end if
        end do
     end do
  end do
  
end subroutine fieldline_cubic


subroutine qkinterp(br, bth, bph, r, th, ph, r1, th1, ph1, &
     br1, bth1, bph1, nr, nth, nph)
  ! Fast trilinear interpolation for vector on regular grid
  implicit none

  integer, intent(in) :: nr
  integer, intent(in) :: nth
  integer, intent(in) :: nph
  double precision, intent(in), dimension(nr,nth,nph) :: br
  double precision, intent(in), dimension(nr,nth,nph) :: bth
  double precision, intent(in), dimension(nr,nth,nph) :: bph
  double precision, intent(in), dimension(nr) :: r
  double precision, intent(in), dimension(nth) :: th
  double precision, intent(in), dimension(nph) :: ph
  double precision, intent(in) :: r1
  double precision, intent(in) :: th1
  double precision, intent(in) :: ph1
  double precision, intent(out) :: br1
  double precision, intent(out) :: bth1
  double precision, intent(out) :: bph1

  double precision :: dr, dth, dph, fr, fth, fph
  double precision :: a000, a001, a010, a011, a100, a101, a110, a111
  integer :: ir, ith, iph

  dr = r(2) - r(1)
  ir = floor((r1 - r(1))/dr) + 1
  fr = (r1 - r(ir))/dr
  dth = th(2) - th(1)
  ith = floor((th1 - th(1))/dth) + 1
  fth = (th1 - th(ith))/dth
  dph = ph(2) - ph(1)
  iph = floor((ph1 - ph(1))/dph) + 1
  fph = (ph1 - ph(iph))/dph

  a000 = (1.0d0-fr)*(1.0d0-fth)*(1.0d0-fph)
  a001 = (1.0d0-fr)*(1.0d0-fth)*fph
  a010 = (1.0d0-fr)*fth*(1.0d0-fph)
  a011 = (1.0d0-fr)*fth*fph
  a100 = fr*(1.0d0-fth)*(1.0d0-fph)
  a101 = fr*(1.0d0-fth)*fph
  a110 = fr*fth*(1.0d0-fph)
  a111 = fr*fth*fph

  if ((ir.eq.size(r)).and.(fr.eq.0.0d0)) then
     ir = ir-1
     fr = 1.0d0
  end if
  if ((ith.eq.size(th)).and.(fth.eq.0.0d0)) then
     ith = ith-1
     fth = 1.0d0
  end if

  br1 = a000*br(ir,ith,iph) + a100*br(ir+1,ith,iph) &
       + a001*br(ir,ith,iph+1) + a101*br(ir+1,ith,iph+1) &
       + a010*br(ir,ith+1,iph) + a110*br(ir+1,ith+1,iph) &
       + a011*br(ir,ith+1,iph+1) + a111*br(ir+1,ith+1,iph+1)

  bth1 = a000*bth(ir,ith,iph) + a100*bth(ir+1,ith,iph) &
       + a001*bth(ir,ith,iph+1) + a101*bth(ir+1,ith,iph+1) &
       + a010*bth(ir,ith+1,iph) + a110*bth(ir+1,ith+1,iph) &
       + a011*bth(ir,ith+1,iph+1) + a111*bth(ir+1,ith+1,iph+1)

  bph1 = a000*bph(ir,ith,iph) + a100*bph(ir+1,ith,iph) &
       + a001*bph(ir,ith,iph+1) + a101*bph(ir+1,ith,iph+1) &
       + a010*bph(ir,ith+1,iph) + a110*bph(ir+1,ith+1,iph) &
       + a011*bph(ir,ith+1,iph+1) + a111*bph(ir+1,ith+1,iph+1)

end subroutine qkinterp


subroutine tricubicPrepare(br, bth, bph, db, m, nr, nth, nph)
  ! Precompute derivatives and read in matrix for tricubic interpolation.
  implicit none

  integer, intent(in) :: nr
  integer, intent(in) :: nth
  integer, intent(in) :: nph
  double precision, intent(in), dimension(nr,nth,nph) :: br
  double precision, intent(in), dimension(nr,nth,nph) :: bth
  double precision, intent(in), dimension(nr,nth,nph) :: bph
  double precision, dimension(64,64), intent(out) :: m
  double precision, dimension(8,nr,nth,nph,3), intent(out) :: db
  integer :: i,j

  ! (1) Estimate derivatives at grid points.
  ! ----------------------------------------
  db = 0.0d0
  ! Function itself:
  db(1,:,:,:,1) = br
  db(1,:,:,:,2) = bth
  db(1,:,:,:,3) = bph
  ! d/dr: (1-sided for boundaries)
  db(2,2:nr-1,:,:,:) = 0.5d0*(db(1,3:nr,:,:,:) - db(1,1:nr-2,:,:,:))
  db(2,1,:,:,:) = 0.5d0*(-3.0d0*db(1,1,:,:,:) + 4.0d0*db(1,2,:,:,:) - db(1,3,:,:,:))
  db(2,nr,:,:,:) = 0.5d0*(3.0d0*db(1,nr,:,:,:) - 4.0d0*db(1,nr-1,:,:,:) + db(1,nr-2,:,:,:))
  ! d/dth: (1-sided for boundaries)
  db(3,:,2:nth-1,:,:) = 0.5d0*(db(1,:,3:nth,:,:) - db(1,:,1:nth-2,:,:)) 
  db(3,:,1,:,:) = 0.5d0*(-3.0d0*db(1,:,1,:,:) + 4.0d0*db(1,:,2,:,:) - db(1,:,3,:,:))
  db(3,:,nth,:,:) = 0.5d0*(3.0d0*db(1,:,nth,:,:)-4.0d0*db(1,:,nth-1,:,:) + db(1,:,nth-2,:,:))
  ! d/dph: (periodic for boundaries)
  db(4,:,:,2:nph-1,:) = 0.5d0*(db(1,:,:,3:nph,:) - db(1,:,:,1:nph-2,:))
  db(4,:,:,1,:) = 0.5d0*(db(1,:,:,2,:) - db(1,:,:,nph-1,:))
  db(4,:,:,nph,:) = 0.5d0*(db(1,:,:,2,:) - db(1,:,:,nph-1,:))
  ! d^2/drdth: (1-sided on boundaries)
  db(5,2:nr-1,:,:,:) = 0.5d0*(db(3,3:nr,:,:,:) - db(3,1:nr-2,:,:,:))
  db(5,1,:,:,:) = 0.5d0*(-3.0d0*db(3,1,:,:,:) + 4.0d0*db(3,2,:,:,:) - db(3,3,:,:,:))
  db(5,nr,:,:,:) = 0.5d0*(3.0d0*db(3,nr,:,:,:) - 4.0d0*db(3,nr-1,:,:,:) + db(3,nr-2,:,:,:))
  ! d^2/drdph: (1-sided on boundaries)
  db(6,2:nr-1,:,:,:) = 0.5d0*(db(4,3:nr,:,:,:) - db(4,1:nr-2,:,:,:))
  db(6,1,:,:,:) = 0.5d0*(-3.0d0*db(4,1,:,:,:) + 4.0d0*db(4,2,:,:,:) - db(4,3,:,:,:))
  db(6,nr,:,:,:) = 0.5d0*(3.0d0*db(4,nr,:,:,:) - 4.0d0*db(4,nr-1,:,:,:) + db(4,nr-2,:,:,:))
  ! d^2/dthdph: (1-sided on boundaries)
  db(7,:,2:nth-1,:,:) = 0.5d0*(db(4,:,3:nth,:,:) - db(4,:,1:nth-2,:,:))
  db(7,:,1,:,:) = 0.5d0*(-3.0d0*db(4,:,1,:,:) + 4.0d0*db(4,:,2,:,:) - db(4,:,3,:,:))
  db(7,:,nr,:,:) = 0.5d0*(3.0d0*db(4,:,nth,:,:)- 4.0d0*db(4,:,nth-1,:,:) + db(4,:,nth-2,:,:))
  ! d^3/drdthdph: (1-sided on boundaries)
  db(8,2:nr-1,:,:,:) = 0.5d0*(db(7,3:nr,:,:,:) - db(7,1:nr-2,:,:,:))
  db(8,1,:,:,:) = 0.5d0*(-3.0d0*db(7,1,:,:,:) + 4.0d0*db(7,2,:,:,:) - db(7,3,:,:,:))
  db(8,nr,:,:,:) = 0.5d0*(3.0d0*db(7,nr,:,:,:) - 4.0d0*db(7,nr-1,:,:,:) + db(7,nr-2,:,:,:))

  ! (2) Read in (inverse) interpolation matrix.
  ! -------------------------------------------
  open(1,file='inverse.txt')
  do i = 1,64
      read(1,*) (m(j,i),j=1,64)
  end do
  close(1)

end subroutine tricubicPrepare


subroutine tricubicInterp(m, db, r, th, ph, r1, th1, ph1, &
     br1, bth1, bph1, nr, nth, nph)
  ! C1 tricubic interpolation for vector on regular grid
  ! - using method of Lekien & Marsden 2005.
  implicit none

  integer, intent(in) :: nr
  integer, intent(in) :: nth
  integer, intent(in) :: nph
  double precision, intent(in), dimension(nr) :: r
  double precision, intent(in), dimension(nth) :: th
  double precision, intent(in), dimension(nph) :: ph
  double precision, intent(in) :: r1
  double precision, intent(in) :: th1
  double precision, intent(in) :: ph1
  double precision, intent(out) :: br1
  double precision, intent(out) :: bth1
  double precision, intent(out) :: bph1

  double precision :: dr, dth, dph, fr, fth, fph
  integer :: ir, ith, iph

  double precision, dimension(64,64), intent(in) :: m
  double precision, dimension(3) :: al
  double precision :: fact
  double precision, dimension(8,nr,nth,nph,3), intent(in) :: db
  integer :: i,j

  ! Identify required cell.
  ! -----------------------
  dr = r(2) - r(1)
  ir = floor((r1 - r(1))/dr) + 1
  fr = (r1 - r(ir))/dr
  dth = th(2) - th(1)
  ith = floor((th1 - th(1))/dth) + 1
  fth = (th1 - th(ith))/dth
  dph = ph(2) - ph(1)
  iph = floor((ph1 - ph(1))/dph) + 1
  fph = (ph1 - ph(iph))/dph

  if ((ir.eq.size(r)).and.(fr.eq.0.0d0)) then
     ir = ir-1
     fr = 1.0d0
  end if
  if ((ith.eq.size(th)).and.(fth.eq.0.0d0)) then
     ith = ith-1
     fth = 1.0d0
  end if

  ! Compute interpolation coefficients.
  ! -----------------------------------
  ! Loop through each combination of powers of r, th, ph
  br1 = 0.0d0
  bth1 = 0.0d0
  bph1 = 0.0d0
  do i = 1,64
     al = 0.0d0
     do j = 0,7
        al = al + m(i,8*j+1)*db(j+1,ir,ith,iph,:) &
             + m(i,8*j+2)*db(j+1,ir+1,ith,iph,:) &
             + m(i,8*j+3)*db(j+1,ir,ith+1,iph,:) &
             + m(i,8*j+4)*db(j+1,ir+1,ith+1,iph,:) &
             + m(i,8*j+5)*db(j+1,ir,ith,iph+1,:) &
             + m(i,8*j+6)*db(j+1,ir+1,ith,iph+1,:) &
             + m(i,8*j+7)*db(j+1,ir,ith+1,iph+1,:) &
             + m(i,8*j+8)*db(j+1,ir+1,ith+1,iph+1,:)
     end do
     fact = fr**mod(i-1,4)*fth**mod((i-1)/4,4)*fph**((i-1)/16)
     br1 = br1 + al(1)*fact
     bth1 = bth1 + al(2)*fact
     bph1 = bph1 + al(3)*fact
  end do

end subroutine tricubicInterp






subroutine jsq(jj, rl, thl, phl, r, th, ph, jsq1, nfl0, nmax, nr, nth, nph)
  implicit none
  
  integer, intent(in) :: nr
  integer, intent(in) :: nth
  integer, intent(in) :: nph
  integer, intent(in) :: nfl0
  integer, intent(in) :: nmax
  double precision, intent(in), dimension(nr,nth,nph) :: jj
  double precision, intent(in), dimension(nfl0,nmax) :: rl
  double precision, intent(in), dimension(nfl0,nmax) :: thl
  double precision, intent(in), dimension(nfl0,nmax) :: phl
  double precision, intent(in), dimension(nr) :: r
  double precision, intent(in), dimension(nth) :: th
  double precision, intent(in), dimension(nph) :: ph
  double precision, intent(out), dimension(nfl0) :: jsq1
  integer :: i,j,k
  double precision :: jj1,dd,smax
  double precision, dimension(nmax) :: s,jjs

  ! Loop through each fieldline:
  do i=1,nfl0     
     !if ((i/100)*100.eq.i) print*,i,' of ',nfl0
     
     ! Integrate value of jj along this fieldline...
     ! (a) Go along fieldline and store jj and arclength:
     k=0
     jjs=0.0d0
     jj1=0.0d0
     s=0.0d0
     do j=1,nmax
        if (rl(i,j).gt.0.0d0) then
           call qkinterp_jj(jj, r, th, ph, rl(i,j), thl(i,j), phl(i,j), &
                jj1, nr, nth, nph)
           if (k.gt.0) then
              dd=rl(i,j)**2 + rl(i,j-1)**2 - 2.0d0*rl(i,j)*rl(i,j-1)* &
                   (sin(thl(i,j))*sin(thl(i,j-1))*cos(phl(i,j)-phl(i,j-1)) &
                   + cos(thl(i,j))*cos(thl(i,j-1))) + 1.0d-12
              if (dd.lt.0.0d0) print*,dd
              s(j)=s(j-1) + sqrt(dd)
           end if
           if (.not.isnan(jj1)) jjs(j)=jjs(j) + jj1
           k=k+1
        end if
     end do
     ! (b) Integrate using trapezium rule:
     jsq1(i)=0.0d0
     smax=0.0d0
     do j=1,nmax
        if (rl(i,j).gt.0.0d0) then
           jsq1(i)=jsq1(i) + 0.5d0*(s(j)-s(j-1))*(jjs(j) + jjs(j-1))
           smax=s(j)
        end if
     end do
     ! Divide by field-line length:
     jsq1(i)=jsq1(i)/smax
  !   if (smax.lt.10.0d0) then
  !      jsq1(i)=0.0d0
  !   else
  !      jsq1(i)=1.0d0
  !   end if
  !   print*,smax,jsq1(i)

  end do
     
end subroutine jsq

subroutine njsq(jj, rl, thl, phl, r, th, ph, jsq1, nfl0, nmax, nr, nth, nph)
  implicit none
  
  integer, intent(in) :: nr
  integer, intent(in) :: nth
  integer, intent(in) :: nph
  integer, intent(in) :: nfl0
  integer, intent(in) :: nmax
  double precision, intent(in), dimension(nr,nth,nph) :: jj
  double precision, intent(in), dimension(nfl0,nmax) :: rl
  double precision, intent(in), dimension(nfl0,nmax) :: thl
  double precision, intent(in), dimension(nfl0,nmax) :: phl
  double precision, intent(in), dimension(nr) :: r
  double precision, intent(in), dimension(nth) :: th
  double precision, intent(in), dimension(nph) :: ph
  double precision, intent(out), dimension(nfl0) :: jsq1
  integer :: i,j,k
  double precision :: jj1,dd
  double precision, dimension(nmax) :: s,jjs

  ! Loop through each fieldline:
  do i=1,nfl0 
     !if ((i/100)*100.eq.i) print*,i,' of ',nfl0
     ! Integrate value of jj along this fieldline...
     ! (a) Go along fieldline and store jj and arclength:
     k=0
     jjs=0.0d0
     jj1=0.0d0
     s=0.0d0
     do j=1,nmax
        if (rl(i,j).gt.0.0d0) then
           call qkinterp_jj(jj, r, th, ph, rl(i,j), thl(i,j), phl(i,j), &
                jj1, nr, nth, nph)
           if (k.gt.0) then
              dd=rl(i,j)**2 + rl(i,j-1)**2 - 2.0d0*rl(i,j)*rl(i,j-1)* &
                   (sin(thl(i,j))*sin(thl(i,j-1))*cos(phl(i,j)-phl(i,j-1)) &
                   + cos(thl(i,j))*cos(thl(i,j-1))) + 1.0d-12
              if (dd.lt.0.0d0) print*,dd
              s(j)=s(j-1) + sqrt(dd)
           end if
           if (.not.isnan(jj1)) jjs(j)=jjs(j) + jj1
           k=k+1
        end if
     end do
     ! (b) Integrate using trapezium rule:
     jsq1(i)=0.0d0
     do j=1,nmax
        if (rl(i,j).gt.0.0d0) then
           jsq1(i)=jsq1(i) + 0.5d0*(s(j)-s(j-1))*(jjs(j) + jjs(j-1))
        end if
     end do

  end do
     
end subroutine njsq

subroutine qkinterp_jj(jj, r, th, ph, r1, th1, ph1, jj1, nr, nth, nph)
  ! Fast trilinear interpolation for vector on regular grid
  implicit none

  integer, intent(in) :: nr
  integer, intent(in) :: nth
  integer, intent(in) :: nph
  double precision, intent(in), dimension(nr,nth,nph) :: jj
  double precision, intent(in), dimension(nr) :: r
  double precision, intent(in), dimension(nth) :: th
  double precision, intent(in), dimension(nph) :: ph
  double precision, intent(in) :: r1
  double precision, intent(in) :: th1
  double precision, intent(in) :: ph1
  double precision, intent(out) :: jj1
  double precision :: dr, dth, dph, fr, fth, fph
  double precision :: a000, a001, a010, a011, a100, a101, a110, a111
  integer :: ir, ith, iph

  dr = r(2) - r(1)
  ir = floor((r1 - r(1))/dr) + 1
  fr = (r1 - r(ir))/dr
  dth = th(2) - th(1)
  ith = floor((th1 - th(1))/dth) + 1
  fth = (th1 - th(ith))/dth
  dph = ph(2) - ph(1)
  iph = floor((ph1 - ph(1))/dph) + 1
  fph = (ph1 - ph(iph))/dph

  a000 = (1.0d0-fr)*(1.0d0-fth)*(1.0d0-fph)
  a001 = (1.0d0-fr)*(1.0d0-fth)*fph
  a010 = (1.0d0-fr)*fth*(1.0d0-fph)
  a011 = (1.0d0-fr)*fth*fph
  a100 = fr*(1.0d0-fth)*(1.0d0-fph)
  a101 = fr*(1.0d0-fth)*fph
  a110 = fr*fth*(1.0d0-fph)
  a111 = fr*fth*fph

  if ((ir.eq.size(r)).and.(fr.eq.0.0d0)) then
     ir = ir-1
     fr = 1.0d0
  end if
  if ((ith.eq.size(th)).and.(fth.eq.0.0d0)) then
     ith = ith-1
     fth = 1.0d0
  end if

  jj1 = a000*jj(ir,ith,iph) + a100*jj(ir+1,ith,iph) &
       + a001*jj(ir,ith,iph+1) + a101*jj(ir+1,ith,iph+1) &
       + a010*jj(ir,ith+1,iph) + a110*jj(ir+1,ith+1,iph) &
       + a011*jj(ir,ith+1,iph+1) + a111*jj(ir+1,ith+1,iph+1)

end subroutine qkinterp_jj
