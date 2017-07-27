!****************************************************************
module input_output
! Global input and output routines.
!****************************************************************

  use shared
 
  implicit none

contains

!****************************************************************
  subroutine readDims(filename,nr,nth,nph)
    
    use netcdf

    character*(*), intent(IN) :: filename
    integer, intent(INOUT) :: nr,nth,nph
    integer :: ncid,dimid
    
    call check( nf90_open(filename, NF90_NOWRITE, ncid) )

    call check( nf90_inq_dimid(ncid, "r", dimid) )
    call check( nf90_inquire_dimension(ncid, dimid, len = nr) )
    call check( nf90_inq_dimid(ncid, "th", dimid) )
    call check( nf90_inquire_dimension(ncid, dimid, len = nth) )
    call check( nf90_inq_dimid(ncid, "ph", dimid) )
    call check( nf90_inquire_dimension(ncid, dimid, len = nph) )

    call check( nf90_close(ncid) )

  end subroutine readDims

  
!****************************************************************
  subroutine readCoords(filename,r,th,ph,y)
    
    use netcdf

    character*(*), intent(IN) :: filename
    integer :: ncid,varid
    real(dp), dimension(:), intent(INOUT) :: r,th,ph
    integer, dimension(:), intent(INOUT) :: y
    real(dp) :: dEq
    
    call check( nf90_open(filename, NF90_NOWRITE, ncid) )

    call check( nf90_inq_varid(ncid, "r", varid) )
    call check( nf90_get_var(ncid, varid, r) )

    call check( nf90_inq_varid(ncid, "th", varid) )
    call check( nf90_get_var(ncid, varid, th) )

    call check( nf90_inq_varid(ncid, "ph", varid) )
    call check( nf90_get_var(ncid, varid, ph) )
       
    call check( nf90_close(ncid) )
  
    dEq = ph(2) - ph(1)
    y = nint(-log(tan(0.5_dp*th))/dEq)
  
  end subroutine readCoords

  
!****************************************************************
  subroutine uniformTheta(th,ph,th1)
      
    real(dp), dimension(:), intent(INOUT) :: th, ph, th1
    integer :: j, ymin
    real(dp) :: dEq
    
    dEq = ph(2) - ph(1)
    ymin = nint(-log(tan(0.5_dp*th(1)))/dEq)
    do j=1,size(th1)
        th1(j) = 2.0_dp*atan(exp(-dEq*real(j-1+ymin)))
    end do  

  end subroutine uniformTheta
  

!****************************************************************
  subroutine readB(filename,br,bth,bph,y)
    !
    ! This version maps B to a uniform grid in y.
    
    use netcdf

    character*(*), intent(IN) :: filename
    integer :: ncid,varid,nr,nph,j,nth1,yprv,ymin,ymax,j1,j2
    real(dp), dimension(:,:,:), intent(INOUT) :: br,bth,bph
    integer, dimension(:), intent(IN) :: y
    real(dp) :: frc
    
    nr = size(br,1)
    nph = size(br,3)
    
    call check( nf90_open(filename, NF90_NOWRITE, ncid) )

    call check( nf90_inq_varid(ncid, "br", varid) )
    do j=1,size(y)
        call check( nf90_get_var(ncid, varid, br(:,y(j) - y(1) + 1,:),start=(/1, j, 1/), &
                count=(/nr,1,nph/)) )
    end do
    call check( nf90_inq_varid(ncid, "bth", varid) )
    do j=1,size(y)
        call check( nf90_get_var(ncid, varid, bth(:,y(j) - y(1) + 1,:),start=(/1, j, 1/), &
                count=(/nr,1,nph/)) )
    end do
    call check( nf90_inq_varid(ncid, "bph", varid) )
    do j=1,size(y)
        call check( nf90_get_var(ncid, varid, bph(:,y(j) - y(1) + 1,:),start=(/1, j, 1/), &
                count=(/nr,1,nph/)) )
    end do
        
    call check( nf90_close(ncid) )

    ! Interpolate intermediate values in y:
    ymin = minval(y)
    ymax = maxval(y)
    nth1 = ymax - ymin + 1
    yprv = 1
    do j=ymin,ymax-1
        frc = real(j - y(yprv))/real(y(yprv+1) - y(yprv))
        j1 = y(yprv)-ymin+1
        j2 = y(yprv+1)-ymin+1
        br(:,j-ymin+1,:) = (1.0_dp - frc)*br(:,j1,:) + frc*br(:,j2,:)
        bth(:,j-ymin+1,:) = (1.0_dp - frc)*bth(:,j1,:) + frc*bth(:,j2,:)
        bph(:,j-ymin+1,:) = (1.0_dp - frc)*bph(:,j1,:) + frc*bph(:,j2,:)  
        if (j+1.eq.y(yprv+1)) yprv = yprv+1
    end do
    
  end subroutine readB

 
  !****************************************************************
  subroutine addA2file(filename,r,th,ph,ar,ath,aph,y)
    
    use netcdf

    character*(*), intent(in) :: filename
    integer :: ncid,r_varid,th_varid,ph_varid,ar_varid,ath_varid,aph_varid
    integer :: nth_dimid,nph_dimid,nr_dimid
    integer :: nth, nph, nr, istatus, j
    real(dp), dimension(:), intent(in) :: r, th, ph
    integer, dimension(:), intent(in) :: y
    real(dp), dimension(:,:,:), intent(in) :: ar, ath, aph

    nr = size(r)
    nth = size(th)
    nph = size(ph)
    
    call check(nf90_open(filename, nf90_write, ncid))
    call check(nf90_inq_dimid(ncid, 'r', nr_dimid))
    call check(nf90_inq_dimid(ncid, 'th', nth_dimid))
    call check(nf90_inq_dimid(ncid, 'ph', nph_dimid))
    ! Define variables ar, ath, aph if they don't already exist,
    ! otherwise get their variable id's and overwrite them.
    istatus = nf90_inq_varid(ncid, 'ar', ar_varid)
    if (istatus.ne.nf90_noerr) then
        call check(nf90_redef(ncid))
        call check(nf90_def_var(ncid,'ar',NF90_DOUBLE, &
            (/nr_dimid,nth_dimid,nph_dimid/),ar_varid))
        call check(nf90_def_var(ncid,'ath',NF90_DOUBLE, &
            (/nr_dimid,nth_dimid,nph_dimid/),ath_varid))
        call check(nf90_def_var(ncid,'aph',NF90_DOUBLE, &
            (/nr_dimid,nth_dimid,nph_dimid/),aph_varid))        
        call check(nf90_enddef(ncid))
    else
        call check(nf90_inq_varid(ncid, 'ath', ath_varid))
        call check(nf90_inq_varid(ncid, 'aph', aph_varid))
    end if
        
    ! Restrict back to original grid in y-direction:
    do j=1,size(y)
        call check( nf90_put_var(ncid, ar_varid, ar(:,y(j) - y(1) + 1,:),start=(/1, j, 1/), &
                count=(/nr,1,nph/)) )
        call check( nf90_put_var(ncid, ath_varid, ath(:,y(j) - y(1) + 1,:),start=(/1, j, 1/), &
                count=(/nr,1,nph/)) )
        call check( nf90_put_var(ncid, aph_varid, aph(:,y(j) - y(1) + 1,:),start=(/1, j, 1/), &
                count=(/nr,1,nph/)) )              
    end do      
        
    call check( nf90_close(ncid) )
   
  end subroutine addA2file

!****************************************************************
  subroutine writeA2file(filename,r,th,ph,ar,ath,aph,y)

    use netcdf

    character*(*), intent(in) :: filename
    integer :: ncid,r_varid,th_varid,ph_varid,ar_varid,ath_varid,aph_varid
    integer :: nth_dimid,nph_dimid,nr_dimid
    integer :: nth, nph, nr, istatus, j
    real(dp), dimension(:), intent(in) :: r, th, ph
    integer, dimension(:), intent(in) :: y
    real(dp), dimension(:,:,:), intent(in) :: ar, ath, aph

    nr = size(r)
    nth = size(th)
    nph = size(ph)

    call check(nf90_create(filename, nf90_write, ncid))
    ! Check for defined dimensions
    istatus = nf90_inq_dimid(ncid, 'r', nr_dimid)
    if (istatus.ne.nf90_noerr) then
        call check(nf90_def_dim(ncid, 'r', nr, nr_dimid)) 
        call check(nf90_def_dim(ncid, 'th', nth, nth_dimid)) 
        call check(nf90_def_dim(ncid, 'ph', nph, nph_dimid)) 
    else
        call check(nf90_inq_dimid(ncid, 'th', nth_dimid))
        call check(nf90_inq_dimid(ncid, 'ph', nph_dimid))
    endif

    ! Define variables ar, ath, aph if they don't already exist,
    ! otherwise get their variable id's and overwrite them.
    istatus = nf90_inq_varid(ncid, 'ar', ar_varid)
    if (istatus.ne.nf90_noerr) then
        call check(nf90_def_var(ncid,'ar',NF90_DOUBLE, &
            (/nr_dimid,nth_dimid,nph_dimid/),ar_varid))
        call check(nf90_def_var(ncid,'ath',NF90_DOUBLE, &
            (/nr_dimid,nth_dimid,nph_dimid/),ath_varid))
        call check(nf90_def_var(ncid,'aph',NF90_DOUBLE, &
            (/nr_dimid,nth_dimid,nph_dimid/),aph_varid))        
        call check(nf90_enddef(ncid))
    else
        call check(nf90_inq_varid(ncid, 'ath', ath_varid))
        call check(nf90_inq_varid(ncid, 'aph', aph_varid))
    end if

    ! Restrict back to original grid in y-direction:
    do j=1,size(y)
        call check( nf90_put_var(ncid, ar_varid, ar(:,y(j) - y(1) + 1,:),start=(/1, j, 1/), &
                count=(/nr,1,nph/)) )
        call check( nf90_put_var(ncid, ath_varid, ath(:,y(j) - y(1) + 1,:),start=(/1, j, 1/), &
                count=(/nr,1,nph/)) )
        call check( nf90_put_var(ncid, aph_varid, aph(:,y(j) - y(1) + 1,:),start=(/1, j, 1/), &
                count=(/nr,1,nph/)) )              
    end do      
        
    call check( nf90_close(ncid) )

  end subroutine writeA2file

!****************************************************************
  subroutine check(istatus)
!----------------------------------------------------------------
! Check (ever so slightly modified from www.unidata.ucar.edu).
! For netcdf.
!----------------------------------------------------------------
    use netcdf
    integer, intent (IN) :: istatus
    
    if (istatus /= nf90_noerr) then
       write(*,*) trim(adjustl(nf90_strerror(istatus)))
    end if
    
  end subroutine check

end module input_output
