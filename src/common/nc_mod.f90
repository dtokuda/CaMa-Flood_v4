module nc_mod
#ifdef UseCDF_CMF
    use PARKIND1, only: &
    &   JPIM, JPIB, JPRM, JPRB, JPRD
    use YOS_CMF_INPUT, only: &
    &   LOGNAM, LLEAPYR
    use datetime_mod, only: &
    &   DateTime
    use cmf_cf_time_mod, only: &
    &   cf_time_axis, &
    &   cf_read_time_axis, &
    &   cf_find_time_record, &
    &   cf_calendar_matches_lleapyr
    use netcdf
    implicit none

    type NCConfig
        integer :: &
        &   ncid, varid, ndims, time_dimid, time_varid, time_len
        integer, allocatable :: &
        &   shape(:)
        character(len=64) :: &
        &   time_name = ''
        type(cf_time_axis) :: &
        &   time_axis
    end type NCConfig



    interface read_nc
        module procedure :: read_nc_r4_2d
        module procedure :: read_nc_r4_3d
    end interface read_nc

contains

subroutine handle_error(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr .and. status /= -43) then
        ! status = -43: no attribute
        write(LOGNAM, *) status, trim(nf90_strerror(status))
        stop
    endif
end subroutine handle_error


integer function dimid2varid(ncid, dimid)
    integer, intent(in) :: &
    &   ncid, dimid
    character(len=64) :: &
    &   name
    call handle_error(nf90_inquire_dimension(ncid, dimid, name=name))
    call handle_error( &
    &   nf90_inq_varid(ncid, name, dimid2varid))
end function dimid2varid

! ===================================================================================================
type(NCConfig) function init_ncconfig(path, varname)
    character(len=*), intent(in) :: &
    &   path, varname
    integer :: &
    &   idim
    integer(kind=JPIM) :: &
    &   ierr
    character(len=256) :: &
    &   message
    logical :: &
    &   calendar_ok
    integer, allocatable :: &
    &   dimids(:)

    call handle_error( &
    &   nf90_open(path, nf90_nowrite, init_ncconfig%ncid))
    write(LOGNAM, '(a,i0)') '    init_ncconfig, open file: ', init_ncconfig%ncid
    call handle_error( &
    &   nf90_inq_varid(init_ncconfig%ncid, trim(varname), init_ncconfig%varid))

    call handle_error(nf90_inquire_variable( &
    &   init_ncconfig%ncid, init_ncconfig%varid, ndims=init_ncconfig%ndims))
    allocate(dimids(init_ncconfig%ndims), source=0)
    call handle_error(nf90_inquire_variable( &
    &   init_ncconfig%ncid, init_ncconfig%varid, dimids=dimids))
    allocate(init_ncconfig%shape(init_ncconfig%ndims), source=0)
    do idim = 1, init_ncconfig%ndims  ! last dim is time
        call handle_error( &
        &   nf90_inquire_dimension(init_ncconfig%ncid, dimids(idim), len=init_ncconfig%shape(idim)))
    enddo
    init_ncconfig%time_dimid = dimids(init_ncconfig%ndims)
    call handle_error(nf90_inquire_dimension( &
    &   init_ncconfig%ncid, init_ncconfig%time_dimid, &
    &   name=init_ncconfig%time_name, len=init_ncconfig%time_len))
    call handle_error(nf90_inq_varid( &
    &   init_ncconfig%ncid, trim(init_ncconfig%time_name), init_ncconfig%time_varid))
    call cf_read_time_axis( &
    &   init_ncconfig%ncid, trim(init_ncconfig%time_name), &
    &   init_ncconfig%time_axis, ierr, message)
    if (ierr /= 0_JPIM) then
        write(LOGNAM, '(2a)') '[nc_mod/init_ncconfig ERROR] ', trim(message)
        stop 9
    endif
    if (init_ncconfig%time_axis%ntime /= init_ncconfig%time_len) then
        write(LOGNAM, '(a,i0,a,i0)') &
        &   '[nc_mod/init_ncconfig ERROR] time coordinate length=', &
        &   init_ncconfig%time_axis%ntime, ', data time dimension length=', init_ncconfig%time_len
        stop 9
    endif
    calendar_ok = cf_calendar_matches_lleapyr( &
    &   init_ncconfig%time_axis%calendar, LLEAPYR, ierr, message)
    if (ierr /= 0_JPIM .or. .not. calendar_ok) then
        write(LOGNAM, '(3a,l1)') '[nc_mod/init_ncconfig ERROR] calendar=', &
        &   trim(init_ncconfig%time_axis%calendar), ', LLEAPYR=', LLEAPYR
        if (ierr /= 0_JPIM) write(LOGNAM, '(a)') trim(message)
        stop 9
    endif
end function init_ncconfig

! ===================================================================================================
subroutine get_nc_domain( &
&   ncconf, &
&   left, right, top, bottom)
    type(NCConfig), intent(in) :: &
    &   ncconf
    real(kind=JPRB), intent(out) :: &
    &   left, right, top, bottom
    real(kind=JPRB), allocatable :: &
    &   var(:)
    real(kind=JPRB) :: &
    &   res
    integer :: &
    &   status
    integer, allocatable :: &
    &   dimids(:)

    allocate(dimids(ncconf%ndims), source=0)
    status = nf90_inquire_variable( &
    &   ncconf%ncid, ncconf%varid, dimids=dimids)

    allocate(var(ncconf%shape(1)), source=-9999._JPRB)
    call handle_error( &
    &   nf90_get_var(ncconf%ncid, dimid2varid(ncconf%ncid, dimids(1)), var(:)))
    res = abs(var(1) - var(2))
    if (var(1) < var(2)) then
        left = real(nint(var(1) - res * 0.5_JPRB), kind=JPRB)
        right = real(nint(var(ncconf%shape(1)) + res * 0.5_JPRB), kind=JPRB)
    else
        left = real(nint(var(1) + res * 0.5_JPRB), kind=JPRB)
        right = real(nint(var(ncconf%shape(1)) - res * 0.5_JPRB), kind=JPRB)
    endif
    deallocate(var)

    allocate(var(ncconf%shape(2)), source=-9999._JPRB)
    call handle_error( &
    &   nf90_get_var(ncconf%ncid, dimid2varid(ncconf%ncid, dimids(2)), var(:)))
    if (var(1) > var(2)) then
        top = real(nint(var(1) + res * 0.5_JPRB), kind=JPRB)
        bottom = real(nint(var(ncconf%shape(2)) - res * 0.5_JPRB), kind=JPRB)
    else
        top = real(nint(var(1) - res * 0.5_JPRB), kind=JPRB)
        bottom = real(nint(var(ncconf%shape(2)) + res * 0.5_JPRB), kind=JPRB)
    endif
    deallocate(var)
end subroutine get_nc_domain


integer function get_nc_dt(ncconf)
    type(NCConfig), intent(in) :: &
    &   ncconf
    get_nc_dt = int(ncconf%time_axis%dt_minutes * 60_JPIB)
end function get_nc_dt


integer(kind=JPIM) function get_nc_start_record(ncconf, start_dt) result(record)
    type(NCConfig), intent(in) :: &
    &   ncconf
    type(DateTime), intent(in) :: &
    &   start_dt
    integer(kind=JPIM) :: &
    &   ierr
    character(len=256) :: &
    &   message

    call cf_find_time_record( &
    &   ncconf%time_axis, start_dt%yyyymmdd, start_dt%hour, 0_JPIM, &
    &   record, ierr, message, .TRUE.)
    if (ierr /= 0_JPIM) then
        write(LOGNAM, '(2a)') '[nc_mod/get_nc_start_record ERROR] ', trim(message)
        stop 9
    endif
end function get_nc_start_record


!subroutine get_nc_scale_offset(unit, var_id, scale, offset)
!    integer, intent(in)  :: unit, var_id
!    real(kind=JPRD), intent(out) :: scale, offset
!    call handle_error(nf90_get_att(unit, var_id, 'scale_factor', scale ))
!    call handle_error(nf90_get_att(unit, var_id, 'add_offset'  , offset))
!    if (scale == 0.d0) scale = 1.d0
!end subroutine get_nc_scale_offset

! ===================================================================================================
subroutine check_get_var_error( &
&   status, &
&   file_is_end)
    integer, intent(in) :: &
    &   status
    logical, intent(out) :: &
    &   file_is_end
    if (status == nf90_noerr) then
        file_is_end = .FALSE.
    elseif (status == nf90_einvalcoords) then
        file_is_end = .TRUE.
    else
        call handle_error(status)
    endif
end subroutine check_get_var_error


subroutine read_nc_r4_2d( &
&   arr, file_is_end, &
&   ncconf, recnum)
    real(kind=JPRM), intent(out) :: &
    &   arr(:,:)
    logical, intent(out) :: &
    &   file_is_end
    type(NCConfig), intent(in) :: &
    &   ncconf
    integer, intent(in) :: &
    &   recnum
    integer :: &
    &   status
    status = nf90_get_var( &
    &   ncconf%ncid, ncconf%varid, arr, &
    &   start=[1,1,recnum], count=[ncconf%shape(1),ncconf%shape(2),1])
    call check_get_var_error( &
    &   status, &
    &   file_is_end)
end subroutine read_nc_r4_2d


subroutine read_nc_r4_3d( &
&   arr, file_is_end, &
&   ncconf, recnum)
    real(kind=JPRM), intent(out) :: &
    &   arr(:,:,:)
    logical, intent(out) :: &
    &   file_is_end
    type(NCConfig), intent(in) :: &
    &   ncconf
    integer, intent(in) :: &
    &   recnum
    integer :: &
    &   status
    status = nf90_get_var( &
    &   ncconf%ncid, ncconf%varid, arr, &
    &   start=[1,1,1,recnum], count=[ncconf%shape(1),ncconf%shape(2),ncconf%shape(3),1])
    call check_get_var_error( &
    &   status, &
    &   file_is_end)
end subroutine read_nc_r4_3d
#endif
end module nc_mod
