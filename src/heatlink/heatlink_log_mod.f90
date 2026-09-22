module heatlink_log_mod
    use, intrinsic :: iso_fortran_env, only: output_unit
    use heatlink_config_mod, only: CHEAT_LOG, LHEAT_DIAG
    implicit none
    private
    public :: HEAT_LOG_UNIT, init_heatlink_log, fin_heatlink_log, write_heatlink_time

    integer, protected, save :: HEAT_LOG_UNIT = output_unit ! [-] Heatlink log unit; stdout for standalone kernels.
    logical, save :: log_open = .false. ! [-] Whether this module owns an open log file.
contains

subroutine init_heatlink_log(cama_log_unit)
    integer, intent(in) :: cama_log_unit ! [-] Existing CaMa log unit for file-opening errors.
    integer :: ios ! [-] I/O status from opening the heat log.
    logical :: already_open ! [-] Whether the requested file is already in use (including the CaMa log).
    character(len = 512) :: message ! [-] Runtime I/O error description.

    if (log_open) call fin_heatlink_log()
    inquire(file = trim(CHEAT_LOG), opened = already_open, iostat = ios)
    if (ios /= 0 .or. already_open .or. len_trim(CHEAT_LOG) == 0) then
        write(cama_log_unit, '(a,1x,a)') 'ERROR: CHEAT_LOG is empty, inaccessible or already open:', trim(CHEAT_LOG)
        error stop 1
    endif
    open(newunit = HEAT_LOG_UNIT, file = trim(CHEAT_LOG), status = 'replace', action = 'write', &
    &   iostat = ios, iomsg = message)
    if (ios /= 0) then
        write(cama_log_unit, '(a,1x,a,2a)') 'ERROR: cannot open CHEAT_LOG:', trim(CHEAT_LOG), ': ', trim(message)
        error stop 1
    endif
    log_open = .true.
    write(cama_log_unit, '(a,1x,a)') 'HEAT-LINK log:', trim(CHEAT_LOG)
    write(HEAT_LOG_UNIT, '(a,l1)') 'HEAT-LINK detailed monitoring: LHEAT_DIAG = ', LHEAT_DIAG
    write(HEAT_LOG_UNIT, '(a)') 'HEAT_TIME: model calendar; BEGIN/END bracket one outer update. Internal HEAT_STEP times are elapsed seconds from run start.'
    if (.not. LHEAT_DIAG) write(HEAT_LOG_UNIT, '(a)') &
    &   'Detailed heat budgets are disabled; temperature checks and physical safeguards remain active.'
end subroutine init_heatlink_log

subroutine write_heatlink_time(stage, step, date, hhmm)
    character(len = *), intent(in) :: stage ! [-] INIT, BEGIN, LOCAL_END or END marker.
    integer, intent(in) :: step ! [-] CaMa outer time-step counter at this marker.
    integer, intent(in) :: date ! [YYYYMMDD] Model calendar date from YOS_CMF_TIME.
    integer, intent(in) :: hhmm ! [HHMM] Model time from YOS_CMF_TIME; not wall-clock time.
    write(HEAT_LOG_UNIT, '(a,1x,a,1x,i0,1x,i8.8,1x,i4.4)') 'HEAT_TIME', stage, step, date, hhmm
end subroutine write_heatlink_time

subroutine fin_heatlink_log()
    if (log_open) close(HEAT_LOG_UNIT)
    log_open = .false.
    HEAT_LOG_UNIT = output_unit
end subroutine fin_heatlink_log
end module heatlink_log_mod
