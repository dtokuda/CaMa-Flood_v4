program test_heatlink_config
    use heatlink_config_mod, only: &
    &   LICE, NNEWTON_MAX_ICE, LHEAT_DIAG, CHEAT_LOG, init_heatlink_config
    implicit none
    character(len = 512) :: path ! [-] Optional fixture path for expected-error tests.

    call get_command_argument(1, path)
    if (len_trim(path) > 0) then
        call init_heatlink_config(trim(path), 6, .false., .false.)
        stop
    endif

    call init_heatlink_config( &
    &   'test/heatlink_config_defaults.nml', 6, .false., .false.)
    if (LHEAT_DIAG) error stop 'LHEAT_DIAG default is not false'
    if (CHEAT_LOG /= 'log_HEAT-LINK.txt') error stop 'CHEAT_LOG default is incorrect'
    if (LICE) error stop 'LICE default is not false'
    if (NNEWTON_MAX_ICE /= 4) error stop 'NNEWTON_MAX_ICE default is not 4'

    call init_heatlink_config( &
    &   'test/heatlink_config_enabled.nml', 6, .false., .false.)
    if (.not. LICE) error stop 'LICE was not read from NHEATLINK'
    if (NNEWTON_MAX_ICE /= 7) error stop 'NNEWTON_MAX_ICE was not read from NHEATLINK'

    if (.not. LHEAT_DIAG) error stop 'LHEAT_DIAG was not read'
    if (CHEAT_LOG /= 'custom_heat.log') error stop 'CHEAT_LOG was not read'
    call init_heatlink_config('test/heatlink_config_defaults.nml', 6, .false., .false.)
    if (LHEAT_DIAG .or. CHEAT_LOG /= 'log_HEAT-LINK.txt') error stop 'monitor settings were not reset'

    write(*, '(a)') 'test_heatlink_config: PASS'
end program test_heatlink_config
