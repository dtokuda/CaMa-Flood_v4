module heatlink_diagnostics_mod
#ifdef heatlink
    use PARKIND1, only: JPIM, JPRB, JPRD
    use heatlink_log_mod, only: HEAT_LOG_UNIT
    use YOS_CMF_MAP, only: NSEQALL
    use YOS_CMF_PROG, only: P2RIVSTO, P2FLDSTO
    use heatlink_config_mod, only: LICE, LHEAT_DIAG
    use const_mod, only: STO_IGNORE
    use phys_const_mod, only: CW, RW, RI, HFUS, TMELT
    use heat_residual_mod, only: HeatResidualStats, record_heat_residual, write_heat_residual, &
    &   HeatConservationStats, record_heat_exchange, write_heat_conservation
    use heat_step_monitor_mod, only: HeatStepState, HeatStepLedger, HeatStepStats, &
    &   capture_heat_step, measure_heat_step, add_step_heat, step_sum, monitor_heat_step, write_heat_step_extrema
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: init_heatlink_diagnostics, fin_heatlink_diagnostics, check_heatlink_temperature
    public :: begin_advection_diagnostics, record_advection_diagnostics, finish_advection_diagnostics
    public :: begin_local_diagnostics, finish_local_diagnostics

    ! Coupling diagnostics own their bookkeeping, never the river's physical state.
    type(HeatResidualStats), save :: residual_stats(5) ! [mixed] Cause-specific energy [J], counts and cell indices [-].
    type(HeatResidualStats), save :: closure_stats(2) ! [mixed] Advection/local closure energy [J] and counts [-].
    type(HeatConservationStats), save :: conservation_stats ! [mixed] Run-wide heat ledger [J] and process counts [-].
    type(HeatStepState), save :: process_start ! [mixed] Start-of-process volumes [m3] and temperature offsets [K].
    type(HeatStepState), save :: hour_start ! [mixed] Start-of-outer-update volumes [m3] and temperature offsets [K].
    type(HeatStepLedger), save :: process_ledger ! [J] Heat exchanges during one process interval.
    type(HeatStepLedger), save :: hour_ledger ! [J] Heat exchanges during one outer update.
    type(HeatStepStats), save :: step_stats(3) ! [mixed] Advection/local/hour extrema; field units are defined by HeatStepStats.
    logical, save :: hour_monitor_active = .false. ! [-] Whether an outer-update snapshot is active.
    real(kind = JPRD), save :: monitor_seconds = 0.0_JPRD ! [s] Elapsed time accumulated from internal steps.
    real(kind = JPRD), save :: hour_start_seconds = 0.0_JPRD ! [s] Elapsed time at the current outer-update start.
    real(kind = JPRD), save :: advection_dt_seconds = 0.0_JPRD ! [s] Duration of the current advection interval.
    character(len = 24), parameter :: residual_reasons(5) = [character(len = 24) :: & ! [-] Labels for the five unapplied-heat causes.
    &   'advection_dry', 'advection_reconstruction', 'dry_local', 'ice_model', 'no_ice_floor']
contains

subroutine init_heatlink_diagnostics()
    if (.not. LHEAT_DIAG) return
    residual_stats = HeatResidualStats()
    closure_stats = HeatResidualStats()
    conservation_stats = HeatConservationStats()
    step_stats = HeatStepStats()
    process_ledger = HeatStepLedger()
    hour_ledger = HeatStepLedger()
    hour_monitor_active = .false.
    monitor_seconds = 0.0_JPRD
    hour_start_seconds = 0.0_JPRD
    write(HEAT_LOG_UNIT,'(a)') 'HEAT_STEP columns: stage step elapsed_seconds dt_seconds delta[J] net[J] exchange_abs[J] unapplied[J] unapplied_abs[J] raw[J] adjusted[J] storage_scale[J] raw_exchange_ratio adjusted_exchange_ratio unapplied_exchange_ratio adjusted_storage_ratio naive_delta[J] naive_adjusted[J] exchange_ratio_valid storage_ratio_valid'
    write(HEAT_LOG_UNIT,'(a)') 'HEAT_STEP uses cellwise state increments and interval-only compensated sums. Invalid ratios are placeholders; do not interpret zero as a valid ratio.'
    write(HEAT_LOG_UNIT,'(a)') 'HEAT_STEP_EXTREME columns: stage metric min_or_max step followed by the same 16 real fields as HEAT_STEP. Extrema cover this run; errors are not accumulated.'
    write(HEAT_LOG_UNIT, '(a)') 'HEAT_RESIDUAL columns: reason signed[J] positive[J] negative[J] absolute[J] max[J] events cells max_cell throughput[J]'
    write(HEAT_LOG_UNIT, '(a)') 'HEAT_CONSERVATION columns: initial[J] current[J] water_net[J] ice_net[J] local_net[J] water_abs[J] ice_abs[J] local_abs[J] unapplied[J] raw[J] adjusted[J] raw_fraction adjusted_fraction water_steps ice_steps local_steps'
    write(HEAT_LOG_UNIT, '(a)') 'HEAT_CONSERVATION raw = current - initial - net_input; adjusted = raw + signed_unapplied. Ice storage is melting-point latent energy; ice skin is massless.'
    write(HEAT_LOG_UNIT, '(a)') 'HEAT_RESIDUAL totals start at this run; diagnostics are never reinjected. Internal transport is counted at both ends.'
    advection_dt_seconds = 0.0_JPRD
end subroutine init_heatlink_diagnostics

! Independent state-versus-external-input audit, using canonical hydraulic storage.
! No diagnostic value is fed back into water, ice, temperature, or flow.
function represented_domain_energy_j(wattmp, icevol, icevol_excess) result(energy_j)
    real(kind = JPRB), intent(in) :: wattmp(:) ! [K] Current liquid-water temperature.
    real(kind = JPRB), allocatable, intent(in) :: icevol(:) ! [m3] Mobile ice; unallocated when LICE is false.
    real(kind = JPRB), allocatable, intent(in) :: icevol_excess(:) ! [m3] Immobile ice; unallocated when LICE is false.
    real(kind = JPRD) :: energy_j ! [J] Domain sensible heat plus melting-point ice latent energy.
    energy_j = real(CW, JPRD) * real(RW, JPRD) * sum( &
    &   (P2RIVSTO(:NSEQALL,1) + P2FLDSTO(:NSEQALL,1)) * real(wattmp(:NSEQALL) - TMELT, JPRD))
    if (LICE) energy_j = energy_j - real(RI, JPRD) * real(HFUS, JPRD) * &
    &   (sum(real(icevol(:NSEQALL), JPRD)) + sum(real(icevol_excess(:NSEQALL), JPRD)))
end function

! Use the same canonical storage and energy reference as the annual ledger.
subroutine capture_monitor_state(state, wattmp, icevol, icevol_excess)
    real(kind = JPRB), intent(in) :: wattmp(:) ! [K] Current liquid-water temperature.
    real(kind = JPRB), allocatable, intent(in) :: icevol(:) ! [m3] Mobile ice; unallocated when LICE is false.
    real(kind = JPRB), allocatable, intent(in) :: icevol_excess(:) ! [m3] Immobile ice; unallocated when LICE is false.
    type(HeatStepState), intent(inout) :: state ! [m3,K] Read-only physical-state snapshot to capture.
    if (LICE) then
        call capture_heat_step(state,P2RIVSTO(:NSEQALL,1)+P2FLDSTO(:NSEQALL,1), &
        &   real(wattmp(:NSEQALL)-TMELT,JPRD),real(icevol(:NSEQALL),JPRD),real(icevol_excess(:NSEQALL),JPRD))
    else
        call capture_heat_step(state,P2RIVSTO(:NSEQALL,1)+P2FLDSTO(:NSEQALL,1),real(wattmp(:NSEQALL)-TMELT,JPRD))
    endif
end subroutine

subroutine finish_monitor_state(state,ledger,index,stage,dt_seconds, wattmp, icevol, icevol_excess)
    real(kind = JPRB), intent(in) :: wattmp(:) ! [K] Current liquid-water temperature.
    real(kind = JPRB), allocatable, intent(in) :: icevol(:) ! [m3] Mobile ice; unallocated when LICE is false.
    real(kind = JPRB), allocatable, intent(in) :: icevol_excess(:) ! [m3] Immobile ice; unallocated when LICE is false.
    type(HeatStepState), intent(in) :: state ! [m3,K] Physical-state snapshot at interval start.
    type(HeatStepLedger), intent(in) :: ledger ! [J] Expected and unapplied heat during this interval.
    integer, intent(in) :: index ! [-] Stage index: advection, local, or outer update.
    character(len = *), intent(in) :: stage ! [-] Stage label written to the diagnostic log.
    real(kind = JPRD), intent(in) :: dt_seconds ! [s] Duration of the monitored interval.
    real(kind = JPRD) :: delta_j ! [J] Domain energy change from cellwise state increments.
    real(kind = JPRD) :: storage_j ! [J] Non-cancelling represented storage scale.
    real(kind = JPRD) :: naive_delta_j ! [J] Energy change from subtracting domain totals, for comparison.
    if (LICE) then
        call measure_heat_step(state,P2RIVSTO(:NSEQALL,1)+P2FLDSTO(:NSEQALL,1), &
        &   real(wattmp(:NSEQALL)-TMELT,JPRD),real(CW,JPRD)*real(RW,JPRD),real(RI,JPRD)*real(HFUS,JPRD), &
        &   delta_j,storage_j,naive_delta_j,real(icevol(:NSEQALL),JPRD),real(icevol_excess(:NSEQALL),JPRD))
    else
        call measure_heat_step(state,P2RIVSTO(:NSEQALL,1)+P2FLDSTO(:NSEQALL,1), &
        &   real(wattmp(:NSEQALL)-TMELT,JPRD),real(CW,JPRD)*real(RW,JPRD),real(RI,JPRD)*real(HFUS,JPRD), &
        &   delta_j,storage_j,naive_delta_j)
    endif
    call monitor_heat_step(step_stats(index),HEAT_LOG_UNIT,stage,monitor_seconds,dt_seconds, &
    &   delta_j,storage_j,naive_delta_j,ledger)
end subroutine

subroutine add_process_to_hour()
    real(kind = JPRD) :: q(4) ! [J] Compensated net, absolute exchange, signed and absolute unapplied heat.
    q = process_ledger%value+process_ledger%correction
    call add_step_heat(hour_ledger,q(1),q(2),q(3),q(4))
end subroutine

subroutine begin_heat_conservation(wattmp, icevol, icevol_excess)
    real(kind = JPRB), intent(in) :: wattmp(:) ! [K] Current liquid-water temperature.
    real(kind = JPRB), allocatable, intent(in) :: icevol(:) ! [m3] Mobile ice; unallocated when LICE is false.
    real(kind = JPRB), allocatable, intent(in) :: icevol_excess(:) ! [m3] Immobile ice; unallocated when LICE is false.
    if (conservation_stats%initialized) return
    conservation_stats%initial_j = represented_domain_energy_j(wattmp, icevol, icevol_excess)
    conservation_stats%initialized = .true.
end subroutine

subroutine begin_advection_diagnostics(wattmp, icevol, icevol_excess)
    real(kind = JPRB), intent(in) :: wattmp(:) ! [K] Current liquid-water temperature.
    real(kind = JPRB), allocatable, intent(in) :: icevol(:) ! [m3] Mobile ice; unallocated when LICE is false.
    real(kind = JPRB), allocatable, intent(in) :: icevol_excess(:) ! [m3] Immobile ice; unallocated when LICE is false.
    if (.not. LHEAT_DIAG) return
    call begin_heat_conservation(wattmp, icevol, icevol_excess)
    call capture_monitor_state(process_start, wattmp, icevol, icevol_excess)
    process_ledger = HeatStepLedger()
    if (.not. hour_monitor_active) then
        call capture_monitor_state(hour_start, wattmp, icevol, icevol_excess)
        hour_start_seconds = monitor_seconds
        hour_ledger = HeatStepLedger()
        hour_monitor_active = .true.
    endif
end subroutine begin_advection_diagnostics

subroutine record_advection_diagnostics(dt_seconds, boundary_heat_j, boundary_absolute_j, ice_export_m3, &
&   advection_unapplied_sensible_heat_j, advection_throughput_j, &
&   advection_domain_heat_budget_error_j, advection_domain_combined_energy_budget_error_j)
    real(kind = JPRB), intent(in) :: dt_seconds ! [s] Duration of this hydraulic transport step.
    real(kind = JPRD), intent(in) :: boundary_heat_j ! [J] Net sensible heat entering across external boundaries.
    real(kind = JPRD), intent(in) :: boundary_absolute_j ! [J] Absolute sensible-heat exchange across external boundaries.
    real(kind = JPRD), intent(in) :: ice_export_m3 ! [m3] Exported ice volume; zero when ice is disabled.
    real(kind = JPRD), intent(in) :: advection_unapplied_sensible_heat_j(:) ! [J] Signed cellwise unapplied transport heat.
    real(kind = JPRD), intent(in) :: advection_throughput_j(:) ! [J] Cellwise absolute heat-transport scale.
    real(kind = JPRD), intent(in) :: advection_domain_heat_budget_error_j ! [J] Domain water-transport closure error.
    real(kind = JPRD), intent(in) :: advection_domain_combined_energy_budget_error_j ! [J] Domain water-plus-ice closure error.
    real(kind = JPRD) :: volumetric_ice_latent_energy_j_m3 ! [J m-3] Magnitude of melting-point ice latent energy.

    if (.not. LHEAT_DIAG) return
    advection_dt_seconds = real(dt_seconds,JPRD)
    volumetric_ice_latent_energy_j_m3 = real(RI,JPRD) * real(HFUS,JPRD)
    call add_step_heat(process_ledger,boundary_heat_j,boundary_absolute_j, &
    &   step_sum(advection_unapplied_sensible_heat_j(:NSEQALL)), &
    &   step_sum(abs(advection_unapplied_sensible_heat_j(:NSEQALL))))
    if (LICE) call add_step_heat(process_ledger,volumetric_ice_latent_energy_j_m3*ice_export_m3, &
    &   volumetric_ice_latent_energy_j_m3*ice_export_m3,0.0_JPRD,0.0_JPRD)
    call record_heat_exchange(conservation_stats, 1, boundary_heat_j, boundary_absolute_j)
    if (LICE) then
        ! Exporting negative latent energy adds energy to the remaining domain.
        call record_heat_exchange(conservation_stats, 2, &
        &   volumetric_ice_latent_energy_j_m3 * ice_export_m3, volumetric_ice_latent_energy_j_m3 * ice_export_m3)
        call record_heat_residual(closure_stats(1), [advection_domain_combined_energy_budget_error_j], &
        &   [boundary_absolute_j + volumetric_ice_latent_energy_j_m3 * ice_export_m3])
    else
        call record_heat_residual(closure_stats(1), [advection_domain_heat_budget_error_j], [boundary_absolute_j])
    endif

    call record_heat_residual(residual_stats(1), advection_unapplied_sensible_heat_j(:NSEQALL), &
    &   advection_throughput_j(:NSEQALL), P2RIVSTO(:NSEQALL,1) + P2FLDSTO(:NSEQALL,1) <= STO_IGNORE)
    call record_heat_residual(residual_stats(2), advection_unapplied_sensible_heat_j(:NSEQALL), &
    &   advection_throughput_j(:NSEQALL), P2RIVSTO(:NSEQALL,1) + P2FLDSTO(:NSEQALL,1) > STO_IGNORE)

end subroutine record_advection_diagnostics

subroutine finish_advection_diagnostics(wattmp, icevol, icevol_excess)
    real(kind = JPRB), intent(in) :: wattmp(:) ! [K] Current liquid-water temperature.
    real(kind = JPRB), allocatable, intent(in) :: icevol(:) ! [m3] Mobile ice; unallocated when LICE is false.
    real(kind = JPRB), allocatable, intent(in) :: icevol_excess(:) ! [m3] Immobile ice; unallocated when LICE is false.
    if (.not. LHEAT_DIAG) return
    monitor_seconds = monitor_seconds + advection_dt_seconds
    call finish_monitor_state(process_start,process_ledger,1,'advection',advection_dt_seconds, wattmp, icevol, icevol_excess)
    call add_process_to_hour()
end subroutine finish_advection_diagnostics

subroutine begin_local_diagnostics(wattmp, icevol, icevol_excess)
    real(kind = JPRB), intent(in) :: wattmp(:) ! [K] Current liquid-water temperature.
    real(kind = JPRB), allocatable, intent(in) :: icevol(:) ! [m3] Mobile ice; unallocated when LICE is false.
    real(kind = JPRB), allocatable, intent(in) :: icevol_excess(:) ! [m3] Immobile ice; unallocated when LICE is false.
    if (.not. LHEAT_DIAG) return
    call begin_heat_conservation(wattmp, icevol, icevol_excess)
    call capture_monitor_state(process_start, wattmp, icevol, icevol_excess)
    process_ledger = HeatStepLedger()
end subroutine begin_local_diagnostics

subroutine finish_local_diagnostics(dt, wattmp, icevol, icevol_excess, &
&   local_added_energy_j, local_dry_energy_j, local_throughput_j, floor_energy_j, &
&   phase_unapplied_energy, phase_energy_budget_error)
    real(kind = JPRB), intent(in) :: wattmp(:) ! [K] Current liquid-water temperature.
    real(kind = JPRB), allocatable, intent(in) :: icevol(:) ! [m3] Mobile ice; unallocated when LICE is false.
    real(kind = JPRB), allocatable, intent(in) :: icevol_excess(:) ! [m3] Immobile ice; unallocated when LICE is false.
    real(kind = JPRB), intent(in) :: dt ! [s] Duration of the local heat update.
    real(kind = JPRB), intent(in) :: local_added_energy_j(:) ! [J] Expected net local heat input per cell.
    real(kind = JPRB), intent(in) :: local_dry_energy_j(:) ! [J] Signed local heat skipped by dry handling.
    real(kind = JPRB), intent(in) :: local_throughput_j(:) ! [J] Absolute local heat-input scale per cell.
    real(kind = JPRB), intent(in) :: floor_energy_j(:) ! [J] Signed heat omitted by the no-ice melting-point floor.
    real(kind = JPRB), allocatable, intent(in) :: phase_unapplied_energy(:) ! [J] Signed heat unapplied by phase handling; ice only.
    real(kind = JPRB), allocatable, intent(in) :: phase_energy_budget_error(:) ! [J] Local phase-change closure error; ice only.
    integer(kind = JPIM) :: reason ! [-] Index of an unapplied-heat cause.

    if (.not. LHEAT_DIAG) return
    if (LICE) then
        call record_heat_residual(residual_stats(4), &
        &   real(phase_unapplied_energy(:NSEQALL) - local_dry_energy_j(:NSEQALL), JPRD), &
        &   real(local_throughput_j(:NSEQALL), JPRD))
        call record_heat_residual(closure_stats(2), -real(phase_energy_budget_error(:NSEQALL), JPRD), &
        &   abs(real(local_added_energy_j(:NSEQALL), JPRD)))
    else
        call record_heat_residual(residual_stats(5), real(floor_energy_j(:NSEQALL), JPRD), &
        &   real(local_throughput_j(:NSEQALL), JPRD))
    endif

    if (LICE) then
        call add_step_heat(process_ledger,step_sum(real(local_added_energy_j(:NSEQALL),JPRD)), &
        &   step_sum(abs(real(local_added_energy_j(:NSEQALL),JPRD))), &
        &   step_sum(real(phase_unapplied_energy(:NSEQALL),JPRD)), &
        &   step_sum(abs(real(local_dry_energy_j(:NSEQALL),JPRD))) + &
        &   step_sum(abs(real(phase_unapplied_energy(:NSEQALL)-local_dry_energy_j(:NSEQALL),JPRD))))
    else
        call add_step_heat(process_ledger,step_sum(real(local_added_energy_j(:NSEQALL),JPRD)), &
        &   step_sum(abs(real(local_added_energy_j(:NSEQALL),JPRD))), &
        &   step_sum(real(local_dry_energy_j(:NSEQALL),JPRD)) + step_sum(real(floor_energy_j(:NSEQALL),JPRD)), &
        &   step_sum(abs(real(local_dry_energy_j(:NSEQALL),JPRD))) + step_sum(abs(real(floor_energy_j(:NSEQALL),JPRD))))
    endif
    call finish_monitor_state(process_start,process_ledger,2,'local',real(dt,JPRD), wattmp, icevol, icevol_excess)
    call add_process_to_hour()
    if (hour_monitor_active) then
        call finish_monitor_state(hour_start,hour_ledger,3,'hour',monitor_seconds-hour_start_seconds, wattmp, icevol, icevol_excess)
        hour_monitor_active = .false.
    endif
    call record_heat_residual(residual_stats(3), real(local_dry_energy_j(:NSEQALL), JPRD), &
    &   real(local_throughput_j(:NSEQALL), JPRD))
    do reason = 1, size(residual_stats)
        call write_heat_residual(HEAT_LOG_UNIT, residual_reasons(reason), residual_stats(reason))
    enddo
    call record_heat_exchange(conservation_stats, 3, sum(real(local_added_energy_j(:NSEQALL), JPRD)), &
    &   sum(abs(real(local_added_energy_j(:NSEQALL), JPRD))))
    call write_heat_conservation(HEAT_LOG_UNIT, conservation_stats, represented_domain_energy_j(wattmp, icevol, icevol_excess), &
    &   sum(residual_stats%positive_j) + sum(residual_stats%negative_j))
    call write_heat_residual(HEAT_LOG_UNIT, 'advection_domain', closure_stats(1), 'HEAT_CLOSURE')
    if (LICE) call write_heat_residual(HEAT_LOG_UNIT, 'local_phase', closure_stats(2), 'HEAT_CLOSURE')

end subroutine finish_local_diagnostics

subroutine check_heatlink_temperature(wattmp, watsto)
    real(kind = JPRB), intent(in) :: wattmp(:) ! [K] End-of-update liquid-water temperature.
    real(kind = JPRB), intent(in) :: watsto(:) ! [m3] End-of-update liquid-water volume.
    integer(kind = JPIM) :: max_cell(1) ! [-] One-based cell index of the maximum temperature.
    logical :: wet(NSEQALL) ! [-] Cells with liquid volume above STO_IGNORE.

    if (.not. all(ieee_is_finite(wattmp(:NSEQALL)))) then
        write(HEAT_LOG_UNIT, '(a)') 'ERROR: non-finite river water temperature.'
        flush(HEAT_LOG_UNIT)
        error stop 'Non-finite river water temperature; see CHEAT_LOG.'
    endif
    wet = watsto(:NSEQALL) > real(STO_IGNORE, JPRB)
    max_cell = maxloc(wattmp(:NSEQALL), mask = wet)
    if (any(wet)) then
        write(HEAT_LOG_UNIT, '(a,2(1x,i0),3(1x,es24.16))') 'THERMAL_WET', count(wet), max_cell(1), &
        &   minval(wattmp(:NSEQALL), mask = wet), maxval(wattmp(:NSEQALL), mask = wet), watsto(max_cell(1))
    endif
    max_cell = maxloc(wattmp(:NSEQALL), mask = .not. wet)
    if (any(.not. wet)) then
        write(HEAT_LOG_UNIT, '(a,2(1x,i0),3(1x,es24.16))') 'THERMAL_DRY', count(.not. wet), max_cell(1), &
        &   minval(wattmp(:NSEQALL), mask = .not. wet), maxval(wattmp(:NSEQALL), mask = .not. wet), watsto(max_cell(1))
    endif
    if (maxval(wattmp(:NSEQALL)) > 350.0_JPRB) &
    &   write(HEAT_LOG_UNIT, '(a)') 'WARNING: river water temperature exceeds 350 K; inspect THERMAL_WET/THERMAL_DRY.'
end subroutine check_heatlink_temperature

subroutine fin_heatlink_diagnostics()
    if (.not. LHEAT_DIAG) return
    call write_heat_step_extrema(HEAT_LOG_UNIT,'advection',step_stats(1))
    call write_heat_step_extrema(HEAT_LOG_UNIT,'local',step_stats(2))
    call write_heat_step_extrema(HEAT_LOG_UNIT,'hour',step_stats(3))
    process_start = HeatStepState()
    hour_start = HeatStepState()
end subroutine fin_heatlink_diagnostics
#endif
end module heatlink_diagnostics_mod
