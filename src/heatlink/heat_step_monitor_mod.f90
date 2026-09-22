module heat_step_monitor_mod
    use PARKIND1, only: JPIM, JPIB, JPRD
    implicit none
    private
    public :: HeatStepState, HeatStepLedger, HeatStepStats, capture_heat_step, measure_heat_step
    public :: add_step_heat, step_sum, monitor_heat_step, write_heat_step_extrema

    ! Independent diagnostics: snapshots and interval sums never update model state.
    type HeatStepState
        real(kind=JPRD), allocatable :: volume(:), theta(:), ice(:), excess(:)
    end type
    type HeatStepLedger
        real(kind=JPRD) :: value(4) = 0.0_JPRD, correction(4) = 0.0_JPRD
    end type
    integer, parameter :: NMETRIC = 8, NFIELD = 16
    integer, parameter :: metric_field(NMETRIC) = [8,9,6,7,11,12,13,14]
    character(len=24), parameter :: metric_name(NMETRIC) = [character(len=24) :: &
    &   'raw_j','adjusted_j','unapplied_j','absolute_unapplied_j', &
    &   'raw_exchange_ratio','adjusted_exchange_ratio','unapplied_exchange_ratio','adjusted_storage_ratio']
    type HeatStepStats
        integer(kind=JPIB) :: steps = 0_JPIB, samples(NMETRIC) = 0_JPIB
        integer(kind=JPIB) :: no_exchange_ratio = 0_JPIB, no_storage_ratio = 0_JPIB
        integer(kind=JPIB) :: min_step(NMETRIC) = 0_JPIB, max_step(NMETRIC) = 0_JPIB
        real(kind=JPRD) :: min_record(NFIELD,NMETRIC) = 0.0_JPRD, max_record(NFIELD,NMETRIC) = 0.0_JPRD
    end type
contains

! Neumaier summation also handles a small partial sum followed by a larger term.
subroutine add_compensated(total, correction, value)
    real(kind=JPRD), intent(inout) :: total, correction
    real(kind=JPRD), intent(in) :: value
    real(kind=JPRD) :: updated
    updated = total + value
    if (abs(total) >= abs(value)) then
        correction = correction + ((total - updated) + value)
    else
        correction = correction + ((value - updated) + total)
    endif
    total = updated
end subroutine

function step_sum(values) result(total)
    real(kind=JPRD), intent(in) :: values(:)
    real(kind=JPRD) :: total, correction
    integer :: i
    total = 0.0_JPRD
    correction = 0.0_JPRD
    do i = 1, size(values)
        call add_compensated(total, correction, values(i))
    enddo
    total = total + correction
end function

subroutine add_step_heat(ledger, net_j, exchange_j, unapplied_j, absolute_unapplied_j)
    type(HeatStepLedger), intent(inout) :: ledger
    real(kind=JPRD), intent(in) :: net_j, exchange_j, unapplied_j, absolute_unapplied_j
    real(kind=JPRD) :: values(4)
    integer :: i
    values = [net_j,exchange_j,unapplied_j,absolute_unapplied_j]
    do i = 1, 4
        call add_compensated(ledger%value(i),ledger%correction(i),values(i))
    enddo
end subroutine

subroutine capture_heat_step(state, volume, theta, ice, excess)
    type(HeatStepState), intent(inout) :: state
    real(kind=JPRD), intent(in) :: volume(:), theta(:)
    real(kind=JPRD), intent(in), optional :: ice(:), excess(:)
    state%volume = volume
    state%theta = theta
    if (present(ice)) state%ice = ice
    if (present(excess)) state%excess = excess
end subroutine

subroutine measure_heat_step(state, volume, theta, capacity, latent, delta_j, storage_j, naive_delta_j, ice, excess)
    type(HeatStepState), intent(in) :: state
    real(kind=JPRD), intent(in) :: volume(:), theta(:), capacity, latent
    real(kind=JPRD), intent(in), optional :: ice(:), excess(:)
    real(kind=JPRD), intent(out) :: delta_j, storage_j, naive_delta_j
    real(kind=JPRD) :: correction, scale_correction, initial_j, final_j
    integer :: i
    delta_j = 0.0_JPRD
    storage_j = 0.0_JPRD
    correction = 0.0_JPRD
    scale_correction = 0.0_JPRD
    do i = 1, size(volume)
        ! Difference of products without first constructing two large cell energies.
        call add_compensated(delta_j,correction,capacity * state%volume(i) * (theta(i)-state%theta(i)))
        call add_compensated(delta_j,correction,capacity * (volume(i)-state%volume(i)) * theta(i))
        call add_compensated(storage_j,scale_correction,capacity * &
        &   max(abs(state%volume(i)*state%theta(i)),abs(volume(i)*theta(i))))
        if (present(ice)) then
            call add_compensated(delta_j,correction,latent*(state%ice(i)-ice(i)))
            call add_compensated(storage_j,scale_correction,latent*max(abs(state%ice(i)),abs(ice(i))))
        endif
        if (present(excess)) then
            call add_compensated(delta_j,correction,latent*(state%excess(i)-excess(i)))
            call add_compensated(storage_j,scale_correction,latent*max(abs(state%excess(i)),abs(excess(i))))
        endif
    enddo
    delta_j = delta_j + correction
    storage_j = storage_j + scale_correction
    ! Retain the large-total subtraction only as a comparison, never as the main monitor.
    initial_j = capacity*sum(state%volume*state%theta)
    final_j = capacity*sum(volume*theta)
    if (present(ice)) then
        initial_j = initial_j - latent*sum(state%ice)
        final_j = final_j - latent*sum(ice)
    endif
    if (present(excess)) then
        initial_j = initial_j - latent*sum(state%excess)
        final_j = final_j - latent*sum(excess)
    endif
    naive_delta_j = final_j-initial_j
end subroutine

subroutine monitor_heat_step(stats, unit, stage, end_seconds, dt_seconds, delta_j, storage_j, naive_delta_j, ledger)
    type(HeatStepStats), intent(inout) :: stats
    integer, intent(in) :: unit
    character(len=*), intent(in) :: stage
    real(kind=JPRD), intent(in) :: end_seconds, dt_seconds, delta_j, storage_j, naive_delta_j
    type(HeatStepLedger), intent(in) :: ledger
    real(kind=JPRD) :: q(4), v(NFIELD)
    logical :: valid(NMETRIC), exchange_valid, storage_valid
    integer :: i, k
    q = ledger%value + ledger%correction
    v = 0.0_JPRD
    v(1:7) = [end_seconds,dt_seconds,delta_j,q]
    v(8) = step_sum([delta_j,-q(1)])
    v(9) = step_sum([delta_j,-q(1),q(3)])
    v(10) = storage_j
    ! Zero denominators have no ratio. No arbitrary 1 J floor or silent division by zero.
    exchange_valid = q(2) > 0.0_JPRD
    if (exchange_valid .and. q(2) < 1.0_JPRD) &
    &   exchange_valid = max(abs(v(8)),abs(v(9)),abs(q(4))) <= huge(1.0_JPRD)*q(2)
    storage_valid = storage_j > 0.0_JPRD
    if (storage_valid .and. storage_j < 1.0_JPRD) &
    &   storage_valid = abs(v(9)) <= huge(1.0_JPRD)*storage_j
    if (exchange_valid) v(11:13) = [v(8),v(9),q(4)]/q(2)
    if (storage_valid) v(14) = v(9)/storage_j
    v(15) = naive_delta_j
    v(16) = step_sum([naive_delta_j,-q(1),q(3)])
    stats%steps = stats%steps + 1_JPIB
    if (.not. exchange_valid) stats%no_exchange_ratio = stats%no_exchange_ratio + 1_JPIB
    if (.not. storage_valid) stats%no_storage_ratio = stats%no_storage_ratio + 1_JPIB
    valid = .true.
    valid(5:7) = exchange_valid
    valid(8) = storage_valid
    do i = 1, NMETRIC
        if (.not. valid(i)) cycle
        k = metric_field(i)
        if (stats%samples(i) == 0 .or. v(k) < stats%min_record(k,i)) then
            stats%min_record(:,i) = v
            stats%min_step(i) = stats%steps
        endif
        if (stats%samples(i) == 0 .or. v(k) > stats%max_record(k,i)) then
            stats%max_record(:,i) = v
            stats%max_step(i) = stats%steps
        endif
        stats%samples(i) = stats%samples(i) + 1_JPIB
    enddo
    write(unit,'(a,1x,a,1x,i0,16(1x,es24.16),2(1x,l1))') &
    &   'HEAT_STEP',stage,stats%steps,v,exchange_valid,storage_valid
end subroutine

subroutine write_heat_step_extrema(unit, stage, stats)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: stage
    type(HeatStepStats), intent(in) :: stats
    integer :: i
    write(unit,'(a,1x,a,3(1x,i0))') 'HEAT_STEP_COUNTS',stage,stats%steps, &
    &   stats%no_exchange_ratio,stats%no_storage_ratio
    do i = 1, NMETRIC
        if (stats%samples(i) == 0) cycle
        write(unit,'(a,3(1x,a),1x,i0,16(1x,es24.16))') &
        &   'HEAT_STEP_EXTREME',stage,metric_name(i),'min',stats%min_step(i),stats%min_record(:,i)
        write(unit,'(a,3(1x,a),1x,i0,16(1x,es24.16))') &
        &   'HEAT_STEP_EXTREME',stage,metric_name(i),'max',stats%max_step(i),stats%max_record(:,i)
    enddo
end subroutine
end module heat_step_monitor_mod
