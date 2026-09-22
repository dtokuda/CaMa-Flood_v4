module heat_residual_mod
    use PARKIND1, only: JPIM, JPIB, JPRD
    implicit none
    private
    public :: HeatResidualStats, record_heat_residual, write_heat_residual
    public :: HeatConservationStats, record_heat_exchange, write_heat_conservation

    ! Diagnostics only: no value in this type is returned to a physical state.
    ! Accumulation starts at process initialization, including for restart runs.
    type HeatResidualStats
        real(kind=JPRD) :: positive_j = 0.0_JPRD, negative_j = 0.0_JPRD
        real(kind=JPRD) :: maximum_absolute_j = 0.0_JPRD, throughput_j = 0.0_JPRD
        integer(kind=JPIB) :: events = 0_JPIB
        integer(kind=JPIM) :: maximum_cell = 0_JPIM
        logical, allocatable :: affected(:)
    end type HeatResidualStats
    type HeatConservationStats
        logical :: initialized = .false.
        real(kind=JPRD) :: initial_j = 0.0_JPRD
        ! Water boundaries, ice boundaries, then local heat input. Kahan sums across steps.
        real(kind=JPRD) :: net_j(3) = 0.0_JPRD, absolute_j(3) = 0.0_JPRD
        real(kind=JPRD) :: net_correction_j(3) = 0.0_JPRD, absolute_correction_j(3) = 0.0_JPRD
        integer(kind=JPIB) :: steps(3) = 0_JPIB
    end type HeatConservationStats
contains
subroutine record_heat_exchange(stats, process, net_j, absolute_j)
    type(HeatConservationStats), intent(inout) :: stats
    integer(kind=JPIM), intent(in) :: process
    real(kind=JPRD), intent(in) :: net_j, absolute_j
    call compensated_add(stats%net_j(process), stats%net_correction_j(process), net_j)
    call compensated_add(stats%absolute_j(process), stats%absolute_correction_j(process), absolute_j)
    stats%steps(process) = stats%steps(process) + 1_JPIB
end subroutine

subroutine compensated_add(total, correction, value)
    real(kind=JPRD), intent(inout) :: total, correction
    real(kind=JPRD), intent(in) :: value
    real(kind=JPRD) :: adjusted, updated
    adjusted = value - correction
    updated = total + adjusted
    correction = (updated - total) - adjusted
    total = updated
end subroutine

subroutine write_heat_conservation(unit, stats, current_j, unapplied_j)
    integer(kind=JPIM), intent(in) :: unit
    type(HeatConservationStats), intent(in) :: stats
    real(kind=JPRD), intent(in) :: current_j, unapplied_j
    real(kind=JPRD) :: raw_j, adjusted_j, scale_j
    raw_j = (current_j - stats%initial_j) - sum(stats%net_j)
    adjusted_j = raw_j + unapplied_j
    scale_j = max(sum(stats%absolute_j), 1.0_JPRD)
    write(unit, '(a,13(1x,es24.16),3(1x,i0))') 'HEAT_CONSERVATION', &
    &   stats%initial_j, current_j, stats%net_j, stats%absolute_j, unapplied_j, &
    &   raw_j, adjusted_j, raw_j / scale_j, adjusted_j / scale_j, stats%steps
end subroutine

subroutine record_heat_residual(stats, residual_j, throughput_j, mask)
    type(HeatResidualStats), intent(inout) :: stats
    real(kind=JPRD), intent(in) :: residual_j(:) ! [J] Expected minus represented energy.
    real(kind=JPRD), intent(in) :: throughput_j(:) ! [J] Nonnegative absolute energy input/output scale.
    logical, intent(in), optional :: mask(:)
    integer(kind=JPIM) :: i
    if (.not. allocated(stats%affected)) allocate(stats%affected(size(residual_j)), source=.false.)
    do i = 1, size(residual_j)
        if (present(mask)) then
            if (.not. mask(i)) cycle
        endif
        stats%throughput_j = stats%throughput_j + throughput_j(i)
        stats%positive_j = stats%positive_j + max(residual_j(i), 0.0_JPRD)
        stats%negative_j = stats%negative_j + min(residual_j(i), 0.0_JPRD)
        if (residual_j(i) /= 0.0_JPRD) then
            stats%events = stats%events + 1_JPIB
            stats%affected(i) = .true.
        endif
        if (abs(residual_j(i)) > stats%maximum_absolute_j) then
            stats%maximum_absolute_j = abs(residual_j(i))
            stats%maximum_cell = i
        endif
    enddo
end subroutine record_heat_residual

subroutine write_heat_residual(unit, reason, stats, prefix)
    integer(kind=JPIM), intent(in) :: unit
    character(len=*), intent(in) :: reason
    type(HeatResidualStats), intent(in) :: stats
    character(len=*), intent(in), optional :: prefix
    character(len=16) :: tag
    integer(kind=JPIM) :: cells
    tag = 'HEAT_RESIDUAL'
    if (present(prefix)) tag = prefix
    cells = 0
    if (allocated(stats%affected)) cells = count(stats%affected)
    ! Column names are emitted once by the caller. Negative sum retains its sign.
    write(unit, '(a,1x,a,5(1x,es24.16),3(1x,i0),1x,es24.16)') &
    &   trim(tag), trim(reason), stats%positive_j + stats%negative_j, &
    &   stats%positive_j, stats%negative_j, stats%positive_j - stats%negative_j, &
    &   stats%maximum_absolute_j, stats%events, cells, stats%maximum_cell, stats%throughput_j
end subroutine write_heat_residual
end module heat_residual_mod
