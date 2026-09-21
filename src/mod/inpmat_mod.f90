module inpmat_mod
    use PARKIND1, only: &
    &   JPRM, JPRB, JPIM, JPIB
    use YOS_CMF_INPUT, only: &
    &   TMPNAM, LOGNAM, &
    &   RMIS, DMIS, &
    &   CSETFILE
    use CMF_UTILS_MOD, only: &
    &   INQUIRE_FID
    use YOS_CMF_INPUT, only: &
    &   NX, NY, NLFP
    use YOS_CMF_MAP, only: &
    &   NSEQMAX, I1SEQX, I1SEQY

    use mapframe_mod, only: &
    &   MapFrame
    use glob_mod, only: &
    &   MNTSEQ
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: Inpmat, append_inpmat, load_inpmat_files, move_append_inpmat

    type Inpmat
        type(MapFrame) :: &
        &   map
        integer :: &
        &   inpn
        integer, allocatable :: &
        &   inpx(:,:,:), inpy(:,:,:)
        real(kind=JPRB), allocatable :: &
        &   inpa(:,:,:)

        contains

        procedure :: map2vec_intrp => map2vec_intrp
        procedure :: matches_shape

    end type Inpmat

    interface Inpmat
        module procedure init_inpmat
    end interface Inpmat

contains

! ===================================================================================================
! Constructor
! ===================================================================================================
function init_inpmat(item) result(obj)
    type(Inpmat) obj
    character(len=*), intent(in) :: item
    character(len=256) :: dir, prefix
    integer :: nxin, nyin
    real(kind=JPRB) :: west, east, south, north
    logical :: is_n2s
    call read_inpmat_nml( &
    &   item, &
    &   dir, prefix)
    call read_inpmat_dim( &
    &   dir, prefix, nxin, nyin, obj%inpn)
    call obj%map%set_shape( &
    &   nxin, nyin)
    call read_inpmat_area( &
    &   dir, prefix, west, east, south, north, is_n2s)
    if (is_n2s) then
        call obj%map%set_domain(west, east, north, south)
    else
        call obj%map%set_domain(west, east, south, north)
    endif
    call read_inpmat_bin( &
    &   dir, prefix, obj%inpx, obj%inpy, obj%inpa, obj%inpn)
    write(LOGNAM, '(2a,i3)') '          ', trim(obj%map%str()), obj%inpn
!    if (.not. allocated(grda)) then
!        allocate(grda(NSEQMAX)); grda(:) = 0.d0
!    endif
    !call calc_gridarea(ainpmat%grda, &
    !&    ainpmat%inpx, ainpmat%inpy, ainpmat%inpa, ainpmat%nxin, ainpmat%nyin, ainpmat%inpn)

    contains

    subroutine read_inpmat_nml(item_name, dir, prefix)
        character(len=*)  , intent(in)  :: item_name
        character(len=256), intent(out) :: dir, prefix
        character(len=256) :: item
        namelist /nml_inpmat/ item, dir, prefix
        integer :: ios

        TMPNAM = INQUIRE_FID()
        open(TMPNAM, file=trim(CSETFILE), status='old', iostat=ios)
        if (ios /= 0) stop 'read_inpmat_nml ERROR : open namelist'
        do
            read(TMPNAM, nml=nml_inpmat, iostat=ios)
            if (ios < 0) stop 'read_inpmat_nml ERROR : not found item'
            if (trim(item) == trim(item_name)) exit
        enddo
        close(TMPNAM)
        write(LOGNAM, '(6a)') '      ', trim(item), ' : ', trim(dir), '/', trim(prefix)
    end subroutine read_inpmat_nml


    function path_inpmat(dir, prefix, suffix) result(path)
        character(len=256) :: path
        character(len=*)   :: dir, prefix, suffix
        path = trim(dir)//'/'//trim(prefix)//trim(suffix)
    end function path_inpmat


    subroutine read_inpmat_dim(dir, prefix, nxin, nyin, inpn)
        character(len=*), intent(in)  :: dir, prefix
        integer         , intent(out) :: nxin, nyin, inpn
        character(len=256) :: path

        TMPNAM = INQUIRE_FID()
        path = path_inpmat(dir, prefix, '.dim')
        open(TMPNAM, file=trim(path), form='formatted')
        read(TMPNAM, *) nxin
        read(TMPNAM, *) nyin
        read(TMPNAM, *) inpn
        close(TMPNAM)
        !write(LOGNAM, '(a,2i5,i3)') '  nxin, nyin, inpn =', nxin, nyin, inpn
    end subroutine read_inpmat_dim


    subroutine read_inpmat_area(dir, prefix, west, east, south, north, is_n2s)
        character(len=*), intent(in)  :: dir, prefix
        real(kind=JPRB),  intent(out) :: west, east, south, north
        logical,          intent(out) :: is_n2s
        character(len=256) :: path
        real(kind=JPRM) :: west_file, east_file, south_file, north_file

        TMPNAM = INQUIRE_FID()
        path = path_inpmat(dir, prefix, '.cfg')
        open(TMPNAM, file=trim(path), form='formatted')
        read(TMPNAM, *)
        read(TMPNAM, *) west_file
        read(TMPNAM, *) east_file
        read(TMPNAM, *) south_file
        read(TMPNAM, *) north_file
        read(TMPNAM, *) is_n2s
        close(TMPNAM)
        west  = real(west_file , kind=JPRB)
        east  = real(east_file , kind=JPRB)
        south = real(south_file, kind=JPRB)
        north = real(north_file, kind=JPRB)
        !write(LOGNAM, '(a,4f0.2,L)') '  west, east, south, north =', west, east, south, north, is_n2s
    end subroutine read_inpmat_area


    subroutine read_inpmat_bin(dir, prefix, inpx, inpy, inpa, inpn)
        character(len=*),              intent(in)  :: dir, prefix
        integer,          allocatable, intent(out) :: inpx(:,:,:), inpy(:,:,:)
        real(kind=JPRB),  allocatable, intent(out) :: inpa(:,:,:)
        integer         ,              intent(in)  :: inpn
        real(kind=JPRM), allocatable :: rinp(:,:,:)
        integer, parameter :: byte_recl = 4
        character(len=256) :: path

        TMPNAM = INQUIRE_FID()
        path = path_inpmat(dir, prefix, '.bin')
        allocate(inpx(NX,NY,inpn), inpy(NX,NY,inpn), inpa(NX,NY,inpn), rinp(NX,NY,inpn))
        open(TMPNAM, file=trim(path), form='unformatted', access='direct', &
        &    recl=byte_recl*size(inpx))
        read(TMPNAM, rec=1) inpx(:,:,:)
        read(TMPNAM, rec=2) inpy(:,:,:)
        read(TMPNAM, rec=3) rinp(:,:,:)
        close(TMPNAM)
        inpa(:,:,:) = real(rinp(:,:,:), kind=JPRB)
        deallocate(rinp)
    end subroutine read_inpmat_bin

end function init_inpmat

! ===================================================================================================
! Interpolation with inpmat
! ===================================================================================================
subroutine map2vec_intrp(self, map, vec)
    class(Inpmat)   , intent(in)  :: self
    real(kind=JPRM), intent(in)  :: map(:,:)
    real(kind=JPRB), intent(out) :: vec(:)
    real(kind=JPRM), parameter :: fill_val = 1.e16_JPRM
    integer :: iseq, ix, iy, inpi, ixin, iyin
    real(kind=JPRB) :: grda

    !$omp parallel do private(ix, iy, inpi, ixin, iyin, grda)
    do iseq = 1, NSEQMAX
        vec(iseq) = 0.0_JPRB
        ix = I1SEQX(iseq)
        iy = I1SEQY(iseq)
        grda = 0.0_JPRB
        do inpi = 1, self%INPN
            if ( self%inpa(ix,iy,inpi) <= 0.0_JPRB ) exit
            ixin = self%INPX(ix, iy, inpi)
            iyin = self%INPY(ix, iy, inpi)
            if (ixin < 1) exit
            if (isnan(map(ixin,iyin))) cycle
            if ( map(ixin,iyin) < 0.0_JPRM ) cycle
            if ( map(ixin,iyin) >= fill_val ) cycle
            vec(iseq) = vec(iseq) + real(map(ixin,iyin), kind=JPRB) * self%inpa(ix,iy,inpi)
            grda      = grda      +                  self%inpa(ix,iy,inpi)
        enddo
        if (grda > 0.0_JPRB) vec(iseq) = vec(iseq) / grda
    enddo
    !$omp end parallel do
end subroutine map2vec_intrp

! ===================================================================================================
subroutine append_inpmat(inpmats, ainpmat)
    type(Inpmat), allocatable, intent(inout) :: inpmats(:)
    type(Inpmat),              intent(in)    :: ainpmat
    type(Inpmat), allocatable                :: tmp(:)
    integer :: old_size

    if (.not. allocated(inpmats)) then
        allocate(inpmats(1)); inpmats(1) = ainpmat
        return
    endif
    old_size = size(inpmats)
    allocate(tmp, source=inpmats)
    deallocate(inpmats)
    allocate(inpmats(old_size + 1))
    inpmats(1:old_size) = tmp(:)
    inpmats(old_size+1) = ainpmat
    deallocate(tmp)
end subroutine append_inpmat

subroutine load_inpmat_files(obj, diminfo_file, inpmat_file)
    type(Inpmat), intent(out) :: obj
    character(len=*), intent(in) :: diminfo_file, inpmat_file
    integer(kind=JPIM) :: nx_dst, ny_dst, nlfp_dst, nxin, nyin
    real(kind=JPRB) :: bounds(4)

    call read_inpmat_diminfo(diminfo_file, nx_dst, ny_dst, nlfp_dst, nxin, nyin, obj%inpn, bounds)
    if (nx_dst /= NX .or. ny_dst /= NY .or. nlfp_dst /= NLFP) then
        write(LOGNAM, '(a)') '[load_inpmat_files ERROR] destination shape mismatch'
        stop 1
    endif
    if (.not. all(ieee_is_finite(bounds))) then
        write(LOGNAM, '(a)') '[load_inpmat_files ERROR] non-finite domain bounds'
        stop 1
    endif
    call obj%map%set_shape(nxin, nyin)
    call obj%map%set_domain(bounds(1), bounds(2), bounds(3), bounds(4))
    call read_inpmat_file(inpmat_file, obj%inpx, obj%inpy, obj%inpa, obj%inpn)
    call validate_inpmat(obj)
end subroutine load_inpmat_files

subroutine read_inpmat_diminfo(path, nx_dst, ny_dst, nlfp_dst, nxin, nyin, inpn, bounds)
    character(len=*), intent(in) :: path
    integer(kind=JPIM), intent(out) :: nx_dst, ny_dst, nlfp_dst, nxin, nyin, inpn
    real(kind=JPRB), intent(out) :: bounds(4)
    character(len=512) :: binary_name
    integer :: unit, ios

    unit = INQUIRE_FID()
    open(unit, file=trim(path), form='formatted', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(LOGNAM, '(2a)') '[inpmat_mod/read_inpmat_diminfo ERROR] cannot open: ', trim(path)
        stop 1
    endif
    read(unit, *, iostat=ios) nx_dst
    if (ios == 0) read(unit, *, iostat=ios) ny_dst
    if (ios == 0) read(unit, *, iostat=ios) nlfp_dst
    if (ios == 0) read(unit, *, iostat=ios) nxin
    if (ios == 0) read(unit, *, iostat=ios) nyin
    if (ios == 0) read(unit, *, iostat=ios) inpn
    ! The embedded name is descriptive; the caller supplies the binary path.
    if (ios == 0) read(unit, '(a)', iostat=ios) binary_name
    if (ios == 0) read(unit, *, iostat=ios) bounds(1)
    if (ios == 0) read(unit, *, iostat=ios) bounds(2)
    if (ios == 0) read(unit, *, iostat=ios) bounds(3)
    if (ios == 0) read(unit, *, iostat=ios) bounds(4)
    close(unit)
    if (ios /= 0) then
        write(LOGNAM, '(2a)') '[inpmat_mod/read_inpmat_diminfo ERROR] invalid file: ', trim(path)
        stop 1
    endif
    if (nx_dst < 1 .or. ny_dst < 1 .or. nlfp_dst < 1 .or. nxin < 1 .or. nyin < 1 .or. inpn < 1) then
        write(LOGNAM, '(2a)') '[inpmat_mod/read_inpmat_diminfo ERROR] non-positive dimension: ', trim(path)
        stop 1
    endif
end subroutine read_inpmat_diminfo


subroutine read_inpmat_file(path, inpx, inpy, inpa, inpn)
    character(len=*), intent(in) :: path
    integer(kind=JPIM), allocatable, intent(out) :: inpx(:,:,:), inpy(:,:,:)
    real(kind=JPRB), allocatable, intent(out) :: inpa(:,:,:)
    integer(kind=JPIM), intent(in) :: inpn
    real(kind=JPRM), allocatable :: rinp(:,:,:)
    integer, parameter :: byte_recl = 4
    integer :: unit, ios
    integer(kind=JPIB) :: file_size, expected_size
    logical :: exists

    inquire(file=trim(path), exist=exists, size=file_size, iostat=ios)
    if (ios /= 0 .or. .not. exists) then
        write(LOGNAM, '(2a)') '[inpmat_mod/read_inpmat_file ERROR] cannot stat: ', trim(path)
        stop 1
    endif
    expected_size = 3_JPIB * int(NX, JPIB) * int(NY, JPIB) * int(inpn, JPIB) * byte_recl
    if (file_size /= expected_size) then
        write(LOGNAM, '(2a)') '[inpmat_mod/read_inpmat_file ERROR] incorrect file size: ', trim(path)
        write(LOGNAM, '(a,i0)') '    actual:   ', file_size
        write(LOGNAM, '(a,i0)') '    expected: ', expected_size
        stop 1
    endif

    allocate(inpx(NX,NY,inpn), inpy(NX,NY,inpn), inpa(NX,NY,inpn), rinp(NX,NY,inpn))
    unit = INQUIRE_FID()
    open(unit, file=trim(path), form='unformatted', access='stream', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(LOGNAM, '(2a)') '[inpmat_mod/read_inpmat_file ERROR] cannot open: ', trim(path)
        stop 1
    endif
    read(unit, iostat=ios) inpx(:,:,:)
    if (ios == 0) read(unit, iostat=ios) inpy(:,:,:)
    if (ios == 0) read(unit, iostat=ios) rinp(:,:,:)
    close(unit)
    if (ios /= 0) then
        write(LOGNAM, '(2a)') '[inpmat_mod/read_inpmat_file ERROR] cannot read: ', trim(path)
        stop 1
    endif
    inpa(:,:,:) = real(rinp(:,:,:), kind=JPRB)
    deallocate(rinp)
end subroutine read_inpmat_file


subroutine validate_inpmat(self)
    class(Inpmat), intent(in) :: self
    integer(kind=JPIM) :: inpi

    if (any(self%inpx < 0) .or. any(self%inpx > self%map%nx)) then
        write(LOGNAM, '(a)') '[inpmat_mod/validate_inpmat ERROR] inpx is out of range'
        stop 1
    endif
    if (any(self%inpy < 0) .or. any(self%inpy > self%map%ny)) then
        write(LOGNAM, '(a)') '[inpmat_mod/validate_inpmat ERROR] inpy is out of range'
        stop 1
    endif
    if (any(.not. ieee_is_finite(self%inpa)) .or. any(self%inpa < 0.0_JPRB)) then
        write(LOGNAM, '(a)') '[inpmat_mod/validate_inpmat ERROR] inpa must be finite and non-negative'
        stop 1
    endif
    if (any((self%inpx > 0) .neqv. (self%inpy > 0)) .or. &
    &   any((self%inpx > 0) .neqv. (self%inpa > 0.0_JPRB))) then
        write(LOGNAM, '(a)') '[inpmat_mod/validate_inpmat ERROR] inconsistent zero padding'
        stop 1
    endif
    do inpi = 2, self%inpn
        if (any((self%inpa(:,:,inpi) > 0.0_JPRB) .and. (self%inpa(:,:,inpi-1) <= 0.0_JPRB))) then
            write(LOGNAM, '(a)') '[inpmat_mod/validate_inpmat ERROR] non-contiguous input slots'
            stop 1
        endif
    enddo
end subroutine validate_inpmat


logical function matches_shape(self, nxin, nyin) result(matches)
    class(Inpmat), intent(in) :: self
    integer(kind=JPIM), intent(in) :: nxin, nyin
    matches = self%map%nx == nxin .and. self%map%ny == nyin
end function matches_shape

subroutine move_append_inpmat(inpmats, ainpmat)
    type(Inpmat), allocatable, intent(inout) :: inpmats(:)
    type(Inpmat),              intent(inout) :: ainpmat
    type(Inpmat), allocatable                :: tmp(:)
    integer :: i, old_size

    if (.not. allocated(inpmats)) then
        allocate(inpmats(1))
        call move_inpmat(ainpmat, inpmats(1))
        return
    endif
    old_size = size(inpmats)
    allocate(tmp(old_size + 1))
    do i = 1, old_size
        call move_inpmat(inpmats(i), tmp(i))
    enddo
    call move_inpmat(ainpmat, tmp(old_size + 1))
    call move_alloc(tmp, inpmats)
end subroutine move_append_inpmat


subroutine move_inpmat(source, destination)
    type(Inpmat), intent(inout) :: source
    type(Inpmat), intent(out) :: destination

    destination%map = source%map
    destination%inpn = source%inpn
    call move_alloc(source%inpx, destination%inpx)
    call move_alloc(source%inpy, destination%inpy)
    call move_alloc(source%inpa, destination%inpa)
    source%inpn = 0_JPIM
end subroutine move_inpmat

end module inpmat_mod
