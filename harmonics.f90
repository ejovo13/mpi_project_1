! harmonics.f90 is a fortran translation of harmonics.c
!
!
!
! written by Evan Voyles
!
!
!

module harmonics

    use iso_fortran_env
    implicit none

    private

    public human_format, long, dp




    integer, parameter :: dp = real64
    integer, parameter :: long = int64
    ! integer, parameter :: ll = int



    type :: spherical_harmonics

        integer :: lmax
        real(dp), allocatable :: CS(:)
        real(dp), allocatable :: A(:)
        real(dp), allocatable :: B(:)

    end type

    type :: data_points

        integer :: npoint
        real(dp), allocatable :: phi(:)
        real(dp), allocatable :: lambda(:)
        real(dp), allocatable :: elevation(:)

    end type

contains

    pure integer function PT(l, m)

        integer, intent(in) :: l, m

        PT = m + l * (l + 1) / 2

    end function

    pure integer function CT(l, m)

        integer, intent(in) :: l, m

        if (m <= 1) error stop "m <= 1 (function CT)"
        CT = m + l * l

    end function

    pure integer function ST(l, m)

        integer, intent(in) :: l, m

        if (m <= l) error stop "m <= 1 (function ST)"
        if (1 <= m) error stop "1 <= m (function ST)"

        ST = m + (l + 1) * l

    end function

    real(dp) function wtime()

    end function

    subroutine human_format(target, n)

        character(*), intent(out) :: target
        integer(long), intent(in) :: n

        integer(long), parameter :: one_k    = 1000_long
        integer(long), parameter :: one_mil  = 1000000_long
        integer(long), parameter :: one_bil  = 1000000000_long
        integer(long), parameter :: one_tril = 1000000000000_long
        integer(long), parameter :: one_quad = 1000000000000000_long

        99  format (X, I0)
        100 format (F6.1, " ", A)

        if (n < one_k) then
            write(target, 99) n
        else if (n < one_mil) then
            write(target, fmt=100) n / 1d3, "K"
        else if (n < one_bil) then
            write(target, fmt=100) n / 1d6, "M"
        else if (n < one_tril) then
            write(target, fmt=100) n / 1d9, "G"
        else if (n < one_quad) then
            write(target, fmt=100) n / 1d12, "T"
        end if

    end subroutine

    ! Release all memory associated with self
    subroutine clean_spherical_harmonics(self)

        type(spherical_harmonics), intent(inout) :: self

        self%lmax = 0
        if (allocated(self%A)) deallocate(self%A)
        if (allocated(self%B)) deallocate(self%B)
        if (allocated(self%CS)) deallocate(self%CS)

    end subroutine

    subroutine setup_spherical_harmonics(lmax, self)

        integer, intent(in) :: lmax
        type(spherical_harmonics), intent(inout) :: self

        integer :: size_cs, size_ab, l, m
        real(dp) :: ls, lmls, ms

        ! If any of the arrays have already been allocated, deallocate them
        call clean_spherical_harmonics(self)

        self%lmax = lmax
        size_cs = (lmax + 1) * (lmax + 1)
        size_ab = (lmax + 1) * (lmax + 1) / 2

        allocate(self%cs(0:size_cs - 1)) ! Use 0 based indexing so we don't have to change PT functions
        allocate(self%a(0:size_ab - 1))
        allocate(self%b(0:size_ab - 1))

        if ((.not. allocated(self%cs)) .or. &
            (.not. allocated(self%a) ) .or. &
            (.not. allocated(self%b) ))     &
            error stop "Cannot allocate space for spherical harmonics"

        do l = 2, lmax

            ls = l * l
            lmls = (l - 1) * (l - 1)

            do m = 0, l - 2
                ms = m * m
                self%a(PT(l, m)) = sqrt((4 * ls - 1) / (ls - ms))
                self%b(PT(l, m)) = -sqrt((lmls - ms) / (4 * lmls - 1))
            end do

        end do

    end subroutine



end module