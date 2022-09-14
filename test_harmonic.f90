program test_harmonic

    use harmonics
    implicit none

    call print_human(123_long)
    call print_human(234000_long)
    call print_human(563400000_long)
    call print_human(145300000_long)
    call print_human(999900000000_long)
    call print_human(999300000000000_long)

contains

    subroutine print_human(val)

        integer(long), intent(in) :: val
        character(15) :: char

        char = ""

        call human_format(char, val)
        print *, char

    end subroutine

end program