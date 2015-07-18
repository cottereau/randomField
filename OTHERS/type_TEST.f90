module type_TEST

    use mpi

    implicit none

    type :: TEST

        logical :: init = .false.

    end type TEST

    contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_TEST(TEST_a, name)
            type(TEST), intent(in) :: TEST_a
            character(len=*), intent(in), optional :: name

            if(TEST_a%init) then
                write(*,*) "TEST-----------------------"
                if(present(name)) write(*,*) name
                write(*,*) "    init     = ", TEST_a%init
                write(*,*) ""
            else
                write(*,*) "TEST has not been initialized----"
                write(*,*) "    init     = ", TEST_a%init
                write(*,*) ""
            end if

        end subroutine show_TEST

end module type_TEST
