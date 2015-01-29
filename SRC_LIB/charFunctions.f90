module charFunctions

    implicit none

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function string_vec_join(stringVec) result(stringTot)

        implicit none

        !INPUT
        character (len=*), dimension(:), intent(in) :: stringVec;

        !OUTPUT
        character (len=200) :: stringTot;

        !LOCAL
        integer :: i

        !write(*,*) "WRITE Flag string_join"
        stringTot = ""

        do i = 1, size(stringVec)
            stringTot = string_join(stringTot, stringVec(i))
        end do

        stringTot = adjustL(stringTot)

    end function string_vec_join

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function string_join(string1, string2) result(stringTot)

        implicit none

        !INPUT
        character (len=*)  , intent(in) :: string1, string2;

        !OUTPUT
        character (len=100) :: stringTot;

        !write(*,*) "WRITE Flag string_join"

        stringTot = trim(adjustL(string1))//trim(adjustL(string2))
        stringTot = adjustL(stringTot)

    end function string_join

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function stringNumb_join(string, number) result(stringTot)

        implicit none

        !INPUT
        character (len=*)  , intent(in) :: string
        integer            , intent(in) :: number

        !OUTPUT
        character (len=100) :: stringTot;

        !LOCAL
        character (len=30)  :: nmbString

        !write(*,*) "WRITE Flag stringNumb_join"

        write(nmbString, fmt='(I8)') number
        stringTot = string_join(string, nmbString)

        !write(*,*) "WRITE Flag 2 stringNumb_join"

    end function

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function numb2String(number) result(stringTot)

        implicit none

        !INPUT
        integer, intent(in) :: number

        !OUTPUT
        character (len=30) :: stringTot;

        !LOCAL
        character (len=30)  :: nmbString

        write(nmbString, fmt='(I8)') number
        stringTot = adjustL(nmbString)

    end function

end module charFunctions
