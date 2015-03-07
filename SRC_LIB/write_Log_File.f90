module write_Log_File

    use charFunctions

    implicit none

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine init_log_file(filename, rang)
        implicit none
        !INPUT
        character(len=*), intent(in) :: filename
        integer         , intent(in) :: rang
        !LOCAL
        integer :: fileId

        fileId = get_fileId (rang)
        !write(*,*) "fileId INSIDE = ", fileId
        !write(*,*) "filename INSIDE = ", filename

        open (unit = fileId , file = filename, action = 'write')

    end subroutine init_log_file

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine finalize_log_file(rang)
        implicit none
        !INPUT
        integer, intent(in) :: rang
        !LOCAL
        integer :: fileId

        fileId = get_fileId()

        close(fileId)

    end subroutine finalize_log_file

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function get_fileId (rang) result(fileId)
        !INPUT
        integer, intent(in), optional :: rang
        integer, save :: counter = 0, saved_Rang = -1
        !OUTPUT
        integer :: fileId
        !LOCAL
        integer :: cte = 50

        !write(*,*) "Inside get_fileId"

        if(counter == 0 .and. present(rang)) then
            saved_Rang = rang
            counter = counter + 1
            fileId = cte + saved_Rang
        else if(.not. saved_Rang == -1) then
            fileId = cte + saved_Rang
        else if (present(rang)) then
            write(*,*) "WARNING!!! get_fileId tryed to initialize rang twice"
            write(*,*) "old value = ", saved_Rang
            write(*,*) "new_value (not attributed) = ", rang
        else
            write(*,*) "ERROR!!! rang is a mandatory variable on the first call of get_fileId"
            stop
        end if

    end function get_fileId

end module write_Log_File
