module systemUt_RF
    use mpi
    use write_Log_File

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine create_folder(folder, path, rang, comm)

        implicit none
        !INPUT
        character(len = *), intent(in) :: folder, path
        integer, intent(in) :: rang, comm
        !LOCAL
        character(len = 200) command, fullName
        integer :: code
        logical :: dirExists

        if(rang == 0) then
            fullName = trim(adjustL(path)) // "/" // trim(adjustL(folder))
            !write(*,*) "fullName = ", fullName
            write(*,*) "Directory is being created: ", fullName
            command = 'mkdir '// trim(adjustL(fullName))
            call system(command)
        end if

        call MPI_BARRIER (comm ,code)

    end subroutine create_folder

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine delete_folder(folder, path, rang, comm)
        implicit none
        !INPUT
        character(len = *), intent(in) :: folder, path
        integer, intent(in) :: rang, comm
        !LOCAL
        character(len = 200) command, fullName

        integer :: code

        fullName = trim(adjustL(path)) // "/" // trim(adjustL(folder))

        if(rang == 0) then
            write (*,*) "WARNING!!! Directory: '", trim(fullName) ,"' will be deleted"
            command = 'rm -r '// trim(adjustL(fullName))
            call system(command)
        end if

        call MPI_BARRIER (comm ,code)

    end subroutine delete_folder

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine rename_folder(folder, path, newName, rang, comm)
        implicit none
        !INPUT
        character(len = *), intent(in) :: folder, path, newName
        integer, intent(in) :: rang, comm
        !LOCAL
        character(len = 200) command, fullName_new, fullName_old
        integer :: compiler = 2; !1 for gfortran and 2 for ifort
        integer :: code

        fullName_old = trim(adjustL(path)) // "/" // trim(adjustL(folder))
        fullName_new = trim(adjustL(path)) // "/" // trim(adjustL(newName))

        if(rang == 0) then
            write (*,*) "WARNING!!! Directory: '", trim(fullName_old) ,"' will be renamed to: ", trim(fullName_new)
            command = 'mv '// trim(adjustL(fullName_old))//" "//trim(adjustL(fullName_new))
            call system(command)
        end if

        call MPI_BARRIER (comm ,code)

    end subroutine rename_folder

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine check_folder_existence(folder, path, compiler, dirExists)
        implicit none
        !INPUT
        character(len = *), intent(in) :: folder, path
        integer, intent(in) :: compiler !1 for gfortran and 2 for ifort
        !OUTPUT
        logical :: dirExists
        !LOCAL
        character(len = 200) fullName

        fullName = trim(adjustL(path)) // "/" // trim(adjustL(folder))

        if(compiler == 1) inquire( file=trim(fullName)//'/.', exist=dirExists )  ! Works with gfortran, but not ifort
        if(compiler == 2) inquire( directory=fullName, exist=dirExists )         ! Works with ifort, but not gfortran

    end subroutine check_folder_existence

end module systemUt_RF


!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
