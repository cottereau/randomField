module type_MESH

    use mpi

    implicit none

    type :: MESH
        !MPI VARIABLES
        integer :: comm = -1
        integer :: rang = -1
        integer :: nb_procs = -1

        !MESH VARIABLES
            !nDim independent
        character (len=15) :: meshType, meshMod;
        integer :: nDim = -1, xNTotal = -1;
            !nDim dependent
        integer         , dimension(:), allocatable :: xNStep;
        double precision, dimension(:), allocatable :: xMax, xMin, xStep;
        double precision, dimension(:), allocatable :: xMaxGlob, xMinGlob;
        logical :: init = .false.

    end type MESH

    contains
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_MESH(MESH_a, nDim, comm, rang, nb_procs)
            type(MESH) :: MESH_a
            integer :: nDim, comm, rang, nb_procs

            MESH_a%nDim     = nDim
            MESH_a%comm     = comm
            MESH_a%rang     = rang
            MESH_a%nb_procs = nb_procs
            allocate(MESH_a%xMax(nDim))
            allocate(MESH_a%xMin(nDim))
            allocate(MESH_a%xStep(nDim))
            allocate(MESH_a%xNStep(nDim))
            allocate(MESH_a%xMaxGlob(nDim))
            allocate(MESH_a%xMinGlob(nDim))
            MESH_a%xMax     = -1
            MESH_a%xMin     = -1
            MESH_a%xStep    = -1
            MESH_a%xNStep   = -1
            MESH_a%xMaxGlob = -1
            MESH_a%xMinGlob = -1

            MESH_a%init = .true.

        end subroutine init_MESH

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_MESH(MESH_a, name, fmt)
            type(MESH), intent(in) :: MESH_a
            character(len=*), intent(in), optional :: name
            character(len = 20), intent(in), optional :: fmt
            character(len = 20) :: dblFmt

            dblFmt = "T20,F15.5"
            if(present(fmt)) dblFmt = fmt

            write(*,*) "MESH------------------------------------------------------------"
            if(present(name)) write(*,*) "|  ", name

            if(MESH_a%init) then
                write(*,*) "|  init     = ", MESH_a%init
                write(*,*) "|"
                write(*,*) "|  MPI---"
                write(*,*) "|  |rang     = ", MESH_a%rang
                write(*,*) "|  |comm     = ", MESH_a%comm
                write(*,*) "|  |nb_procs = ", MESH_a%nb_procs
                write(*,*) "|"
                write(*,*) "|  Input---"
                write(*,*) "|  |nDim     = ", MESH_a%nDim
                write(*,*) "|  |meshType = ", MESH_a%meshType
                write(*,*) "|  |meshMod  = ", MESH_a%meshMod
                write(*,"(A,("//dblFmt//"))") " |  |xMinGlob = ", MESH_a%xMinGlob
                write(*,"(A,("//dblFmt//"))") " |  |xMaxGlob = ", MESH_a%xMaxGlob
                write(*,"(A,("//dblFmt//"))") " |  |xStep    = ", MESH_a%xStep
                write(*,*) "|"
                write(*,*) "|  |Process--"
                write(*,*) "|  |NStep    = ", MESH_a%xNStep
                write(*,*) "|  |xNTotal  = ", MESH_a%xNTotal
                write(*,"(A,("//dblFmt//"))") " |  |xMin     = ", MESH_a%xMin
                write(*,"(A,("//dblFmt//"))") " |  |xMax     = ", MESH_a%xMax
            else
                write(*,*) "|  init     = ", MESH_a%init
                write(*,*) "|  MESH has not been initialized----"
            end if
            write(*,*) "|---------------------------------------------------------------"
            write(*,*) ""

        end subroutine show_MESH

!---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine finalize_MESH(MESH_a)
            type(MESH) :: MESH_a

            if (allocated(MESH_a%xMax))     deallocate(MESH_a%xMax)
            if (allocated(MESH_a%xMin))     deallocate(MESH_a%xMin)
            if (allocated(MESH_a%xStep))    deallocate(MESH_a%xStep)
            if (allocated(MESH_a%xNStep))   deallocate(MESH_a%xNStep)
            if (allocated(MESH_a%xMaxGlob)) deallocate(MESH_a%xMaxGlob)
            if (allocated(MESH_a%xMinGlob)) deallocate(MESH_a%xMinGlob)

            MESH_a%init = .false.

        end subroutine finalize_MESH

end module type_MESH
