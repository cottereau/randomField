module type_MESH

    use mpi
    use charFunctions

    implicit none

    type :: MESH
        !MPI VARIABLES
        integer :: comm = -1
        integer :: rang = -1
        integer :: nb_procs = -1
        integer :: topComm = -1

        !MESH VARIABLES
            !nDim independent
        character (len=15) :: meshType, meshMod;
        integer :: nDim = -1, xNTotal = -1;
        logical :: independent
            !nDim dependent
        integer         , dimension(:), allocatable :: xNStep, procPerDim;
        integer         , dimension(:), allocatable :: neigh, coords;
        double precision, dimension(:), allocatable :: xMax, xMin; !Exact Values of the division
        double precision, dimension(:), allocatable :: xStep;
        double precision, dimension(:), allocatable :: xMaxLoc, xMinLoc; !Rounded Values of the non-overlapping area
        double precision, dimension(:), allocatable :: xMaxBound, xMinBound; !Bounding box of the domain in this proc
        double precision, dimension(:), allocatable :: xMaxGlob, xMinGlob;
        double precision, dimension(:,:), allocatable :: xMaxNeigh, xMinNeigh; !Rounded Values of the overlapping area
        integer         , dimension(:,:), allocatable :: indexNeigh, neighShift
        double precision :: overlap = 0.5; !Absolute value (should change to percentual of corrL)
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
            allocate(MESH_a%xMaxLoc(nDim))
            allocate(MESH_a%xMinLoc(nDim))
            allocate(MESH_a%procPerDim(nDim))
            allocate(MESH_a%coords(nDim))
            allocate(MESH_a%xMaxBound(nDim))
            allocate(MESH_a%xMinBound(nDim))
            allocate(MESH_a%neigh((3**nDim)-1))
            allocate(MESH_a%xMaxNeigh(nDim,(3**nDim)-1))
            allocate(MESH_a%xMinNeigh(nDim,(3**nDim)-1))
            allocate(MESH_a%indexNeigh(2,(3**nDim)-1))
            allocate(MESH_a%neighShift(nDim,(3**nDim)-1))
            !allocate(MESH_a%intShift(nDim,(3**nDim)-1))

            !allocate(MESH_a%overLap(nDim))
            MESH_a%xMax     = -1
            MESH_a%xMin     = -1
            MESH_a%xStep    = -1
            MESH_a%xNStep   = -1
            MESH_a%xMaxGlob = -1
            MESH_a%xMinGlob = -1
            MESH_a%xMaxLoc  = -1
            MESH_a%xMinLoc  = -1
            MESH_a%procPerDim = -1
            MESH_a%xMaxNeigh(:,:) = 0
            MESH_a%xMinNeigh(:,:) = 0
            MESH_a%indexNeigh(:,:) = -1
            MESH_a%neigh(:) = -2 !-1 is already the default when the proc is in the topology border
            MESH_a%neighShift(:,:) = 0

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
            character(len = 20) :: dblFmt, matDblFmt, intFmt, matIntFmt

            dblFmt = "T20,F15.5"
            intFmt = "T20,I15"
            matDblFmt = string_vec_join(["T20,",trim(numb2String(MESH_a%nDim)),"F15.5"])
            matIntFmt = string_vec_join(["T20,",trim(numb2String(MESH_a%nDim)),"I15"])
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
                write(*,*) "|  |topComm  = ", MESH_a%topComm
                write(*,*) "|"
                write(*,*) "|  Input---"
                write(*,*) "|  |nDim     = ", MESH_a%nDim
                write(*,*) "|  |independent = ", MESH_a%independent
                write(*,*) "|  |meshType = ", MESH_a%meshType
                write(*,*) "|  |meshMod  = ", MESH_a%meshMod
                write(*,"(A,("//dblFmt//"))") " |  |xMinGlob = ", MESH_a%xMinGlob
                write(*,"(A,("//dblFmt//"))") " |  |xMaxGlob = ", MESH_a%xMaxGlob
                write(*,"(A,("//dblFmt//"))") " |  |xStep    = ", MESH_a%xStep
                write(*,"(A,F15.5)")" |  |overLap  = ", MESH_a%overLap
                write(*,*) "|"
                write(*,*) "|  |Process--"
                write(*,*) "|  |NStep      = ", MESH_a%xNStep
                write(*,"(A,("//dblFmt//"))") " |  |xMin       = ", MESH_a%xMin
                write(*,"(A,("//dblFmt//"))") " |  |xMax       = ", MESH_a%xMax
                write(*,*) "|  |xNTotal    = ", MESH_a%xNTotal
                write(*,*) "|  |procPerDim = ", MESH_a%procPerDim
                write(*,*) "|  |coords     = ", MESH_a%coords
                !write(*,*) "|  |neigh(1:2nDim)      = ", MESH_a%neigh(1:2*MESH_a%nDim)
                write(*,"(A,("//dblFmt//"))") " |  |xMinLoc    = ", MESH_a%xMinLoc
                write(*,"(A,("//dblFmt//"))") " |  |xMaxLoc    = ", MESH_a%xMaxLoc
                write(*,"(A,("//matDblFmt//"))") " |  |xMinNeigh  = ", MESH_a%xMinNeigh
                write(*,"(A,("//matDblFmt//"))") " |  |xMaxNeigh  = ", MESH_a%xMaxNeigh
                write(*,"(A,("//dblFmt//"))") " |  |xMinBound  = ", MESH_a%xMinBound
                write(*,"(A,("//dblFmt//"))") " |  |xMaxBound  = ", MESH_a%xMaxBound
                write(*,"(A,("//intFmt//"))") " |  |neigh      = ", MESH_a%neigh
                write(*,"(A,("//matIntFmt//"))") " |  |neighShift = ", MESH_a%neighShift(:,:)

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

            if (allocated(MESH_a%xMax))       deallocate(MESH_a%xMax)
            if (allocated(MESH_a%xMin))       deallocate(MESH_a%xMin)
            if (allocated(MESH_a%xStep))      deallocate(MESH_a%xStep)
            if (allocated(MESH_a%xNStep))     deallocate(MESH_a%xNStep)
            if (allocated(MESH_a%xMaxGlob))   deallocate(MESH_a%xMaxGlob)
            if (allocated(MESH_a%xMinGlob))   deallocate(MESH_a%xMinGlob)
            if (allocated(MESH_a%procPerDim)) deallocate(MESH_a%procPerDim)
            if (allocated(MESH_a%coords))     deallocate(MESH_a%coords)
            if (allocated(MESH_a%neigh))      deallocate(MESH_a%neigh)
            if (allocated(MESH_a%xMaxNeigh))  deallocate(MESH_a%xMaxNeigh)
            if (allocated(MESH_a%xMinNeigh))  deallocate(MESH_a%xMinNeigh)
            if (allocated(MESH_a%xMaxLoc))    deallocate(MESH_a%xMaxLoc)
            if (allocated(MESH_a%xMinLoc))    deallocate(MESH_a%xMinLoc)
            if (allocated(MESH_a%indexNeigh)) deallocate(MESH_a%indexNeigh)
            if (allocated(MESH_a%xMaxBound))  deallocate(MESH_a%xMaxBound)
            if (allocated(MESH_a%xMinBound))  deallocate(MESH_a%xMinBound)
            if (allocated(MESH_a%neighShift)) deallocate(MESH_a%neighShift)

            MESH_a%init = .false.

        end subroutine finalize_MESH

end module type_MESH
