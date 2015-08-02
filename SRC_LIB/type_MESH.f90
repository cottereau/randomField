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
        double precision, dimension(:), allocatable :: xMaxGlob, xMinGlob;
        double precision, dimension(:), allocatable :: xStep;
        double precision, dimension(:), allocatable :: xMaxLoc, xMinLoc; !Rounded Values of the non-overlapping area
        double precision, dimension(:), allocatable :: xMaxBound, xMinBound; !Bounding box of the domain in this proc
        double precision, dimension(:,:), allocatable :: xMaxNeigh, xMinNeigh; !Rounded Values of the overlapping area
        integer         , dimension(2,1) :: indexLocal
        integer         , dimension(:,:), allocatable :: indexNeigh, neighShift
        double precision :: overlap !Size of the overlap (in corrL)
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
            MESH_a%coords(:)  = -1
            MESH_a%xMaxNeigh(:,:) = 0
            MESH_a%xMinNeigh(:,:) = 0
            MESH_a%indexNeigh(:,:) = -1
            MESH_a%overlap = 0.0D0
            MESH_a%neigh(:) = -2 !-1 is already the default when the proc is in the topology border
            MESH_a%neighShift(:,:) = 0

            MESH_a%init = .true.

        end subroutine init_MESH

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_MESH(MESH_a, name, fmt, unit_in)
            !INPUT
            type(MESH), intent(in) :: MESH_a
            character(len=*), intent(in), optional :: name
            character(len = 20), intent(in), optional :: fmt
            integer, intent(in), optional :: unit_in
            !LOCAL
            character(len = 20) :: dblFmt, matDblFmt, intFmt, matIntFmt
            integer :: unit

            unit = 6 !Screen
            if(present(unit_in)) unit = unit_in

            dblFmt = "T20,F15.5"
            intFmt = "T20,I15"
            matDblFmt = string_vec_join(["T20,",trim(numb2String(MESH_a%nDim)),"F15.5"])
            matIntFmt = string_vec_join(["T20,",trim(numb2String(MESH_a%nDim)),"I15"])
            if(present(fmt)) dblFmt = fmt

            write(unit,*) "MESH------------------------------------------------------------"
            if(present(name)) write(unit,*) "|  ", name

            if(MESH_a%init) then
                write(unit,*) "|  init     = ", MESH_a%init
                write(unit,*) "|"
                write(unit,*) "|  MPI---"
                write(unit,*) "|  |rang     = ", MESH_a%rang
                write(unit,*) "|  |comm     = ", MESH_a%comm
                write(unit,*) "|  |nb_procs = ", MESH_a%nb_procs
                write(unit,*) "|  |topComm  = ", MESH_a%topComm
                write(unit,*) "|"
                write(unit,*) "|  Input---"
                write(unit,*) "|  |nDim     = ", MESH_a%nDim
                write(unit,*) "|  |independent = ", MESH_a%independent
                write(unit,*) "|  |meshType = ", MESH_a%meshType
                write(unit,*) "|  |meshMod  = ", MESH_a%meshMod
                write(unit,"(A,("//dblFmt//"))") " |  |xMinGlob = ", MESH_a%xMinGlob
                write(unit,"(A,("//dblFmt//"))") " |  |xMaxGlob = ", MESH_a%xMaxGlob
                write(unit,"(A,("//dblFmt//"))") " |  |xStep    = ", MESH_a%xStep
                write(unit,"(A,F15.5)")" |  |overLap  = ", MESH_a%overLap
                write(unit,*) "|"
                write(unit,*) "|  |Process--"
                write(unit,*) "|  |xNStep     = ", MESH_a%xNStep
                write(unit,"(A,("//dblFmt//"))") " |  |xMin       = ", MESH_a%xMin
                write(unit,"(A,("//dblFmt//"))") " |  |xMax       = ", MESH_a%xMax
                write(unit,*) "|  |xNTotal    = ", MESH_a%xNTotal
                write(unit,*) "|  |procPerDim = ", MESH_a%procPerDim
                write(unit,*) "|  |coords     = ", MESH_a%coords
                write(unit,"(A,("//dblFmt//"))") " |  |xMinLoc    = ", MESH_a%xMinLoc
                write(unit,"(A,("//dblFmt//"))") " |  |xMaxLoc    = ", MESH_a%xMaxLoc
                write(unit,"(A,("//dblFmt//"))") " |  |xMinBound  = ", MESH_a%xMinBound
                write(unit,"(A,("//dblFmt//"))") " |  |xMaxBound  = ", MESH_a%xMaxBound
                call show_MESHneigh(MESH_a, onlyExisting = .false., unit_in = unit)

            else
                write(unit,*) "|  init     = ", MESH_a%init
                write(unit,*) "|  MESH has not been initialized----"
            end if
            write(unit,*) "|---------------------------------------------------------------"
            write(unit,*) ""

        end subroutine show_MESH

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_MESHneigh(MESH_a, name, onlyExisting, unit_in)
            implicit none
            !INPUT
            type(MESH), intent(in) :: MESH_a
            character(len=*), intent(in), optional :: name
            logical, intent(in) :: onlyExisting
            integer, intent(in), optional :: unit_in
            !LOCAL
            character(len = 3) :: nDim, space
            character(len = 50) :: fmtNum, fmtChar, fmtDble
            integer :: i
            integer :: unit


            unit = 6 !Screen
            if(present(unit_in)) unit = unit_in

            nDim = trim(numb2String(MESH_a%nDim))
            space = trim(numb2String(MESH_a%nDim*6 + 1))

            if(present(name)) write(unit,*) "|  ", name

            fmtNum = "(I10, A1, "//nDim//"I6, A1, "//nDim//"F6.2, A1, "//nDim//"F6.2)"
            fmtDble = "(I10, A1, "//nDim//"F6.2, A1, "//nDim//"F6.2, "//nDim//"F6.2, A1, "//nDim//"F6.2)"
            fmtChar = "(A10, A"//space//", A"//space//", A"//space//")"

            if(MESH_a%init) then

                write(unit,fmtChar) "Neighbour", "|Shift                 ", "|xMin                    ", "|xMax                  "

                do i = 1, size(MESH_a%neigh)
                    if(onlyExisting .and. MESH_a%neigh(i)<0) cycle
                    write(unit,fmtNum) MESH_a%neigh(i), "|", MESH_a%neighShift(:,i), "|", MESH_a%xMinNeigh(:,i), "|", MESH_a%xMaxNeigh(:,i)
                end do

            else
                write(unit,*) "|  init     = ", MESH_a%init
                write(unit,*) "|  MESH has not been initialized----"
            end if

            write(unit,*) ""

        end subroutine show_MESHneigh

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
