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
        integer :: log_ID = -1

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
        integer         , dimension(:), allocatable :: pointsPerCorrL;
        double precision, dimension(:), allocatable :: xMaxLoc, xMinLoc; !Rounded Values of the non-overlapping area
        double precision, dimension(:), allocatable :: xMaxExt, xMinExt; !Bounding box of the domain in this proc
        double precision, dimension(:,:), allocatable :: xMaxNeigh, xMinNeigh; !Rounded Values of the overlapping area
        integer         , dimension(2,1) :: indexLocal
        integer         , dimension(:,:), allocatable :: indexNeigh, neighShift
        double precision, dimension(:), allocatable :: overlap !Size of the overlap (in corrL)
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
            allocate(MESH_a%xMaxExt(nDim))
            allocate(MESH_a%xMinExt(nDim))
            allocate(MESH_a%neigh((3**nDim)-1))
            allocate(MESH_a%xMaxNeigh(nDim,(3**nDim)-1))
            allocate(MESH_a%xMinNeigh(nDim,(3**nDim)-1))
            allocate(MESH_a%indexNeigh(2,(3**nDim)-1))
            allocate(MESH_a%neighShift(nDim,(3**nDim)-1))
            allocate(MESH_a%overLap(nDim))
            allocate(MESH_a%pointsPerCorrL(nDim))

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
            MESH_a%overlap(:) = -1.0D0
            MESH_a%neigh(:) = -2 !-1 is already the default when the proc is in the topology border
            MESH_a%neighShift(:,:) = 0

            MESH_a%init = .true.

        end subroutine init_MESH

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_MESH(MESH_a, name, fmt, forLog_in, unit_in)
            !INPUT
            type(MESH), intent(in) :: MESH_a
            character(len=*), intent(in), optional :: name
            character(len = 20), intent(in), optional :: fmt
            integer, intent(in), optional :: unit_in
            logical, intent(in), optional :: forLog_in
            !LOCAL
            character(len = 20) :: dblFmt, matDblFmt, intFmt, matIntFmt
            integer :: unit
            logical :: forLog
            logical :: active

            active = .true.
            unit = 6 !Screen
            if(present(unit_in)) unit = unit_in

            forLog = .false.
            if(present(forLog_in)) then
#ifdef MAKELOG
                if(forLog_in) unit = MESH_a%log_ID
                forLog = forLog_in
#else
                if(forLog_in) write(*,*) "WARNING!!! Inside show_MESH, forLog_in = .true. but MAKELOG was not defined"
                active = .false.
#endif
            end if

            if(active) then

                if(unit <= 0) then
                    write(*,*) "ERROR!!! Inside show_MESH unit = ", unit
                    stop("")
                end if

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
                    write(unit,"(A,("//dblFmt//"))") " |  |xMinExt  = ", MESH_a%xMinExt
                    write(unit,"(A,("//dblFmt//"))") " |  |xMaxExt  = ", MESH_a%xMaxExt
                    call show_MESHneigh(MESH_a, onlyExisting = .false., forLog = forLog, unit_in = unit)

                else
                    write(unit,*) "|  init     = ", MESH_a%init
                    write(unit,*) "|  MESH has not been initialized----"
                end if
                write(unit,*) "|---------------------------------------------------------------"
                write(unit,*) ""

            end if

        end subroutine show_MESH

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_MESHneigh(MESH_a, name, onlyExisting, forLog, unit_in)
            implicit none
            !INPUT
            type(MESH), intent(in) :: MESH_a
            character(len=*), intent(in), optional :: name
            logical, intent(in) :: onlyExisting
            logical, intent(in), optional :: forLog
            integer, intent(in), optional :: unit_in
            !LOCAL
            character(len = 3) :: nDim, space
            character(len = 50) :: fmtNum, fmtChar, fmtDble
            integer :: i
            integer :: unit
            logical :: active

            active = .true.
            unit = 6 !Screen
            if(present(unit_in)) unit = unit_in

            if(present(forLog)) then
#ifdef MAKELOG
                if(forLog) unit = MESH_a%log_ID
#else
                !if(forLog) write(*,*) "WARNING!!! Inside show_MESHneigh, forLog = .true. but MAKELOG was not defined"
                if(forLog) active = .false.
#endif
            end if

            if(active) then

                if(unit <= 0) then
                    write(*,*) "ERROR!!! Inside show_MESHneigh unit = ", unit
                    stop("")
                end if

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
                        write(unit,fmtNum) MESH_a%neigh(i), "|", MESH_a%neighShift(:,i), "|", &
                            MESH_a%xMinNeigh(:,i), "|", MESH_a%xMaxNeigh(:,i)
                    end do

                else
                    write(unit,*) "|  init     = ", MESH_a%init
                    write(unit,*) "|  MESH has not been initialized----"
                end if

                write(unit,*) ""

            end if

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
            if (allocated(MESH_a%xMaxExt))  deallocate(MESH_a%xMaxExt)
            if (allocated(MESH_a%xMinExt))  deallocate(MESH_a%xMinExt)
            if (allocated(MESH_a%neighShift)) deallocate(MESH_a%neighShift)
            if (allocated(MESH_a%overlap))    deallocate(MESH_a%overlap)
            if (allocated(MESH_a%pointsPerCorrL)) deallocate(MESH_a%pointsPerCorrL)

            MESH_a%init = .false.

        end subroutine finalize_MESH

end module type_MESH
