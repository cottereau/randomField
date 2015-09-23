module mesh_RF

    !use mpi
    use math_RF
    use write_Log_File
    use type_RF
    use type_MESH

    implicit none

contains
    !--------------²---------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_XPoints(MSH, RDF, xPoints)
        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH
        type(RF)   :: RDF
        double precision, dimension(:, :), allocatable, intent(out), target :: xPoints;

        !LOCAL
        integer :: i, j, counterXPoints
        integer, dimension(MSH%nDim) :: tempXNStep
!        double precision, dimension(MSH%nDim) :: xMinForStep, xMaxForStep



!        xMinForStep = MSH%xMinLoc
!        xMaxForStep = MSH%xMaxLoc
!
!        if(MSH%independent) then
!            do i = 1, size(MSH%xMaxNeigh, 2)
!
!                if(MSH%neigh(i)<0) cycle
!
!                call snap_to_grid(MSH, MSH%xMinNeigh(:,i), MSH%xMaxNeigh(:,i))
!                do j = 1, MSH%nDim
!                    if (MSH%xMinNeigh(j,i) < MSH%xMinExt(j)) MSH%xMinExt(j) = MSH%xMinNeigh(j,i)
!                    if (MSH%xMaxNeigh(j,i) > MSH%xMaxExt(j)) MSH%xMaxExt(j) = MSH%xMaxNeigh(j,i)
!
!                    if (MSH%neighShift(j,i) == 1) then
!                        if (MSH%xMinNeigh(j,i) < xMinForStep(j)) xMinForStep(j) = MSH%xMinNeigh(j,i)
!                        if (MSH%xMaxNeigh(j,i) > xMaxForStep(j)) xMaxForStep(j) = MSH%xMaxNeigh(j,i)
!                    end if
!
!                end do
!            end do
!        end if
!
!        MSH%xNStep = find_xNStep(xMinForStep, xMaxForStep, MSH%xStep)

        MSH%xNStep  = find_xNStep(MSH%xMinExt, MSH%xMaxExt, MSH%xStep)
        MSH%xNTotal = product(MSH%xNStep)
        RDF%xNTotal = MSH%xNTotal

        allocate(xPoints(MSH%nDim, MSH%xNTotal))

        !Internal Points
        counterXPoints = 0;
        tempXNStep = find_xNStep(MSH%xMinLoc, MSH%xMaxLoc, MSH%xStep)

        do i = 1, product(tempXNStep)
            call get_Permutation(i, MSH%xMaxLoc, tempXNStep, xPoints(:,i), MSH%xMinLoc, snapExtremes = .true.);
        end do

        counterXPoints = counterXPoints + product(tempXNStep);
        MSH%indexLocal(1,1) = 1
        MSH%indexLocal(2,1) = counterXPoints

        !Border Points
        if(MSH%independent) then
            do j = 1, size(MSH%xMaxNeigh, 2)

                !Check if the neighbour exist
                if(MSH%neigh(j)<0) cycle

                !Take into account only positive neighbours
                if(minval(MSH%neighShift(:,j)) == -1) cycle

                tempXNStep = find_xNStep(MSH%xMinNeigh(:,j), MSH%xMaxNeigh(:,j), MSH%xStep)
                do i = 1, product(tempXNStep)
                     call get_Permutation(i, MSH%xMaxNeigh(:,j), tempXNStep, &
                                          xPoints(:,counterXPoints + i), MSH%xMinNeigh(:,j), &
                                          snapExtremes = .true.);
                end do
                MSH%indexNeigh(1, j) = counterXPoints + 1
                counterXPoints = counterXPoints + product(tempXNStep);
                MSH%indexNeigh(2,j) = counterXPoints
            end do
        end if

        RDF%xPoints => xPoints

    end subroutine set_XPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_XPoints_globCoords(RDF, MSH)
        implicit none

        !INPUT AND OUTPUT
        type(RF)  , intent(inout) :: RDF
        type(MESH), intent(inout) :: MSH

        RDF%origin = find_xNStep(MSH%xMinGlob, MSH%xMinLoc, MSH%xStep)

    end subroutine get_XPoints_globCoords
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine snap_to_grid(MSH, extInf,extSup)
        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH

        !OUTPUT
        double precision, dimension(:), intent(inout) :: extSup, extInf

        !LOCAL VARIABLES
        integer :: i

        do i = 1, MSH%nDim
            extInf(i) = (MSH%xStep(i) * dble(nint((extInf(i) - MSH%xMinGlob(i))/MSH%xStep(i)))) &
                         + MSH%xMinGlob(i)
            extSup(i) = (MSH%xStep(i) * dble(nint((extSup(i) - MSH%xMinGlob(i))/MSH%xStep(i)))) &
                         + MSH%xMinGlob(i)
        end do

    end subroutine snap_to_grid

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine allocate_xPoints(MSH, RDF, xPoints)
        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH
        type(RF)   :: RDF
        double precision, dimension(:, :), allocatable, intent(out), target :: xPoints;

        !LOCAL VARIABLES

        !write(get_fileId(),*) "-> Allocating xPoints";

        !write(get_fileId(),*) "-> Finding xNStep";
        MSH%xNStep = find_xNStep(MSH%xMinLoc, MSH%xMaxLoc, MSH%xStep)
        MSH%xNTotal = product(MSH%xNStep)
        RDF%xNTotal = MSH%xNTotal

        allocate(xPoints(MSH%nDim, MSH%xNTotal))

        RDF%xPoints => xPoints

    end subroutine allocate_xPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_neighbours (MSH)
        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH

        !LOCAL VARIABLES
        integer :: i, j, code, delta;
        integer, dimension(:), allocatable :: shift

        allocate(shift(MSH%nDim))

        !write(*,*) "set_neighbours"

        !Defining lateral neighbours
        if(MSH%nDim == 1) then
            shift = [-1]
            call find_rank (MSH, shift, 1)
            shift = [1]
            call find_rank (MSH, shift, 2)
        else if(MSH%nDim == 2) then
            shift = [0, -1]
            call find_rank (MSH, shift, 1)
            shift = [0, 1]
            call find_rank (MSH, shift, 2)
            shift = [-1, 0]
            call find_rank (MSH, shift, 3)
            shift = [1, 0]
            call find_rank (MSH, shift, 4)
        else if(MSH%nDim == 3) then
            shift = [0, 0, -1]
            call find_rank (MSH, shift, 1)
            shift = [0, 0, 1]
            call find_rank (MSH, shift, 2)
            shift = [0, -1, 0]
            call find_rank (MSH, shift, 3)
            shift = [0, 1, 0]
            call find_rank (MSH, shift, 4)
            shift = [-1, 0, 0]
            call find_rank (MSH, shift, 5)
            shift = [1, 0, 0]
            call find_rank (MSH, shift, 6)
        end if
!        do i = 0, MSH%nDim - 1
!            !write(*,*) "i = ", i
!            call MPI_CART_SHIFT (MSH%topComm,i,1,MSH%neigh(2*i +1),MSH%neigh(2*(i+1)),code)
!        end do

        !Defining corner neighbours
        if(MSH%nDim == 2) then
            shift = [-1, -1]
            call find_rank (MSH, shift, 5)
            shift = [1, -1]
            call find_rank (MSH, shift, 6)
            shift = [-1, 1]
            call find_rank (MSH, shift, 7)
            shift = [1, 1]
            call find_rank (MSH, shift, 8)
        else if(MSH%nDim == 3) then
            shift = [-1, -1, -1]
            call find_rank (MSH, shift, 19)
            shift = [1, -1, -1]
            call find_rank (MSH, shift, 20)
            shift = [-1, 1, -1]
            call find_rank (MSH, shift, 21)
            shift = [1, 1, -1]
            call find_rank (MSH, shift, 22)

            shift = [-1, -1, 1]
            call find_rank (MSH, shift, 23)
            shift = [1, -1, 1]
            call find_rank (MSH, shift, 24)
            shift = [-1, 1, 1]
            call find_rank (MSH, shift, 25)
            shift = [1, 1, 1]
            call find_rank (MSH, shift, 26)
        end if

        !Defining vertex neighbours
        if(MSH%nDim == 3) then
            shift = [-1, -1, 0]
            call find_rank (MSH, shift, 7)
            shift = [1, -1, 0]
            call find_rank (MSH, shift, 8)
            shift = [-1, 1, 0]
            call find_rank (MSH, shift, 9)
            shift = [1, 1, 0]
            call find_rank (MSH, shift, 10)

            shift = [-1, 0, -1]
            call find_rank (MSH, shift, 11)
            shift = [1, 0, -1]
            call find_rank (MSH, shift, 12)
            shift = [-1, 0, 1]
            call find_rank (MSH, shift, 13)
            shift = [1, 0, 1]
            call find_rank (MSH, shift, 14)

            shift = [0, -1, -1]
            call find_rank (MSH, shift, 15)
            shift = [0, 1, -1]
            call find_rank (MSH, shift, 16)
            shift = [0, -1, 1]
            call find_rank (MSH, shift, 17)
            shift = [0, 1, 1]
            call find_rank (MSH, shift, 18)
        end if

        call set_shift(MSH)

        deallocate(shift)

    end subroutine set_neighbours

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine find_rank (MSH, shift, neighPos)

        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH

        !INPUT
        integer, dimension(:), intent(in) :: shift
        integer, intent(in) :: neighPos

        !LOCAL VARIABLES
        integer :: i, code, neigh;
        integer, dimension(:), allocatable :: pos
        logical :: possible


        allocate(pos(MSH%nDim))
        possible = .true.
        pos = MSH%coords+shift

        !MSH%neighShift(:,neighPos) = shift
        !MSH%intShift(:,neighPos) = shift

        do i = 1, MSH%nDim
            if (pos(i) < 0 .or. pos(i) > MSH%procPerDim(i) - 1) then
                MSH%neigh(neighPos) = -1
                !MSH%neighShift(:,neighPos) = 0
                possible = .false.
                exit
            end if
        end do

        if(possible) then
            call MPI_CART_RANK (MSH%topComm,pos,MSH%neigh(neighPos),code)
        end if

        deallocate(pos)

    end subroutine find_rank

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_shift (MSH)

        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH
        !LOCAL VARIABLES
        integer :: code, neighPos;


        do neighPos = 1, size(MSH%neigh)
            !if(MSH%rang == 0) write (*,*) "neighPos = ", neighPos
            if(MSH%neigh(neighPos) < 0) cycle

            call MPI_CART_COORDS (MSH%topComm, MSH%neigh(neighPos), MSH%nDim, MSH%neighShift(:,neighPos), code)

            MSH%neighShift(:,neighPos) = MSH%neighShift(:,neighPos) - MSH%coords

            !if(MSH%rang == 0) write (*,*) "MSH%neighShift(:,neighPos) = ", MSH%neighShift(:,neighPos)

        end do

    end subroutine set_shift

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_overlap_geometry (MSH, corrL)

        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH

        !INPUT
        double precision, dimension(:) :: corrL

        !LOCAL VARIABLES
        integer :: i, d
        integer :: code, neighPos;
        double precision, parameter :: notPresent = -1.0D0
        double precision, dimension(MSH%nDim) :: tempExtreme

        call wLog("-> Redimensioning for Overlap")

        !Redimensioning the internal part
        do neighPos = 1, 2*MSH%nDim

            if(MSH%neigh(neighPos) < 0) cycle

            where(MSH%neighShift(:,neighPos) < 0)
                MSH%xMinLoc = MSH%xMin + MSH%overlap*corrL/2.0D0 + MSH%xStep
            elsewhere(MSH%neighShift(:,neighPos) > 0)
                MSH%xMaxLoc = MSH%xMax - MSH%overlap*corrL/2.0D0 - MSH%xStep
            end where

        end do

        call wLog("MSH%xMinLoc")
        call wLog(MSH%xMinLoc)
        call wLog("MSH%xMaxLoc")
        call wLog(MSH%xMaxLoc)

        !Dimensioning overlapping area
        do neighPos = 1, size(MSH%neigh)
            if(MSH%neigh(neighPos) < 0) cycle

            !write(get_fileId(),*)

            where(MSH%neighShift(:,neighPos) > 0)
                MSH%xMaxNeigh(:,neighPos) = MSH%xMax + MSH%overlap*corrL/2
                MSH%xMinNeigh(:,neighPos) = MSH%xMax - MSH%overlap*corrL/2
            elsewhere(MSH%neighShift(:,neighPos) < 0)
                MSH%xMaxNeigh(:,neighPos) = MSH%xMin + MSH%overlap*corrL/2
                MSH%xMinNeigh(:,neighPos) = MSH%xMin - MSH%overlap*corrL/2
            elsewhere
                MSH%xMaxNeigh(:,neighPos) = MSH%xMaxLoc
                MSH%xMinNeigh(:,neighPos) = MSH%xMinLoc
            end where

            !Redimensioning Bounding Box
            where(MSH%xMaxNeigh(:,neighPos) > MSH%xMaxExt(:)) MSH%xMaxExt(:) = MSH%xMaxNeigh(:,neighPos)
            where(MSH%xMinNeigh(:,neighPos) > MSH%xMinExt(:)) MSH%xMinExt(:) = MSH%xMinNeigh(:,neighPos)
        end do

        call show_MESHneigh(MSH, " ", onlyExisting = .false., forLog = .true.)

    end subroutine get_overlap_geometry

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_Global_Extremes_Mesh(xPoints, xMinGlob, xMaxGlob, communicator)
        !Find the extremes of the mesh
        implicit none

        !INPUT
        double precision, dimension(1:, 1:), intent(in) :: xPoints;
        integer         , optional         , intent(in) :: communicator
        !OUTPUT
        double precision, dimension(1:), intent(out) :: xMinGlob, xMaxGlob;
        !LOCAL
        integer :: code, nDim, i;
        double precision, dimension(:), allocatable :: xMinLoc, xMaxLoc;
        integer :: effectComm

        if(present(communicator)) then
            effectComm = communicator
        else
            effectComm = MPI_COMM_WORLD
        end if

        nDim = size(xPoints, 1)


        allocate(xMinLoc(nDim))
        allocate(xMaxLoc(nDim))

        xMinLoc = minval(xPoints, 2)
        xMaxLoc = maxval(xPoints, 2)

        do i = 1, nDim
            call MPI_ALLREDUCE (xMinLoc(i), xMinGlob(i), 1, MPI_DOUBLE_PRECISION, &
                MPI_MIN, effectComm, code)
            call MPI_ALLREDUCE (xMaxLoc(i), xMaxGlob(i), 1, MPI_DOUBLE_PRECISION, &
                MPI_MAX, effectComm, code)
        end do

        deallocate(xMinLoc)
        deallocate(xMaxLoc)

    end subroutine get_Global_Extremes_Mesh

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function find_xNStep(xMin, xMax, xStep) result (xNStep)
        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: xMin, xMax, xStep;

        !OUTPUT
        integer, dimension(:), allocatable :: xNStep
        integer :: i

        allocate(xNStep(size(xStep)))

        xNStep = 1 + nint((xMax-xMin)/(xStep));

        do i = 1, size(xStep)
            if(xNStep(i) < 1) then
                write(*,*) "ERROR!!! Inside find_xNStep, xMin is greater than xMax"
                !write(get_fileId(),*) "ERROR!!! Inside find_xNStep, xMin is greater than xMax"
                !write(get_fileId(),*) " xMax = ", xMax
                !write(get_fileId(),*) " xMin = ", xMin
                !write(get_fileId(),*) " xStep = ", xStep
                !write(get_fileId(),*) " xNStep = ", xNStep
                stop(" ")
            end if
        end do

    end function find_xNStep

    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    subroutine define_topography(RDF, MSH, coordList)

        implicit none
        !INPUT OUTPUT
        type(RF)   :: RDF
        type(MESH) :: MSH
        !INPUT
        double precision, dimension(:,:), target, intent(in), optional :: coordList
        !LOCAL
        logical, dimension(MSH%nDim) :: periods
        integer :: code

        periods(:) = .false.


        if (MSH%meshMod == "unv") then
            RDF%xPoints => coordList
            call wLog("-> defining_UNV_extremes")
            call get_Global_Extremes_Mesh(RDF%xPoints, MSH%xMinGlob, MSH%xMaxGlob, RDF%comm)
            RDF%xNTotal = MSH%xNTotal
            RDF%xMinExt = MSH%xMinGlob
            MSH%xMinExt = MSH%xMinGlob
            RDF%xMaxExt = MSH%xMaxGlob
            MSH%xMaxExt = MSH%xMaxGlob
            MSH%xMin    = minval(RDF%xPoints, 2)
            MSH%xMax    = maxval(RDF%xPoints, 2)
            MSH%xNTotal = size(RDF%xPoints, 2)
            RDF%xNTotal = MSH%xNTotal

        else
            call wLog("-> set_procPerDim")
            call set_procPerDim (MSH%nb_procs, MSH%nDim, MSH%procPerDim)
            call wLog("-> MPI_CART_CREATE")
            call MPI_CART_CREATE (MSH%comm, MSH%nDim, MSH%procPerDim, periods, .false., MSH%topComm, code)
            call wLog("-> MPI_CART_COORDS")
            call MPI_CART_COORDS (MSH%topComm, MSH%rang, MSH%nDim, MSH%coords, code)
            call wLog("-> define_generation_geometry")
            call define_generation_geometry (MSH, RDF)

            if(RDF%independent) then
                call wLog("-> set_neighbours")
                call set_neighbours (MSH)
                call wLog("-> get_overlap_geometry")
                call get_overlap_geometry (MSH, RDF%corrL)
            end if

            !write(get_fileId(),*) "-> Getting Global Matrix Reference"
            call get_XPoints_globCoords(RDF, MSH)
            !write(get_fileId(),*) "     RDF%origin = ", RDF%origin
            !write(get_fileId(),*) " "

        end if

    end subroutine define_topography

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine define_generation_geometry (MSH, RDF)

        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH
        type(RF), intent(in) :: RDF

        !LOCAL
        integer :: i
        double precision, dimension(MSH%nDim) :: delta, half, localSpace, localBase, xLocReq
        double precision, dimension(MSH%nDim) :: locA, locO, start
        !Defining MSH%xStep
        !write(get_fileId(), *) " 	IN  RDF%corrL          = ", RDF%corrL
        !write(get_fileId(), *) "        IN  MSH%pointsPerCorrL = ", MSH%pointsPerCorrL
        MSH%xStep = RDF%corrL/dble(MSH%pointsPerCorrL)
        !write(get_fileId(), *) "        OUT MSH%xStep          = ", MSH%xStep
        !write(get_fileId(), *) " "

        !write(get_fileId(), *) "        IN   MSH%overlap =  ", MSH%overlap

        !Rounding overlap
        if(MSH%independent) then
            MSH%overlap = ceiling(MSH%overlap*RDF%corrL/(2*MSH%xStep)) * 2*MSH%xStep/RDF%corrL
            !!write(get_fileId(), *) "   Rounded MSH%overlap = ", MSH%overlap
        else
            MSH%overlap = 0
        end if

        !write(get_fileId(), *) "        OUT  MSH%overlap =  ", MSH%overlap
        !write(get_fileId(), *) " "

        !Local areas
        localSpace = (MSH%xMaxGlob - MSH%xMinGlob) - MSH%overlap*RDF%corrL*dble(MSH%procPerDim-1)

        !!write(get_fileId(), *) " localSpace  BEFORE= ", localSpace
        if(MSH%independent) then
            localBase  = 4*MSH%xStep*dble(MSH%procPerDim)
        else
            localBase  = MSH%xStep*dble(MSH%procPerDim)
        end if

        !!write(get_fileId(), *) " localBase = ", localBase

        where(localSpace < localBase)
            localSpace = localBase
        elsewhere
             localSpace = dble(ceiling(localSpace/localBase)) * localBase
        end where

        !!write(get_fileId(), *) " localSpace AFTER = ", localSpace

        !Redefining global extremes
!        call wLog("  MSH%procPerDim= ")
!        call wLog(MSH%procPerDim)
!        call wLog("  RDF%corrL= ")
!        call wLog(RDF%corrL)
!        call wLog("  MSH%pointsPerCorrL= ")
!        call wLog(MSH%pointsPerCorrL)
!        call wLog("  MSH%xStep= ")
!        call wLog(MSH%xStep)
!        call wLog("  MSH%independent= ")
!        call wLog(MSH%independent)
!        call wLog("  MSH%overlap= ")
!        call wLog(MSH%overlap)
!        call wLog("  procPerDim= ")
!        call wLog(MSH%procPerDim)
        call wLog(" ")
        call wLog("        IN  MSH%xMinGlob = ")
        call wLog(MSH%xMinGlob)
        call wLog("        IN  MSH%xMaxGlob = ")
        call wLog(MSH%xMaxGlob)
        call wLog("        IN  delta        = ")
        call wLog(MSH%xMaxGlob - MSH%xMinGlob)

        half  = (MSH%xMaxGlob + MSH%xMinGlob)/2.0D0
        delta = localSpace + MSH%overlap*RDF%corrL*dble(MSH%procPerDim-1)
        MSH%xMinGlob = half - delta/2.0D0
        MSH%xMaxGlob = half + delta/2.0D0

        call wLog("        OUT MSH%xMinGlob = ")
        call wLog(MSH%xMinGlob)
        call wLog("        OUT MSH%xMaxGlob = ")
        call wLog(MSH%xMaxGlob)
        call wLog("        OUT delta        = ")
        call wLog(MSH%xMaxGlob - MSH%xMinGlob)
        call wLog(" ")

        !Setting Local Extremes
        locO =  MSH%overlap*RDF%corrL/2.0D0
        locA =  (delta - (MSH%overlap * RDF%corrL * dble(MSH%procPerDim-1)))/dble(MSH%procPerDim)

!        call wLog(" locO = ")
!        call wLog(locO)
!        call wLog(" locA = ")
!        call wLog(locA)

        delta = locA

        where(MSH%procPerDim /= 1) delta = delta + locO
        where(MSH%coords /= 0 .and. MSH%coords /= MSH%procPerDim-1)
            delta = delta + locO
        end where

!        call wLog(" delta = ")
!        call wLog(delta)

        MSH%xMin = MSH%xMinGlob
        where(MSH%coords > 0) MSH%xMin = MSH%xMin + (locA + locO) + (locA + 2.0D0*locO)*(dble(MSH%coords-1))

        MSH%xMax = MSH%xMin + delta

        if(.not. MSH%independent) then
            where(MSH%coords /= MSH%procPerDim-1) MSH%xMax = MSH%xMax-MSH%xStep
        end if

        MSH%xMinLoc   = MSH%xMin
        MSH%xMaxLoc   = MSH%xMax

        !Snaping points to the grid
        call snap_to_grid(MSH, MSH%xMinLoc, MSH%xMaxLoc)

        !Bounding Box
        MSH%xMaxExt = MSH%xMaxLoc
        MSH%xMinExt = MSH%xMinLoc

!        call wLog("        OUT MSH%xMin    = ")
!        call wLog(MSH%xMin)
!        call wLog("        OUT MSH%xMax    = ")
!        call wLog(MSH%xMax)
        call wLog("        OUT MSH%xMinLoc = ")
        call wLog(MSH%xMinLoc)
        call wLog("        OUT MSH%xMaxLoc = ")
        call wLog(MSH%xMaxLoc)

    end subroutine define_generation_geometry

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine set_Local_Extremes_From_Coords (MSH)
!        implicit none
!
!        !INPUT AND OUTPUT
!        type(MESH) :: MSH
!
!        !LOCAL
!        double precision, dimension(MSH%nDim) :: procDelta
!        integer :: i
!
!        procDelta = (MSH%xMaxGlob - MSH%xMinGlob)/MSH%procPerDim
!
!        MSH%xMin = procDelta*MSH%coords + MSH%xMinGlob
!        MSH%xMax = MSH%xMin + procDelta
!
!        if(.not. MSH%independent) then
!            where(MSH%coords /= MSH%procPerDim-1) MSH%xMax = MSH%xMax-MSH%xStep
!        end if
!
!        MSH%xMinLoc = MSH%xMin
!        MSH%xMaxLoc = MSH%xMax
!
!    end subroutine set_Local_Extremes_From_Coords
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine get_Permutation_from_Mat(pos, matrix, nDim, pVec);
!
!        implicit none
!
!        !INPUT
!        integer                        , intent(in)           :: pos;
!        double precision, dimension(1:, 1:), intent(in)           :: matrix;
!        integer, intent(in) :: nDim
!
!        !OUTPUT
!        double precision, dimension(1:), intent(out) :: pVec;
!        !LOCAL VARIABLES
!        integer :: i, j;
!        integer :: seedStep;
!        integer, dimension(:), allocatable :: nStep;
!
!        allocate(nStep(nDim))
!
!        nStep = size(matrix,1)
!
!        do j = 1, nDim
!            if (j /= nDim) seedStep = product(nStep(j+1:));
!            if (j == nDim) seedStep = 1;
!            i = cyclicMod(int((pos-0.9)/seedStep)+1, nStep(j))
!            pVec(j) = matrix(i, j);
!        end do
!
!        deallocate(nStep)
!
!    end subroutine get_Permutation_from_Mat

end module mesh_RF
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
