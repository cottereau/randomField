module mesh_RF

    !use mpi
    use math_RF
    use write_Log_File
    use type_RF
    use type_MESH

    implicit none

contains
    !-----------------------------------------------------------------------------------------------
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
        double precision, dimension(MSH%nDim) :: xMinForStep, xMaxForStep

        !Snaping points to the grid and discover the bounding box
        call snap_to_grid(MSH, MSH%xMinLoc, MSH%xMaxLoc)
        MSH%xMaxBound = MSH%xMaxLoc;
        MSH%xMinBound = MSH%xMinLoc;
        xMinForStep = MSH%xMinLoc
        xMaxForStep = MSH%xMaxLoc

        if(MSH%independent) then
            do i = 1, size(MSH%xMaxNeigh, 2)

                if(MSH%neigh(i)<0) cycle

                call snap_to_grid(MSH, MSH%xMinNeigh(:,i), MSH%xMaxNeigh(:,i))
                do j = 1, MSH%nDim
                    if (MSH%xMinNeigh(j,i) < MSH%xMinBound(j)) MSH%xMinBound(j) = MSH%xMinNeigh(j,i)
                    if (MSH%xMaxNeigh(j,i) > MSH%xMaxBound(j)) MSH%xMaxBound(j) = MSH%xMaxNeigh(j,i)

                    if (MSH%neighShift(j,i) == 1) then
                        if (MSH%xMinNeigh(j,i) < xMinForStep(j)) xMinForStep(j) = MSH%xMinNeigh(j,i)
                        if (MSH%xMaxNeigh(j,i) > xMaxForStep(j)) xMaxForStep(j) = MSH%xMaxNeigh(j,i)
                    end if

                end do
            end do
        end if

        MSH%xNStep = find_xNStep(xMinForStep, xMaxForStep, MSH%xStep)
        MSH%xNTotal = product(MSH%xNStep)
        RDF%xNTotal = MSH%xNTotal

        !write(get_fileId(),*) "      xMinForStep    = ", xMinForStep
        !write(get_fileId(),*) "      xMaxForStep    = ", xMaxForStep
        !write(get_fileId(),*) "      MSH%xMinBound  = ", MSH%xMinBound
        !write(get_fileId(),*) "      MSH%xMaxBound  = ", MSH%xMaxBound
        !write(get_fileId(),*) "      MSH%xNStep     = ", MSH%xNStep

        allocate(xPoints(MSH%nDim, MSH%xNTotal))

        !Internal Points
        counterXPoints = 0;
        tempXNStep = find_xNStep(MSH%xMinLoc, MSH%xMaxLoc, MSH%xStep)
        !write(get_fileId(),*) "   Internal Points"
        !write(get_fileId(),*) "      MSH%xMinLoc = ", MSH%xMinLoc
        !write(get_fileId(),*) "      MSH%xMaxLoc = ", MSH%xMaxLoc
        !write(get_fileId(),*) "      MSH%xStep   = ", MSH%xStep
        !write(get_fileId(),*) "      tempXNStep  = ", tempXNStep

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
                     call get_Permutation(i, MSH%xMaxNeigh(:,j), tempXNStep, xPoints(:,counterXPoints + i), MSH%xMinNeigh(:,j), snapExtremes = .true.);
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

        write(get_fileId(),*) "-> Allocating xPoints";

        write(get_fileId(),*) "-> Finding xNStep";
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
    subroutine set_XPoints_perCorrL(xMin, xMax, pointsPerCorrL, corrL, xPoints, xMinGlob, xMaxGlob, cutExtremes)
        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: corrL;
        double precision, dimension(:), intent(inout) :: xMin, xMax;
        integer, intent(in) :: pointsPerCorrL
        double precision, dimension(:), intent(in), optional :: xMinGlob, xMaxGlob
        logical, optional, intent(in) :: cutExtremes

        !OUTPUT
        double precision, dimension(:,:), allocatable, intent(OUT) :: xPoints;


        !LOCAL VARIABLES
        integer :: nDim, i, xNTotal;
        double precision, dimension(:), allocatable :: xStep

        nDim    = size(xMax)
        allocate(xStep(nDim))

        xStep = corrL/dble(pointsPerCorrL-1)
        !call recalculate_xStep(xMin, xMax, xStep)

        write(*,*) "xStep = ", xStep
        write(*,*) "xMax  = ", xMax
        write(*,*) "xMin  = ", xMin
        write(*,*) "corrL = ", corrL
        write(*,*) "pointsPerCorrL = ", pointsPerCorrL

!        if(present(cutExtremes)) then
!            if(cutExtremes) call cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
!        end if

        !call set_XPoints(xMin, xMax, xStep, xPoints)

!        if(present(cutExtremes)) then
!            if(cutExtremes) call restoreBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
!        end if

        deallocate(xStep)

    end subroutine set_XPoints_perCorrL

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine set_Local_Extremes_Mesh (xMin, xMax, xMinGlob, xMaxGlob, rang, nb_procs)
!        !Find the Boundaries for the box in each processor when using automatic mesh
!        implicit none
!
!        !INPUT
!        integer                        , intent(in) :: rang, nb_procs;
!        double precision, dimension(1:), intent(in) :: xMaxGlob;
!        double precision, dimension(1:), intent(in) :: xMinGlob;
!        !OUTPUT
!        double precision, dimension(1:), intent(out) :: xMax;
!        double precision, dimension(1:), intent(out) :: xMin;
!
!        !LOCAL VARIABLES
!        integer :: i, j, testRang = 0;
!        integer :: baseStep
!        integer :: seedStep, nDim, basicStep;
!        integer, allocatable, dimension(:) :: bStepVec
!        double precision, dimension(:), allocatable :: xProcDelta;
!        double precision :: procIterations
!
!        nDim = size(xMin);
!        baseStep = nint(dble(nb_procs)**(1.0d0/nDim))
!        write(get_fileId(),*) "Setting Local Extremes"
!        allocate (xProcDelta(nDim))
!
!
!        if (nb_procs == baseStep**nDim) then
!            write(get_fileId(),*) "Exact Division"
!            basicStep  = nint(dble(nb_procs)**(1.0d0/nDim))
!            xProcDelta = (xMaxGlob-xMinGlob)/basicStep
!
!            !        if(rang == testRang) write(*,*) "nDim = ", nDim
!            !        if(rang == testRang) write(*,*) "basicStep = ", basicStep
!            !        if(rang == testRang) write(*,*) "xProcDelta = ", xProcDelta
!            !        if(rang == testRang) write(*,*) "nb_procs = ", nb_procs
!            !        if(rang == testRang) write(*,*) "dble(nb_procs)**(1/nDim) = ", dble(nb_procs)**(1/nDim)
!            !        if(rang == testRang) write(*,*) "nint(dble(nb_procs)**(1/nDim)) = ", nint(dble(nb_procs)**(1/nDim))
!
!            do j = 1, nDim
!                seedStep = basicStep**(nDim-j);
!                i = cyclicMod(int(rang/seedStep) + 1, basicStep)
!                !if(rang == testRang) write(*,*) "i = ", i, " seedStep = ", seedStep, "seedStep", seedStep
!                xMin(j) = (dble(i-1))*xProcDelta(j)+xMinGlob(j);
!                xMax(j) = xMin(j) + xProcDelta(j)
!                !if(rang == testRang) write(*,*) "AFTER xMin ", j, " = ", xMin(j)
!            end do
!
!
!
!            !write(*,*) "Proc = ", rang, "xMin = ", xMin!, "xMax = ", xMax;
!
!        !write(*,*) "Proc = ", rang, "xMin = ", xMin!, "xMax = ", xMax;
!        !write(*,*) "Proc = ", rang, "xMin = ", xMin!, "xMax = ", xMax;
!
!        else if (isPowerOf(dble(nb_procs), 2)) then
!            write(get_fileId(),*) "Power of two"
!            allocate (bStepVec(nDim))
!
!            !Defining the basic Step for each dimension
!            bStepVec(:) = 1
!            if(nb_procs /= 1) then
!                procIterations = log(dble(nb_procs))/log(2.0D0)
!                do j = 1, nint(procIterations)
!                    i = cyclicMod(j, nDim)
!                    bStepVec(i) = bStepVec(i)*2
!                    !write(*,*) "i = ", i
!                    !write(*,*) "bStepVec = ", bStepVec
!                end do
!            end if
!            !Defining coordinates in each proc
!            xProcDelta = (xMaxGlob-xMinGlob)/bStepVec
!
!            write(get_fileId(),*) "xMaxGlob = ", xMaxGlob
!            write(get_fileId(),*) "xMinGlob = ", xMinGlob
!            write(get_fileId(),*) "nb_procs = ", nb_procs
!            write(get_fileId(),*) "bStepVec = ", bStepVec
!            write(get_fileId(),*) "xProcDelta = ", xProcDelta
!
!            do j = 1, nDim
!                seedStep = product(bStepVec(j+1:));
!                if (j == nDim) seedStep = 1;
!                i = cyclicMod(int(rang/seedStep) + 1, bStepVec(j))
!                !if(rang == testRang) write(*,*) "i = ", i, " seedStep = ", seedStep, "seedStep", seedStep
!                xMin(j) = (dble(i-1))*xProcDelta(j)+xMinGlob(j);
!                xMax(j) = xMin(j) + xProcDelta(j)
!                !if(rang == testRang) write(*,*) "AFTER xMin ", j, " = ", xMin(j)
!            end do
!            deallocate (bStepVec)
!        else
!            stop "ERROR, no mesh division algorithm for this number of procs"
!        end if
!
!        do i = 1, nDim
!            write(get_fileId(),*) "Dim ", i
!            write(get_fileId(),fmt="(2A30)") "xMin", "xMax"
!            write(get_fileId(),fmt="(2F30.15)") xMin(i), xMax(i)
!        end do
!
!        deallocate (xProcDelta)
!
!    end subroutine set_Local_Extremes_Mesh

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_Local_Extremes_From_Coords (MSH)
        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH

        !LOCAL
        double precision, dimension(MSH%nDim) :: procDelta
        integer :: i

        procDelta = (MSH%xMaxGlob - MSH%xMinGlob)/MSH%procPerDim

        MSH%xMin = procDelta*MSH%coords + MSH%xMinGlob
        MSH%xMax = MSH%xMin + procDelta

        if(.not. MSH%independent) then
            where(MSH%coords /= MSH%procPerDim-1) MSH%xMax = MSH%xMax-MSH%xStep
        end if

        MSH%xMinLoc = MSH%xMin
        MSH%xMaxLoc = MSH%xMax

    end subroutine set_Local_Extremes_From_Coords

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_procPerDim (MSH)
        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH

        !LOCAL VARIABLES
        integer :: i, j;
        double  precision :: procRootDim, logProc2;

        procRootDim = dble(MSH%nb_procs)**(1/dble(MSH%nDim))
        logProc2   = log(dble(MSH%nb_procs))/log(2.0D0)

        if (areEqual(procRootDim, dble(nint(procRootDim)))) then
            write(get_fileId(),*) "    Exact Division"
            !write(*,*) "Exact Division"
            MSH%procPerDim(:) = nint(dble(MSH%nb_procs)**(1.0d0/MSH%nDim))
        else if(areEqual(logProc2, dble(nint(logProc2)))) then
            write(get_fileId(),*) "    Power of two"
            !write(*,*) "Power of two"

            MSH%procPerDim(:) = 1
            if(MSH%nb_procs /= 1) then
                do j = 1, nint(logProc2)
                    i = cyclicMod(j, MSH%nDim)
                    MSH%procPerDim(i) = MSH%procPerDim(i)*2
                end do
            end if
        else
            stop "ERROR, no mesh division algorithm for this number of procs"
        end if

    end subroutine set_procPerDim

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
    subroutine define_generation_geometry (MSH, RDF)

        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH
        type(RF), intent(in) :: RDF

        !LOCAL
        integer :: i
        double precision, dimension(MSH%nDim) :: delta, half, localSpace, localBase
        double precision :: constBase
        !Defining MSH%xStep
        MSH%xStep = RDF%corrL/dble(MSH%pointsPerCorrL)
        !write(get_fileId(), *) " 	MSH%xStep = ", MSH%xStep

        !Rounding overlap        
        if(MSH%independent) then
            MSH%overlap = ceiling(MSH%overlap*RDF%corrL/(2*MSH%xStep)) * 2*MSH%xStep/RDF%corrL
            !write(get_fileId(), *) " 	Rounded MSH%overlap = ", MSH%overlap
        else
            MSH%overlap = 0
        end if

        !Local areas
        localSpace = (MSH%xMaxGlob - MSH%xMinGlob) - MSH%overlap*dble(MSH%procPerDim-1)
        localBase  = MSH%xStep*MSH%procPerDim
        constBase = 1.0D0
        if(MSH%independent) constBase = 2.0D0
        
        where(localSpace < constBase*localBase) 
            localSpace = constBase*localBase
        elsewhere
            localSpace = dble(ceiling(localSpace/localBase)) * localBase
        end where
        !write(get_fileId(), *) " 	localSpace = ", localSpace
        
        !Redefining global extremes
        half  = (MSH%xMaxGlob + MSH%xMinGlob)/2.0D0
        !Defining rounded delta between max and min
        delta = localSpace + MSH%overlap*(dble(MSH%procPerDim)-1) 
        !write(get_fileId(), *) " delta  = ", MSH%xMaxGlob - MSH%xMinGlob
        MSH%xMinGlob = half - delta/2.0D0
        MSH%xMaxGlob = half + delta/2.0D0

        !write(get_fileId(), *) " 	MSH%xMinGlob = ", MSH%xMinGlob
        !write(get_fileId(), *) " 	MSH%xMaxGlob = ", MSH%xMaxGlob
        write(get_fileId(), *) " "
        write(get_fileId(), *) "       nProcsPerDim      "
        write(get_fileId(), *) "   ", MSH%procPerDim
        write(get_fileId(), *) " 	Area For Overlap  "
        write(get_fileId(), *) "   ", MSH%overlap * dble(MSH%procPerDim-1)
        write(get_fileId(), *) " 	Area For Local    "
        write(get_fileId(), *) "   ", delta - (MSH%overlap * dble(MSH%procPerDim-1))
        write(get_fileId(), *) " "

    end subroutine define_generation_geometry

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

        !Redimensioning the internal part
        do neighPos = 1, 2*MSH%nDim

            if(MSH%neigh(neighPos) < 0) cycle

            where(MSH%neighShift(:,neighPos) < 0)
                MSH%xMinLoc = MSH%xMin + MSH%overlap*corrL/2.0D0 + MSH%xStep
            elsewhere(MSH%neighShift(:,neighPos) > 0)
                MSH%xMaxLoc = MSH%xMax - MSH%overlap*corrL/2.0D0 - MSH%xStep
            end where

!            do i = 1, MSH%nDim
!                if(MSH%neighShift(i,neighPos) < 0) then
!                    MSH%xMinLoc(i) = MSH%xMin(i) + MSH%overlap(i)*corrL(i)/2 + MSH%xStep(i)
!                else if (MSH%neighShift(i,neighPos) > 0) then
!                    MSH%xMaxLoc(i) = MSH%xMax(i) - MSH%overlap(i)*corrL(i)/2 - MSH%xStep(i)
!                end if
!
!                !Checking if the points are in the range
!                if(MSH%xMinLoc(i) < MSH%xMinGlob(i)) then
!                    write(get_fileId(),*)"WARNING xMinLoc exceeded global dimensions (consider changing the overlap size or the global size)"
!                    write(get_fileId(),*)"neigh = ",neighPos," dim = ", i
!                    write(get_fileId(),*)"MSH%xMinLoc(i) = ", MSH%xMinLoc(i)
!                    write(get_fileId(),*)"MSH%xMinGlob(i)    = ", MSH%xMinGlob(i)
!                    MSH%xMinLoc(i) = MSH%xMinGlob(i)
!                end if
!                if(MSH%xMaxLoc(i) > MSH%xMaxGlob(i)) then
!                    write(get_fileId(),*)"WARNING xMaxLoc exceeded global dimensions (consider changing the overlap size or the global size)"
!                    write(get_fileId(),*)"neigh = ",neighPos," dim = ", i
!                    write(get_fileId(),*)"MSH%xMaxLoc(i) = ", MSH%xMaxLoc(i)
!                    write(get_fileId(),*)"MSH%xMaxGlob(i)    = ", MSH%xMaxGlob(i)
!                    MSH%xMaxLoc(i) = MSH%xMaxGlob(i)
!                end if
!
!            end do

        end do

        !Dimensioning overlapping area
        do neighPos = 1, size(MSH%neigh)
            if(MSH%neigh(neighPos) < 0) cycle

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

!            do i = 1, MSH%nDim
!                if(MSH%neighShift(i,neighPos) > 0) then
!                    MSH%xMaxNeigh(i,neighPos) = MSH%xMax(i) + MSH%overlap(i)*corrL(i)/2
!                    MSH%xMinNeigh(i,neighPos) = MSH%xMax(i) - MSH%overlap(i)*corrL(i)/2
!                else if (MSH%neighShift(i,neighPos) < 0) then
!                    MSH%xMaxNeigh(i,neighPos) = MSH%xMin(i) + MSH%overlap(i)*corrL(i)/2
!                    MSH%xMinNeigh(i,neighPos) = MSH%xMin(i) - MSH%overlap(i)*corrL(i)/2
!                else
!                    MSH%xMaxNeigh(i,neighPos) = MSH%xMaxLoc(i)
!                    MSH%xMinNeigh(i,neighPos) = MSH%xMinLoc(i)
!
!                end if
!
!                !Checking if the points are in the range
!                if(MSH%xMinNeigh(i,neighPos) < MSH%xMinGlob(i)) then
!                    write(get_fileId(),*)"WARNING!! MSH%xMinNeigh exceeded global dimensions (consider changing the overlap size or the global size)"
!                    write(get_fileId(),*)"neigh = ",neighPos," dim = ", i
!                    write(get_fileId(),*)"MSH%xMinNeigh(i,neighPos) = ", MSH%xMinNeigh(i,neighPos)
!                    write(get_fileId(),*)"MSH%xMinGlob(i)    = ", MSH%xMinGlob(i)
!                    MSH%xMinNeigh(i,neighPos) = MSH%xMinGlob(i)
!                end if
!                if(MSH%xMaxNeigh(i,neighPos) > MSH%xMaxGlob(i)) then
!                    write(get_fileId(),*)"WARNING!! MSH%xMaxNeigh exceeded global dimensions (consider changing the overlap size or the global size)"
!                    write(get_fileId(),*)"neigh = ",neighPos," dim = ", i
!                    write(get_fileId(),*)"MSH%xMaxNeigh(i,neighPos) = ", MSH%xMaxNeigh(i,neighPos)
!                    write(get_fileId(),*)"MSH%xMaxGlob(i)    = ", MSH%xMaxGlob(i)
!                    MSH%xMaxNeigh(i,neighPos) = MSH%xMaxGlob(i)
!                end if
!            end do
        end do

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

        !write(*,*) ">>>>>>>>> Communicating Extremes (unstructured) "

        if(present(communicator)) then
            effectComm = communicator
        else
            effectComm = MPI_COMM_WORLD
        end if

        nDim = size(xPoints, 1)

        !write(*,*) "nDim = ", nDim

        allocate(xMinLoc(nDim))
        allocate(xMaxLoc(nDim))

        xMinLoc = minval(xPoints, 2)
        xMaxLoc = maxval(xPoints, 2)

        !call dispCarvalhol(xMinLoc, "xMinLoc")
        !call dispCarvalhol(xMaxLoc, "xMaxLoc")

        do i = 1, nDim
            !write(*,*) "i = ", i
            call MPI_ALLREDUCE (xMinLoc(i), xMinGlob(i), 1, MPI_DOUBLE_PRECISION, &
                MPI_MIN, effectComm, code)
            call MPI_ALLREDUCE (xMaxLoc(i), xMaxGlob(i), 1, MPI_DOUBLE_PRECISION, &
                MPI_MAX, effectComm, code)
        end do

        !call dispCarvalhol(xMinGlob, "xMinGlob")
        !call dispCarvalhol(xMaxGlob, "xMaxGlob")

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

!        write(*,*) "HERE find_xNStep"
!        write(*,*) " (xMax-xMin) = ", (xMax-xMin)
!        write(*,*) " (xStep) = ", (xStep)
!        write(*,*) " nint((xMax-xMin)/(xStep)) = ", nint((xMax-xMin)/(xStep))
!        write(*,*) " 1 + nint((xMax-xMin)/(xStep)) = ", 1 + nint((xMax-xMin)/(xStep))

        xNStep = 1 + nint((xMax-xMin)/(xStep));

        do i = 1, size(xStep)
            if(xNStep(i) < 1) stop "ERROR!!! Inside find_xNStep, xMin is greater than xMax"
        end do

    end function find_xNStep

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_Permutation_from_Mat(pos, matrix, nDim, pVec);

        implicit none

        !INPUT
        integer                        , intent(in)           :: pos;
        double precision, dimension(1:, 1:), intent(in)           :: matrix;
        integer, intent(in) :: nDim

        !OUTPUT
        double precision, dimension(1:), intent(out) :: pVec;
        !LOCAL VARIABLES
        integer :: i, j;
        integer :: seedStep;
        integer, dimension(:), allocatable :: nStep;

        allocate(nStep(nDim))

        nStep = size(matrix,1)

        do j = 1, nDim
            if (j /= nDim) seedStep = product(nStep(j+1:));
            if (j == nDim) seedStep = 1;
            i = cyclicMod(int((pos-0.9)/seedStep)+1, nStep(j))
            pVec(j) = matrix(i, j);
        end do

        deallocate(nStep)

    end subroutine get_Permutation_from_Mat

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
