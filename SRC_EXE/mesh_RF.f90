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
    subroutine set_XPoints_independent(MSH, RDF, xPoints)
        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH
        type(RF)   :: RDF
        double precision, dimension(:, :), allocatable, intent(out), target :: xPoints;

        !LOCAL
        integer :: i, j, counterXPoints
        integer         , dimension(:), allocatable :: tempXNStep

        allocate(tempXNStep(MSH%nDim))

        !Snaping points to the grid and discover the bounding box
        call snap_to_grid(MSH, MSH%xMinLoc, MSH%xMaxLoc)
        MSH%xMaxBound = MSH%xMaxLoc;
        MSH%xMinBound = MSH%xMinLoc;
        do i = 1, size(MSH%xMaxNeigh, 2)

            if(MSH%neigh(i)<0) cycle

            call snap_to_grid(MSH, MSH%xMinNeigh(:,i), MSH%xMaxNeigh(:,i))
            do j = 1, MSH%nDim
                if (MSH%xMinNeigh(j,i) < MSH%xMinBound(j)) MSH%xMinBound(j) = MSH%xMinNeigh(j,i)
                if (MSH%xMaxNeigh(j,i) > MSH%xMaxBound(j)) MSH%xMaxBound(j) = MSH%xMaxNeigh(j,i)
            end do
        end do

        MSH%xNStep = find_xNStep(MSH%xMinBound, MSH%xMaxBound, MSH%xStep)
        MSH%xNTotal = product(MSH%xNStep)
        RDF%xNTotal = MSH%xNTotal

        allocate(xPoints(MSH%nDim, MSH%xNTotal))

        write(get_fileId(),*) "shape(xPoints) = ", shape(xPoints)


        !Internal Points
        counterXPoints = 0;
        tempXNStep = find_xNStep(MSH%xMinLoc, MSH%xMaxLoc, MSH%xStep)
        do i = 1, product(tempXNStep)
            call get_Permutation(i, MSH%xMaxLoc, tempXNStep, xPoints(:,i), MSH%xMinLoc, snapExtremes = .true.);
        end do

        counterXPoints = counterXPoints + product(tempXNStep);

        !Border Points
        do j = 1, size(MSH%xMaxNeigh, 2)
            if(MSH%neigh(j)<0) cycle
            tempXNStep = find_xNStep(MSH%xMinNeigh(:,j), MSH%xMaxNeigh(:,j), MSH%xStep)
            do i = 1, product(tempXNStep)
                call get_Permutation(i, MSH%xMaxNeigh(:,j), tempXNStep, xPoints(:,counterXPoints + i), MSH%xMinNeigh(:,j), snapExtremes = .true.);
            end do
            MSH%indexNeigh(1, j) = counterXPoints + 1
            counterXPoints = counterXPoints + product(tempXNStep);
            MSH%indexNeigh(2,j) = counterXPoints
        end do

        RDF%xPoints => xPoints

        deallocate(tempXNStep)

    end subroutine set_XPoints_independent

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

        !LOCAL VARIABLES
        double precision, dimension(:), allocatable :: xBottom, xTop
        integer :: i
        double precision :: norm

        write(get_fileId(),*) "-> Creating xPoints";

        allocate(xBottom(MSH%nDim))
        allocate(xTop(MSH%nDim))

        do i = 1, MSH%nDim
            norm = (MSH%xMaxGlob(i)-MSH%xMinGlob(i))/MSH%xStep(i)
            if(MSH%rang == 0 .and. (.not. areEqual(norm, dble(nint(norm))))) then
                write(*,*) "WARNING!! In dimension ", i, "extremes and step are incompatible"
                write(*,*) "   xMinGlob = ", MSH%xMinGlob(i)
                write(*,*) "   xMaxGlob = ", MSH%xMaxGlob(i)
                write(*,*) "   xStep    = ", MSH%xStep(i)
            end if
            xBottom(i) = (MSH%xStep(i) * intDivisor((MSH%xMin(i) - MSH%xMinGlob(i)), MSH%xStep(i), up = .true.)) &
                         + MSH%xMinGlob(i)
            xTop(i)    = (MSH%xStep(i) * intDivisor((MSH%xMax(i) - MSH%xMinGlob(i)), MSH%xStep(i), up = .false.)) &
                         + MSH%xMinGlob(i)
            if(areEqual(xTop(i), MSH%xMax(i))) then
                if(.not. areEqual(xTop(i), MSH%xMaxGlob(i))) then
                    xTop(i) = xTop(i) - MSH%xStep(i)
                end if
            end if
        end do

        write(get_fileId(),*) "xBottom = ", xBottom
        write(get_fileId(),*) "   xTop = ", xTop


        write(get_fileId(),*) "-> Finding xNStep";
        MSH%xNStep = find_xNStep(xBottom, xTop, MSH%xStep)
        MSH%xNTotal = product(MSH%xNStep)
        RDF%xNTotal = MSH%xNTotal

        allocate(xPoints(MSH%nDim, MSH%xNTotal))

        write(get_fileId(),*) "-> Filling xPoints";
        do i = 1, MSH%xNTotal
            call get_Permutation(i, xTop, MSH%xNStep, xPoints(:,i), xBottom, snapExtremes = .true.);
        end do

        RDF%xPoints => xPoints

        !call show_MESH(MSH)
        !call dispCarvalhol(transpose(RDF%xPoints), "transpose(RDF%xPoints)")

        deallocate(xBottom)
        deallocate(xTop)

    end subroutine set_XPoints

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

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_Local_Extremes_Mesh (xMin, xMax, xMinGlob, xMaxGlob, rang, nb_procs)
        !Find the Boundaries for the box in each processor when using automatic mesh
        implicit none

        !INPUT
        integer                        , intent(in) :: rang, nb_procs;
        double precision, dimension(1:), intent(in) :: xMaxGlob;
        double precision, dimension(1:), intent(in) :: xMinGlob;
        !OUTPUT
        double precision, dimension(1:), intent(out) :: xMax;
        double precision, dimension(1:), intent(out) :: xMin;

        !LOCAL VARIABLES
        integer :: i, j, testRang = 0;
        integer :: baseStep
        integer :: seedStep, nDim, basicStep;
        integer, allocatable, dimension(:) :: bStepVec
        double precision, dimension(:), allocatable :: xProcDelta;
        double precision :: procIterations

        nDim = size(xMin);
        baseStep = nint(dble(nb_procs)**(1.0d0/nDim))
        write(get_fileId(),*) "Setting Local Extremes"
        allocate (xProcDelta(nDim))


        if (nb_procs == baseStep**nDim) then
            write(get_fileId(),*) "Exact Division"
            basicStep  = nint(dble(nb_procs)**(1.0d0/nDim))
            xProcDelta = (xMaxGlob-xMinGlob)/basicStep

            !        if(rang == testRang) write(*,*) "nDim = ", nDim
            !        if(rang == testRang) write(*,*) "basicStep = ", basicStep
            !        if(rang == testRang) write(*,*) "xProcDelta = ", xProcDelta
            !        if(rang == testRang) write(*,*) "nb_procs = ", nb_procs
            !        if(rang == testRang) write(*,*) "dble(nb_procs)**(1/nDim) = ", dble(nb_procs)**(1/nDim)
            !        if(rang == testRang) write(*,*) "nint(dble(nb_procs)**(1/nDim)) = ", nint(dble(nb_procs)**(1/nDim))

            do j = 1, nDim
                seedStep = basicStep**(nDim-j);
                i = cyclicMod(int(rang/seedStep) + 1, basicStep)
                !if(rang == testRang) write(*,*) "i = ", i, " seedStep = ", seedStep, "seedStep", seedStep
                xMin(j) = (dble(i-1))*xProcDelta(j);
                xMax(j) = xMin(j) + xProcDelta(j)
                !if(rang == testRang) write(*,*) "AFTER xMin ", j, " = ", xMin(j)
            end do



            !write(*,*) "Proc = ", rang, "xMin = ", xMin!, "xMax = ", xMax;

        !write(*,*) "Proc = ", rang, "xMin = ", xMin!, "xMax = ", xMax;
        !write(*,*) "Proc = ", rang, "xMin = ", xMin!, "xMax = ", xMax;

        else if (isPowerOf(dble(nb_procs), 2)) then
            write(get_fileId(),*) "Power of two"
            allocate (bStepVec(nDim))

            !Defining the basic Step for each dimension
            bStepVec(:) = 1
            if(nb_procs /= 1) then
                procIterations = log(dble(nb_procs))/log(2.0D0)
                do j = 1, nint(procIterations)
                    i = cyclicMod(j, nDim)
                    bStepVec(i) = bStepVec(i)*2
                    !write(*,*) "i = ", i
                    !write(*,*) "bStepVec = ", bStepVec
                end do
            end if
            !Defining coordinates in each proc
            xProcDelta = (xMaxGlob-xMinGlob)/bStepVec

            write(get_fileId(),*) "procIterations = ", procIterations
            write(get_fileId(),*) "xMaxGlob = ", xMaxGlob
            write(get_fileId(),*) "xMinGlob = ", xMinGlob
            write(get_fileId(),*) "nb_procs = ", nb_procs
            write(get_fileId(),*) "bStepVec = ", bStepVec
            write(get_fileId(),*) "xProcDelta = ", xProcDelta

            do j = 1, nDim
                seedStep = product(bStepVec(j+1:));
                if (j == nDim) seedStep = 1;
                i = cyclicMod(int(rang/seedStep) + 1, bStepVec(j))
                !if(rang == testRang) write(*,*) "i = ", i, " seedStep = ", seedStep, "seedStep", seedStep
                xMin(j) = (dble(i-1))*xProcDelta(j);
                xMax(j) = xMin(j) + xProcDelta(j)
                !if(rang == testRang) write(*,*) "AFTER xMin ", j, " = ", xMin(j)
            end do
            deallocate (bStepVec)
        else
            stop "ERROR, no mesh division algorithm for this number of procs"
        end if

        do i = 1, nDim
            write(get_fileId(),*) "Dim ", i
            write(get_fileId(),fmt="(2A30)") "xMin", "xMax"
            write(get_fileId(),fmt="(2F30.15)") xMin(i), xMax(i)
        end do

        deallocate (xProcDelta)

    end subroutine set_Local_Extremes_Mesh

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_Local_Extremes_From_Coords (MSH)
        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH

        !LOCAL
        double precision, dimension(:), allocatable :: procDelta

        allocate(procDelta(MSH%nDim))

        procDelta = (MSH%xMaxGlob - MSH%xMinGlob)/MSH%procPerDim

        MSH%xMin = procDelta*MSH%coords + MSH%xMinGlob
        MSH%xMax = MSH%xMin + procDelta
        MSH%xMinLoc = MSH%xMin
        MSH%xMaxLoc = MSH%xMax

        write(get_fileId(),*) "MSH%xMin = ", MSH%xMin
        write(get_fileId(),*) "MSH%xMax = ", MSH%xMax

        deallocate(procDelta)

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


        write(get_fileId(),*) "Setting Local Extremes"
        !write(*,*) "Setting Local Extremes"

        procRootDim = dble(MSH%nb_procs)**(1/dble(MSH%nDim))
        logProc2   = log(dble(MSH%nb_procs))/log(2.0D0)

        if (areEqual(procRootDim, dble(nint(procRootDim)))) then
            write(get_fileId(),*) "Exact Division"
            !write(*,*) "Exact Division"
            MSH%procPerDim(:) = nint(dble(MSH%nb_procs)**(1.0d0/MSH%nDim))
        else if(areEqual(logProc2, dble(nint(logProc2)))) then
            write(get_fileId(),*) "Power of two"
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
    subroutine redefine_Global_Inputs (MSH, RDF, pointsPerCorrL)

        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH
        type(RF)   :: RDF
        integer, intent(in) :: pointsPerCorrL

        !LOCAL
        integer :: i
        integer :: testrank = 4

        MSH%xStep = RDF%corrL / dble(pointsPerCorrL);

        write(get_fileId(),*) "BEFORE:"
        write(get_fileId(),*) "MSH%xMinGlob = ", MSH%xMinGlob
        write(get_fileId(),*) "MSH%xMaxGlob = ", MSH%xMaxGlob
        write(get_fileId(),*) "MSH%overlap  = ", MSH%overlap

        !Redefining global extremes
        do i = 1, MSH%nDim
            MSH%xMinGlob(i) = dble(MSH%procPerDim(i))*MSH%xStep(i) * &
                              intDivisor(MSH%xMinGlob(i), MSH%xStep(i)*dble(MSH%procPerDim(i)), up = .false.)
            MSH%xMaxGlob(i) = dble(MSH%procPerDim(i))*MSH%xStep(i) * &
                              intDivisor(MSH%xMaxGlob(i), MSH%xStep(i)*dble(MSH%procPerDim(i)), up = .true.)
        end do


        !Redefining overlap
        MSH%overlap = 2*MSH%xStep(1) * intDivisor(MSH%overlap, MSH%xStep(1)*2, up = .true.)

        write(get_fileId(),*) "AFTER:"
        write(get_fileId(),*) "MSH%xMinGlob = ", MSH%xMinGlob
        write(get_fileId(),*) "MSH%xMaxGlob = ", MSH%xMaxGlob
        write(get_fileId(),*) "MSH%overlap  = ", MSH%overlap

    end subroutine redefine_Global_Inputs

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
    subroutine redefine_extremes (MSH, corrL)

        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH

        !INPUT
        double precision, dimension(:) :: corrL

        !LOCAL VARIABLES
        integer :: i, d
        integer :: code, neighPos;
        !integer, dimension(:), allocatable :: shift
        integer :: testrank = 4

        !allocate(shift(MSH%nDim))



        !call MPI_CART_COORDS (MSH%topComm, neigh, MSH%nDim, shift, code)

        !if(MSH%rang == testrank) write(*,*) "Inside redefine_extremes"

        !if(MSH%rang == testrank) call show_MESH(MSH, "BEFORE")

        !Redimensioning the internal part
        do neighPos = 1, 2*MSH%nDim
            !if(MSH%rang == testrank) write(*,*) "neighPos = ", neighPos
            if(MSH%neigh(neighPos) < 0) cycle

            !call MPI_CART_COORDS (MSH%topComm, MSH%neigh(neighPos), MSH%nDim, shift, code)

            !shift = shift - MSH%coords

            do i = 1, MSH%nDim
                if(MSH%neighShift(i,neighPos) < 0) then
                    MSH%xMinLoc(i) = MSH%xMin(i) + MSH%overlap*corrL(i)/2 + MSH%xStep(i)
                else if (MSH%neighShift(i,neighPos) > 0) then
                    MSH%xMaxLoc(i) = MSH%xMax(i) - MSH%overlap*corrL(i)/2 - MSH%xStep(i)
                end if

                !Checking if the points are in the range
                if(MSH%xMinLoc(i) < MSH%xMinGlob(i)) then
                    write(get_fileId(),*)"WARNING xMinLoc exceeded global dimensions (consider changing the overlap size or the global size)"
                    write(get_fileId(),*)"neigh = ",neighPos," dim = ", i
                    write(get_fileId(),*)"MSH%xMinLoc(i) = ", MSH%xMinLoc(i)
                    write(get_fileId(),*)"MSH%xMinGlob(i)    = ", MSH%xMinGlob(i)
                    MSH%xMinLoc(i) = MSH%xMinGlob(i)
                end if
                if(MSH%xMaxLoc(i) > MSH%xMaxGlob(i)) then
                    write(get_fileId(),*)"WARNING xMaxLoc exceeded global dimensions (consider changing the overlap size or the global size)"
                    write(get_fileId(),*)"neigh = ",neighPos," dim = ", i
                    write(get_fileId(),*)"MSH%xMaxLoc(i) = ", MSH%xMaxLoc(i)
                    write(get_fileId(),*)"MSH%xMaxGlob(i)    = ", MSH%xMaxGlob(i)
                    MSH%xMaxLoc(i) = MSH%xMaxGlob(i)
                end if

            end do
        end do

        !Dimensioning overlapping area
        do neighPos = 1, size(MSH%neigh)
            if(MSH%neigh(neighPos) < 0) cycle

!            call MPI_CART_COORDS (MSH%topComm, MSH%neigh(neighPos), MSH%nDim, shift, code)
!
!            shift = shift - MSH%coords

            do i = 1, MSH%nDim
                if(MSH%neighShift(i,neighPos) > 0) then
                    MSH%xMaxNeigh(i,neighPos) = MSH%xMax(i) + MSH%overlap*corrL(i)/2
                    MSH%xMinNeigh(i,neighPos) = MSH%xMax(i) - MSH%overlap*corrL(i)/2
                else if (MSH%neighShift(i,neighPos) < 0) then
                    MSH%xMaxNeigh(i,neighPos) = MSH%xMin(i) + MSH%overlap*corrL(i)/2
                    MSH%xMinNeigh(i,neighPos) = MSH%xMin(i) - MSH%overlap*corrL(i)/2
                else
                    MSH%xMaxNeigh(i,neighPos) = MSH%xMaxLoc(i)
                    MSH%xMinNeigh(i,neighPos) = MSH%xMinLoc(i)
                end if

                !Checking if the points are in the range
                if(MSH%xMinNeigh(i,neighPos) < MSH%xMinGlob(i)) then
                    write(get_fileId(),*)"WARNING!! MSH%xMinNeigh exceeded global dimensions (consider changing the overlap size or the global size)"
                    write(get_fileId(),*)"neigh = ",neighPos," dim = ", i
                    write(get_fileId(),*)"MSH%xMinNeigh(i,neighPos) = ", MSH%xMinNeigh(i,neighPos)
                    write(get_fileId(),*)"MSH%xMinGlob(i)    = ", MSH%xMinGlob(i)
                    MSH%xMinNeigh(i,neighPos) = MSH%xMinGlob(i)
                end if
                if(MSH%xMaxNeigh(i,neighPos) > MSH%xMaxGlob(i)) then
                    write(get_fileId(),*)"WARNING!! MSH%xMaxNeigh exceeded global dimensions (consider changing the overlap size or the global size)"
                    write(get_fileId(),*)"neigh = ",neighPos," dim = ", i
                    write(get_fileId(),*)"MSH%xMaxNeigh(i,neighPos) = ", MSH%xMaxNeigh(i,neighPos)
                    write(get_fileId(),*)"MSH%xMaxGlob(i)    = ", MSH%xMaxGlob(i)
                    MSH%xMaxNeigh(i,neighPos) = MSH%xMaxGlob(i)
                end if

            end do
        end do

        !if(MSH%rang == testrank) call show_MESH(MSH, "AFTER")

        !deallocate(shift)

    end subroutine redefine_extremes

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

        allocate(xNStep(size(xStep)))

!        write(*,*) "HERE find_xNStep"
!        write(*,*) " (xMax-xMin) = ", (xMax-xMin)
!        write(*,*) " (xStep) = ", (xStep)
!        write(*,*) " nint((xMax-xMin)/(xStep)) = ", nint((xMax-xMin)/(xStep))
!        write(*,*) " 1 + nint((xMax-xMin)/(xStep)) = ", 1 + nint((xMax-xMin)/(xStep))

        xNStep = 1 + nint((xMax-xMin)/(xStep));

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
