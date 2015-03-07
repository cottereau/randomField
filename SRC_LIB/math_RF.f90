module math_RF

    use displayCarvalhol
    use mpi
    use write_Log_File
    use constants_RF

    !All logic and math routines

!    interface set_Extremes
!       module procedure set_ExtremesStructured,   &
!           set_ExtremesUnstruct
!    end interface set_Extremes

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function isPowerOf(num, power) result(isPower)

        implicit none
        !INPUT
        double precision, intent(in) :: num
        integer, intent(in) :: power

        !OUTPUT
        logical :: isPower
        !integer :: i
        double precision :: compNum

        compNum = 1.0D0
        isPower = .false.
        !i = 0;

        !write(*,*) "FLAG1"
        do while ((num - compNum) >= TOLERANCE)
            !write(*,*) "i = ", i
            !write(*,*) "compNum = ", compNum
            compNum = compNum * dble(power)
            !i = i+1
        end do
        !write(*,*) "FLAG2"

        if(abs(compNum-num) < TOLERANCE) isPower = .true.

    end function isPowerOf

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function cyclicMod(pos, base) result(resPos)

        implicit none
        !INPUT
        integer, intent(in) :: pos, base; !desired position
        !OUTPUT
        integer :: resPos;!result Position

        resPos = mod(pos,base);
        if(resPos == 0) resPos = base;

    end function cyclicMod

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function int_roundDown(dividend, divisor) result(quotient)

        implicit none

        !INPUT
        integer, intent(in)           :: dividend, divisor;
        !OUTPUT
        integer:: quotient;

        quotient = dividend/divisor
        if (dividend == divisor*(dividend/divisor)) quotient = quotient - 1

    end function int_roundDown

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine init_random_seed(seedIn, seedStart)
       ! POST: The seed for the random number generation method random_number() has been reset

        implicit none
        !INPUT
        integer, dimension(1:), optional, intent(in) :: seedIn
        integer, optional, intent(in) :: seedStart
        !LOCAL
        integer :: clock
        integer, dimension(:), allocatable :: seed

        if(present(seedIn)) then
            call random_seed(PUT = seedIn)
        else if (present(seedStart)) then
            call calculate_random_seed(seed, seedStart)
            call random_seed(PUT = seed)
            deallocate(seed)
        else
            call calculate_random_seed(seed)
            call random_seed(PUT = seed)
            deallocate(seed)
        end if

    end subroutine init_random_seed

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine reorder_vector(vector)
       ! POST: The seed for the random number generation method random_number() has been reset

        implicit none
        !INPUT
        double precision, dimension(1:), intent(inout) :: vector

        !LOCAL
        integer :: i
        integer :: minPos
        double precision :: temp
        logical, dimension(:), allocatable :: mask

        allocate(mask(size(vector)))

        minPos = MINLOC(vector, dim=1 )
        mask = .true.

        do i = 1, size(vector)
            minPos = MINLOC(vector, dim=1, mask =  mask)
            temp   =  vector(i)
            vector(i) = vector(minPos)
            vector(minPos) = temp
            mask(i) = .false.
        end do

        deallocate(mask)

    end subroutine reorder_vector


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine calculate_random_seed(seed, seedStart)

        implicit none
        !INPUT
        integer, optional, intent(in) :: seedStart
        !OUTPUT
        integer, dimension(:), allocatable, intent(out) :: seed
        !LOCAL
        integer :: i
        integer :: n
        integer :: clock, tempClock

        if(.not.allocated(seed)) then
            call random_seed(size = n)
            allocate(seed(n))
        end if
        call system_clock(COUNT=clock)
        tempClock = clock
        do while (clock == tempClock)
            call system_clock(COUNT=clock)
        end do

        if(present(seedStart)) then
            seed = 72 + seedStart*18 + 37*(/ (i - 1, i = 1, n) /)
        else
            seed = clock + 37*(/ (i - 1, i = 1, n) /)
        end if

    end subroutine calculate_random_seed

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_SequenceParam(axis, nPlane, nStep, beg, step, ends)

        implicit none

        !INPUT
        integer,               intent(in) :: axis, nPlane
        integer, dimension(1:), intent(in) :: nStep;

        !OUTPUT
        integer, intent(out) :: beg, step, ends;

        !LOCAL VARIABLES
        integer :: i, nDim, position, posPlane, gap;

        beg  = 0
        nDim = size(nStep)

        step = product(nStep(axis+1:));
        if (i == nDim) step = 1;
        gap = step*nStep(axis)
        beg = cyclicMod(nPlane, step) &
            + int((nPlane-0.9)/step)*gap
        ends = beg + step*(nStep(axis)-1)

    end subroutine get_SequenceParam

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_DeltaMatrix(xPoints, deltaMatrix)

        implicit none

        !INPUT
        double precision, dimension(:, :), intent(in)  :: xPoints;

        !OUTPUT
        double precision, dimension(:,:), intent(out) :: deltaMatrix;

        !LOCAL VARIABLES
        logical         , dimension(:, :)   , allocatable :: uniqMask;
        logical         , dimension(:)      , allocatable :: minMask;
        integer          :: pos, posNeigh;
        integer          :: i, j, k, nDim, nPoints, nFactors;
        double precision :: tolerance = 1.0d-6, increment;

        !Matrix of distances between points: distMatrix (point 1, point 2, dimension)
        nDim        = size(xPoints, 1)
        nPoints     = size(xPoints, 2)
        deltaMatrix = 0;

        allocate(minMask    (nPoints))
        allocate(uniqMask   (nDim, nPoints))

        !call DispCarvalhol(transpose(xPoints),"transpose(xPoints)", nColumns = 15)

        !Building the unicity mask
        do i = 1, nDim
            uniqMask(i,:) = .true.
            do j = 1, nPoints
                if(uniqMask(i,j)) then
                    do k = j + 1, nPoints
                        if(xPoints(i, k) == xPoints(i, j)) then
                            uniqMask(i,k) = .false.
                        end if
                    end do
                end if
            end do
        end do

        !call DispCarvalhol(transpose(uniqMask(:,1:20)),"transpose(uniqMask(:,1:20))", nColumns = 15)
        !call DispCarvalhol(transpose(uniqMask(:,:)),"transpose(uniqMask(:,:))", nColumns = 15)

        !Building the delta matrix
        deltaMatrix = 0
        do i = 1, nDim
            minMask  = uniqMask(i,:)
            nFactors = count(minMask)
            do j = 1, nFactors - 1
                pos          = minloc(xPoints(i,:), dim = 1, mask = minMask)
                minMask(pos) = .false.
                posNeigh     = minloc(xPoints(i,:), dim = 1, mask = minMask)
                increment    = (xPoints(i, posNeigh) - xPoints(i, pos)) / 2
                deltaMatrix(i, pos)      = deltaMatrix(i, pos     ) + increment
                deltaMatrix(i, posNeigh) = deltaMatrix(i, posNeigh) + increment

                !!Supposing the extremes symetrics
                !if (pos == 1) deltaMatrix(pos     , i) = deltaMatrix(pos     , i) + increment
                !if (pos == (nFactors - 1)) deltaMatrix(posNeigh, i) = deltaMatrix(posNeigh, i) + increment
            end do
        end do

        !call DispCarvalhol(deltaMatrix,"deltaMatrix", nColumns = 15)

        !Filling the rest o the delta matrix (repeated values)
        do i = 1, nDim
            uniqMask(i, :) = .true.
            do j = 1, nPoints
                if(uniqMask(i, j)) then
                    do k = j + 1, nPoints
                        if(xPoints(i, k) == xPoints(i, j)) then
                            deltaMatrix(i,k) = deltaMatrix(i,j)
                            uniqMask(i,k) = .false.
                        end if
                    end do
                end if
            end do
        end do


        if(allocated(minMask))  deallocate(minMask)
        if(allocated(uniqMask)) deallocate(uniqMask)

    end subroutine set_DeltaMatrix

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_DistMatrix(xPoints, distMatrix)

        implicit none

        !INPUT
        double precision, dimension(:, :), intent(in)  :: xPoints;

        !OUTPUT
        double precision, dimension(:,:,:), intent(out) :: distMatrix;

        !LOCAL VARIABLES
        integer          :: i, nDim, nPoints;

        nDim    = size(xPoints, 1)
        nPoints = size(xPoints, 2)

        !Building the distance Matrix
        do i = 1, nPoints
            !write(*,*) "i = ", i
            distMatrix(i, :, :) = transpose(xPoints);
        end do
        do i = 1, nDim
            distMatrix(:, :, i) = transpose(distMatrix(:, :, i)) - distMatrix(:, :, i)
        end do
        distMatrix(:, :, :) = abs(distMatrix(:, :, :))

    end subroutine set_DistMatrix

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine reorderToGlobal(all_RandField, all_xPoints, mapping)

        implicit none

        !INPUT

        !OUTPUT
        double precision, dimension(:, :), intent(inout) :: all_RandField;
        double precision, dimension(:, :), intent(inout) :: all_xPoints;
        integer         , dimension(:)   , intent(out), allocatable, optional :: mapping;

        !LOCAL VARIABLES
        integer          :: i, j, nPoints, nDim, Nmc, coefStart, pos;
        double precision :: delta, tol = 1d-10, minStep;
        integer , dimension(:), allocatable :: globIndex;
        logical , dimension(:), allocatable :: minMask;
        double precision, dimension(:, :), allocatable :: tempRandField;
        double precision, dimension(:, :), allocatable :: tempXPoints;

        nPoints  = size(all_xPoints, 1)
        nDim     = size(all_xPoints, 2)
        Nmc      = size(all_RandField , 2)

        allocate(globIndex(nPoints))
        allocate(minMask(nPoints))
        allocate(tempRandField(size(all_RandField,1), size(all_RandField,2)))
        allocate(tempXPoints(size(all_xPoints,1), size(all_xPoints,2)))

        minStep = 1000;
        coefStart = 1
        globIndex = 0

        do i = nDim, 1, -1
            minStep = 1000;
            do j = 1, nPoints-1
                delta = all_xPoints(j + 1, i) - all_xPoints(j, i)
                if ((abs(delta) < minStep) .and. (abs(delta) > tol)) minStep =  abs(delta)
            end do
            if (i == nDim) then
                globIndex(:) = nint((all_xPoints(:, i) - minval(all_xPoints(:,i)))/minStep) + 1
                coefStart    = maxval(globIndex(:))
            else
                globIndex(:) = coefStart &
                    * (nint((all_xPoints(:, i) - minval(all_xPoints(:,i)))/minStep)) &
                    + globIndex(:)
                coefStart    = coefStart * maxval(nint(all_xPoints(:, i)/minStep))
            end if
        end do


        !Compacting globIndex and reordering
        tempRandField = all_RandField
        tempXPoints   = all_xPoints

        if(present(mapping)) allocate(mapping(nPoints))
        minMask = .true.
        do i = 1, nPoints
            pos = minloc(globIndex, 1, minMask)
            minMask(pos)   = .false.
            globIndex(pos) = i
            if(present(mapping)) mapping(i) = pos
            all_RandField(i,:) = tempRandField(pos,:)
            all_xPoints(i,:)   = tempXPoints(pos,:)
        end do

        deallocate(globIndex)
        deallocate(minMask)
        deallocate(tempRandField)
        deallocate(tempXPoints)

    end subroutine reorderToGlobal

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine calculateAverageCorrL(randField, xRange, xNStep, averageCorrL)
        implicit none
        !INPUT
        double precision, dimension(:, :), intent(in) :: randField;
        double precision, dimension(:),    intent(in) :: xRange;
        integer,          dimension(:),    intent(in) :: xNStep;
        !OUTPUT
        double precision, dimension(:),    intent(out) :: averageCorrL;
        !LOCAL VARIABLES
        double precision, dimension(:),    allocatable :: ptAvg, ptStdDev;
        integer :: i, j, k, Nmc, beg, step, end;
        integer :: nPoints, nPlanes, nDim;
        double precision :: corrL
        !double precision, dimension(:), allocatable :: Ra, Rb, Rc, R

        write(*,*) ">>>> Calculating comparation corrL (global)"

        !        call DispCarvalhol(randField, "randField")
        !call DispCarvalhol(xRange, "xRange")
        !call DispCarvalhol(xNStep, "xNStep")

        nPoints      = size(randField, 1);
        Nmc          = size(randField, 2);
        nDim         = size(xNStep);
        averageCorrL = 0;

        !        write(*,*) "nPoints = ", nPoints
        !        write(*,*) "Nmc = ", Nmc
        !        write(*,*) "nDim = ", nDim
        !        write(*,*) "averageCorrL = ", averageCorrL

        !        write(*,*) "Before Allocation"
        allocate(ptAvg(nPoints))
        allocate(ptStdDev(nPoints))

        ptAvg    = sum(randField, 2)/Nmc
        ptStdDev = sqrt(sum(randField**2, 2)/Nmc - (ptAvg)**2)

        do i = 1, nDim
            nPlanes = nPoints/xNStep(i)
            do j = 1, nPlanes
                call get_SequenceParam(i, j, xNStep, beg, step, end)
                !
                if(end > nPoints) end = nPoints

                averageCorrL(i) = sum(((1d0/dble(Nmc) * matmul(randField(beg:end:step, :),randField(beg,:))) &
                    - (ptAvg(beg:end:step)    * ptAvg(beg)))                               &
                    / (ptStdDev(beg:end:step) * ptStdDev(beg)))                            &
                    + averageCorrL(i)

                !                averageCorrL(i) = sum(((1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),randField(beg,:))) &
                !                                      - (ptAvg(beg+step:end:step)    * ptAvg(beg)))                               &
                !                                      / (ptStdDev(beg+step:end:step) * ptStdDev(beg)))                            &
                !                                       + averageCorrL(i)

                !Ra = 1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),&
                !                            randField(beg,:))
                !Rb = average(beg+step:end:step)*average(beg)
                !Rc = stdDeviation(beg+step:end:step)*stdDeviation(beg)
                !R = (Ra - Rb) / Rc
                !averageCorrL(i) = sum(R) * xMax(i,1)/xNStep(i,1) + averageCorrL(i)

                !call dispCarvalhol(1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),&
                !                            randField(beg,:)),"Ra")
                !call dispCarvalhol(average(beg+step:end:step)*average(beg),"Rb")
                !call dispCarvalhol(stdDeviation(beg+step:end:step)*stdDeviation(beg),"Rc")
                !call dispCarvalhol(((1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),randField(beg,:))) - &
                !                      (average(beg+step:end:step)*average(beg))) / &
                !                      (stdDeviation(beg+step:end:step)*stdDeviation(beg)),"R")
                !write(*,*) "R"
                !write(*,*)(((1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),randField(beg,:))) - &
                !                      (average(beg+step:end:step)*average(beg))) / &
                !                      (stdDeviation(beg+step:end:step)*stdDeviation(beg)))
                !write(*,*) "averageCorrL(",i,") = ", averageCorrL(i)

            end do
            averageCorrL(i) = (xRange(i)/xNStep(i))*averageCorrL(i) / dble(nPlanes)
        end do
        averageCorrL = 2*averageCorrL !Symmetry

        if(allocated(ptAvg)) deallocate(ptAvg)
        if(allocated(ptStdDev)) deallocate(ptStdDev)

    end subroutine calculateAverageCorrL

end module math_RF
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