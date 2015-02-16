module mesh_RF

    !use mpi
    use math_RF
    use write_Log_File

contains
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_XPoints(xMinInd, xMaxInd, xStep, xPoints, xMinGlob, xMaxGlob, cutExtremes)
        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: xMinInd, xMaxInd
        double precision, dimension(:), intent(in), optional :: xMinGlob, xMaxGlob
        logical, optional, intent(in) :: cutExtremes

        !OUTPUT
        double precision, dimension(:), intent(in) :: xStep
        double precision, dimension(:,:), allocatable, intent(OUT) :: xPoints;

        !LOCAL VARIABLES
        integer :: nDim, i, xNTotal;
        integer , dimension(:) , allocatable :: xNStep;
        double precision, dimension(:), allocatable :: xMin, xMax

        nDim = size(xStep)
        allocate(xNStep(nDim))
        allocate(xMin(nDim))
        allocate(xMax(nDim))

        write(get_fileId(),*) "-> Creating xPoints";

        xMin = xMinInd
        xMax = xMaxInd

        write(get_fileId(),*) " xMin = ", xMin;
        write(get_fileId(),*) " xMax = ", xMax;

        if(present(cutExtremes)) then
            if(cutExtremes) call cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
        end if

        do i = 1, size(xMax)
            if(dble(ceiling((xMax(i) - xMin(i))/xStep(i))) /= (xMax(i) - xMin(i))/xStep(i)) then
                write(*,*) "WARNING! xStep is not compatible with the extremes"
                write(*,*) "    It will create an incomplete mesh"
                write(*,*) "xMin = ", xMin
                write(*,*) "xMax = ", xMax
                write(*,*) "xStep = ", xStep
            end if
        end do

        !write(*,*) "-> Finding xNStep";
        xNStep = find_xNStep(xMin, xMax, xStep)

        !write(*,*) "  xStep = ", xStep;

        xNTotal = product(xNStep)
        allocate(xPoints(nDim, xNTotal))

        !write(*,*) "-> Filling xNTotal";
        do i = 1, xNTotal
            call get_Permutation(i, xMax, xNStep, xPoints(:,i), xMin, snapExtremes = .true.);
        end do

        deallocate(xNStep)
        deallocate(xMin)
        deallocate(xMax)

    end subroutine set_XPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine recalculate_xStep(xMin, xMax, xStep)
        implicit none
        !INPUT
        double precision, dimension(:), intent(in) :: xMin, xMax
        !OUTPUT
        double precision, dimension(:), intent(inout) :: xStep
        !LOCAL
        integer :: i

        do i = 1, size(xMax)
            !write(*,*) "dble(ceiling((xMax(i) - xMin(i))/xStep(i))) = ", dble(ceiling((xMax(i) - xMin(i))/xStep(i)))
            !write(*,*) "(xMax(i) - xMin(i))/xStep(i) = ", (xMax(i) - xMin(i))/xStep(i)
            if(dble(ceiling((xMax(i) - xMin(i))/xStep(i))) /= (xMax(i) - xMin(i))/xStep(i)) then
                write(get_fileId(),*) "WARNING! xStep (", i, ") was not compatible with the extremes"
                write(get_fileId(),*) "    It will be recalculated"
                write(get_fileId(),*) "xMin(", i, ") = ", xMin(i)
                write(get_fileId(),*) "xMax(", i, ") = ", xMax(i)
                write(get_fileId(),*) "old xStep (", i, ") = ", xStep(i)
                xStep(i) = (xMax(i) - xMin(i))/dble(ceiling((xMax(i) - xMin(i))/xStep(i)))
                write(get_fileId(),*) "new xStep (", i, ") = ", xStep(i)
            end if
        end do

    end subroutine recalculate_xStep

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine recalculate_corrL(xMin, xMax, corrL)
        implicit none
        !INPUT
        double precision, dimension(:), intent(in) :: xMin, xMax
        !OUTPUT
        double precision, dimension(:), intent(inout) :: corrL
        !LOCAL
        integer :: i

        do i = 1, size(xMax)
            !write(*,*) "dble(ceiling((xMax(i) - xMin(i))/xStep(i))) = ", dble(ceiling((xMax(i) - xMin(i))/xStep(i)))
            !write(*,*) "(xMax(i) - xMin(i))/xStep(i) = ", (xMax(i) - xMin(i))/xStep(i)
            if(dble(ceiling((xMax(i) - xMin(i))/corrL(i))) /= (xMax(i) - xMin(i))/corrL(i)) then
                write(get_fileId(),*) "WARNING! corrL (", i, ") was not compatible with the extremes"
                write(get_fileId(),*) "    It will be recalculated"
                write(get_fileId(),*) "xMin(", i, ") = ", xMin(i)
                write(get_fileId(),*) "xMax(", i, ") = ", xMax(i)
                write(get_fileId(),*) "old corrL (", i, ") = ", corrL(i)
                corrL(i) = (xMax(i) - xMin(i))/dble(ceiling((xMax(i) - xMin(i))/corrL(i)))
                write(get_fileId(),*) "new corrL (", i, ") = ", corrL(i)
            end if
        end do

    end subroutine recalculate_corrL
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
        call recalculate_xStep(xMin, xMax, xStep)

        write(*,*) "xStep = ", xStep
        write(*,*) "xMax  = ", xMax
        write(*,*) "xMin  = ", xMin
        write(*,*) "corrL = ", corrL
        write(*,*) "pointsPerCorrL = ", pointsPerCorrL

        if(present(cutExtremes)) then
            if(cutExtremes) call cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
        end if

        call set_XPoints(xMin, xMax, xStep, xPoints)

        if(present(cutExtremes)) then
            if(cutExtremes) call restoreBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
        end if

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
        integer                       , intent(in)    :: rang, nb_procs;
        double precision, dimension(1:), intent(in) :: xMaxGlob;
        double precision, dimension(1:), intent(in) :: xMinGlob;
        !OUTPUT
        double precision, dimension(1:), intent(out) :: xMax;
        double precision, dimension(1:), intent(out) :: xMin;

        !LOCAL VARIABLES
        integer :: i, j, testRang = 0;
        integer :: seedStep, nDim, basicStep;
        double precision, dimension(:), allocatable :: xProcDelta;

        nDim = size(xMin);
        allocate (xProcDelta(nDim))
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

        deallocate (xProcDelta)

    end subroutine set_Local_Extremes_Mesh

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
        !To adjust the extremes in each proc so we don't have the same point in two processors on the interfaces
        implicit none
        !INPUT
        double precision, dimension(1:), intent(in)  :: xMinGlob, xMaxGlob, xStep
        !OUTPUT
        double precision, dimension(1:), intent(inout) :: xMin, xMax;
        !LOCAL
        integer :: i

        write(get_fileId(),*) "-> Cutting Borders";
        do i = 1, size(xMax)
            if(.not. xMax(i) == xMaxGlob(i)) xMax(i) = xMax(i) - xStep(i)
        end do

    end subroutine cutBorders

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine restoreBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
        !To reverse the operation of "cutBorders"
        implicit none
        !INPUT
        double precision, dimension(1:), intent(in)  :: xMinGlob, xMaxGlob, xStep
        !OUTPUT
        double precision, dimension(1:), intent(inout) :: xMin, xMax;
        !LOCAL
        integer :: i

        write(get_fileId(),*) "-> Restoring Borders";
        do i = 1, size(xMax)
            if(.not. xMax(i) == xMaxGlob(i)) xMax(i) = xMax(i) + xStep(i)
        end do

    end subroutine restoreBorders

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

        xNStep = 1 + nint((xMax-xMin)/(xStep));

    end function

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_Permutation(pos, qmax, nStep, pVec, qmin, snapExtremes)

        implicit none

        !INPUT
        integer                        , intent(in)           :: pos;
        double precision, dimension(1:), intent(in)           :: qmax;
        double precision, dimension(1:), intent(in), optional :: qmin;
        integer,          dimension(1:), intent(in)           :: nStep;
        logical, optional, intent(in) :: snapExtremes
        !OUTPUT
        double precision, dimension(1:), intent(out) :: pVec;
        !LOCAL VARIABLES
        integer :: i, j;
        integer :: seedStep, nDim;
        double precision :: contrib

        nDim = size(nStep);
        contrib = 0.0d0

        if (present(snapExtremes)) then
            if(snapExtremes) then
                contrib = 1.0d0
            end if
        end if

        if (present(qmin)) then
            do j = 1, nDim
                seedStep = product(nStep(j+1:));
                if (j == nDim) seedStep = 1;
                i = cyclicMod(int((pos-0.9)/seedStep)+1, nStep(j))
                pVec(j) = (dble(i)-0.5d0-contrib/2.0d0)         &
                          *(qmax(j)-qmin(j))/(nStep(j)-contrib) &
                          + qmin(j);
            end do
        else
            do j = 1, nDim
                seedStep = product(nStep(j+1:));
                if (j == nDim) seedStep = 1;
                i = cyclicMod(int((pos-0.9)/seedStep)+1, nStep(j))
                pVec(j) = (dble(i)-0.5d0-contrib/2.0d0) &
                          *(qmax(j))/(nStep(j)-contrib); !qmin = 0
            end do
        end if

    end subroutine get_Permutation

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
