module mesh_RF

    !use mpi
    use math_RF

contains
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_XPoints(xMin, xMax, xStep, xPoints)
        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: xMin, xMax, xStep;

        !OUTPUT
        double precision, dimension(:,:), allocatable, intent(OUT) :: xPoints;

        !LOCAL VARIABLES
        integer :: nDim, i, xNTotal;
        integer , dimension(:) , allocatable :: xNStep;

        nDim    = size(xMax)
        allocate(xNStep(nDim))
        xNStep  = 1+ceiling((xMax-xMin)/(xStep));
        xNTotal = product(xNStep)
        allocate(xPoints(nDim, xNTotal))

        do i = 1, xNTotal
            call get_Permutation(i, xMax, xNStep, xPoints(:,i), xMin, snapExtremes = .true.);
        end do

        deallocate(xNStep)

    end subroutine set_XPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_Local_Extremes_Mesh (xMin, xMax, rang, nb_procs)
        !Find the Boundaries for the box in each processor when using automatic mesh
        implicit none

        !INPUT
        integer                       , intent(in)    :: rang, nb_procs;
        double precision, dimension(1:), intent(inout) :: xMax;
        double precision, dimension(1:), intent(inout) :: xMin;
        !OUTPUT

        !LOCAL VARIABLES
        integer :: i, j, testRang = 0;
        integer :: seedStep, nDim, basicStep;
        double precision, dimension(:), allocatable :: xProcDelta;

        nDim = size(xMin);
        allocate (xProcDelta(nDim))
        basicStep  = nint(dble(nb_procs)**(1.0d0/nDim))
        xProcDelta = (xMax-xMin)/basicStep

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
