module mesh_RF

    !use mpi
    use math_RF
    use write_Log_File
    use type_RF
    use type_MESH
    use fftw3

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
        integer :: d

        allocate(xPoints(MSH%nDim, MSH%xNTotal))

        call wLog("MSH%xMinBound = ")
        call wLog(MSH%xMinBound)
        call wLog("MSH%xStep = ")
        call wLog(MSH%xStep)
        call wLog("MSH%xNStep = ")
        call wLog(MSH%xNStep)
        call wLog("MSH%xNTotal = ")
        call wLog(MSH%xNTotal)

        call setGrid(xPoints, MSH%xMinBound, MSH%xStep, MSH%xNStep)

        RDF%xPoints => xPoints

        if(RDF%nDim == 2) then
            RDF%xPoints_2D(1:MSH%nDim, 1:MSH%xNStep(1), 1:MSH%xNStep(2)) => xPoints
            !call DispCarvalhol(RDF%xPoints_2D(1,:,:), "RDF%xPoints_2D 1", unit_in=RDF%log_ID)
            !call DispCarvalhol(RDF%xPoints_2D(2,:,:), "RDF%xPoints_2D 2", unit_in=RDF%log_ID)
        else if(RDF%nDim == 3) then
            RDF%xPoints_3D(1:MSH%nDim, 1:MSH%xNStep(1), 1:MSH%xNStep(2), 1:MSH%xNStep(3)) => xPoints
            !call wLog("Point in minimal position = ")
            !call wLog(RDF%xPoints_3D(:,22,1,1))
        end if

        !call DispCarvalhol(transpose(RDF%xPoints), "transpose(RDF%xPoints)", unit_in=RDF%log_ID)

    end subroutine set_XPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine setGrid(xPoints, xMinBound, xStep, xNStep, inverse)

        implicit none

        !INPUT
        double precision, dimension(:)  , intent(in) :: xMinBound
        integer, dimension(:), intent(in) :: xNStep
        double precision, dimension(:), intent(in) :: xStep
        logical, intent(in), optional :: inverse
        !OUTPUT
        double precision, dimension(:,:), intent(out) :: xPoints

        !LOCAL
        integer(kind=8) :: totalSize
        integer :: nDim
        integer(kind=8) :: sizePattern, unityMult, patternMult
        integer :: start, end, i, d
        logical :: effec_Inverse

        totalSize = product(xNStep)
        nDim = size(xNStep)

        if(product(xNStep) /= size(xPoints,2)) then
            write(*,*) "ERROR, inside set Grid shape(xPoints) and xNStep are different"
            write(*,*) "shape(xPoints) = ", shape(xPoints)
            write(*,*) "xNStep         = ", xNStep
            stop(" ")
        end if

        effec_Inverse = .false.
        if(present(inverse)) effec_Inverse =inverse

        if(effec_Inverse) then
            do d = nDim, 1, -1
                sizePattern = product(xNStep(d:nDim))
                unityMult   = sizePattern/xNStep(d)
                patternMult = totalSize/sizePattern

                !Building the basic pattern
                do i=1, xNStep(d)
                    start = (i-1)*unityMult + 1
                    end   = start + unityMult - 1
                    xPoints(d, start:end) = xMinBound(d) + xStep(d)*dble(i -1)
                end do

                !Replicating the pattern
                do i=2, patternMult
                    start = (i-1)*sizePattern + 1
                    end   = start + sizePattern - 1
                    xPoints(d, start:end) = xPoints(d, 1:sizePattern)
                end do
            end do
        else
            do d = 1, nDim
                sizePattern = product(xNStep(1:d))
                unityMult   = sizePattern/xNStep(d)
                patternMult = totalSize/sizePattern

                !Building the basic pattern
                do i=1, xNStep(d)
                    start = (i-1)*unityMult + 1
                    end   = start + unityMult - 1
                    xPoints(d, start:end) = xMinBound(d) + xStep(d)*dble(i -1)
                end do

                !Replicating the pattern
                do i=2, patternMult
                    start = (i-1)*sizePattern + 1
                    end   = start + sizePattern - 1
                    xPoints(d, start:end) = xPoints(d, 1:sizePattern)
                end do
            end do
        end if


    end subroutine setGrid


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
