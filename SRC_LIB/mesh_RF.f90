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

        xPoints = dble(meshGridInt(MSH%xNStep,0))

        do d = 1, MSH%nDim
            xPoints(d,:) = MSH%xMinGlob(d) + MSH%xStep(d)*xPoints(d,:)
        end do

        RDF%xPoints => xPoints

    end subroutine set_XPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function meshGridInt(xNStep, nStart) result(intMatrix)

        implicit none

        !INPUT
        integer, dimension(:), intent(in)   :: xNStep
        !OUTPUT
        integer, dimension(:,:), allocatable :: intMatrix
        !LOCAL
        integer :: totalSize, nDim
        integer :: sizePattern, unityMult, patternMult
        integer :: start, end, i, d
        integer :: nStart

        totalSize = product(xNStep)
        nDim      = size(xNStep)

        allocate(intMatrix(nDim, totalSize))

        do d = 1, nDim
            sizePattern = product(xNStep(1:d))
            unityMult   = sizePattern/xNStep(d)
            patternMult = totalSize/sizePattern

            !Building the basic pattern
            do i=1, xNStep(d)
                start = (i-1)*unityMult + 1
                end   = start + unityMult - 1
                intMatrix(d, start:end) = i -1 + nStart
            end do

            !Replicating the pattern
            do i=2, patternMult
                start = (i-1)*sizePattern + 1
                end   = start + sizePattern - 1
                intMatrix(d, start:end) = intMatrix(d, 1:sizePattern)
            end do
        end do


    end function meshGridInt

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine snap_to_grid(MSH, extInf,extSup)
!        implicit none
!
!        !INPUT AND OUTPUT
!        type(MESH) :: MSH
!
!        !OUTPUT
!        double precision, dimension(:), intent(inout) :: extSup, extInf
!
!        !LOCAL VARIABLES
!        integer :: i
!
!        do i = 1, MSH%nDim
!            extInf(i) = (MSH%xStep(i) * dble(nint((extInf(i) - MSH%xMinGlob(i))/MSH%xStep(i)))) &
!                         + MSH%xMinGlob(i)
!            extSup(i) = (MSH%xStep(i) * dble(nint((extSup(i) - MSH%xMinGlob(i))/MSH%xStep(i)))) &
!                         + MSH%xMinGlob(i)
!        end do
!
!    end subroutine snap_to_grid

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
        MSH%xNStep = find_xNStep(MSH%xMinInt, MSH%xMaxInt, MSH%xStep)
        MSH%xNTotal = product(MSH%xNStep)
        RDF%xNTotal = MSH%xNTotal

        allocate(xPoints(MSH%nDim, MSH%xNTotal))

        RDF%xPoints => xPoints

    end subroutine allocate_xPoints


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
