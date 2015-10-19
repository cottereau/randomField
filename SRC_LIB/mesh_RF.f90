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
        integer, dimension(MSH%nDim) :: tempXNStep

        !LOCAL
        integer :: i, j, counterXPoints
        !integer, dimension(MSH%nDim) :: tempXNStep

        if(MSH%independent) then

            !ALLOCATING (Internal + External)

            call wLog("GENERATION BY VOLUME (Independent)")
            call wLog("MSH%xNStep ")
            call wLog(int(MSH%xNStep))
            call wLog("MSH%xNTotal ")
            call wLog(MSH%xNTotal )
            call wLog("RDF%xNTotal ")
            call wLog(int(RDF%xNTotal))

            call wLog("TOTAL")
            call wLog("MSH%xNStep")
            call wLog(MSH%xNStep)

            allocate(xPoints(MSH%nDim, MSH%xNTotal))

            call wLog("shape(xPoints)    ")
            call wLog(shape(xPoints))
            call wLog(" ")

            counterXPoints = 0;

            !Internal Points
            tempXNStep = find_xNStep(MSH%xMinInt, MSH%xMaxInt, MSH%xStep)
            do i = 1, product(tempXNStep)
                call get_Permutation(i, MSH%xMaxInt, tempXNStep, &
                    xPoints(:,i+counterXPoints), MSH%xMinInt,  &
                    snapExtremes = .true.);
            end do
            counterXPoints = counterXPoints + product(tempXNStep);

            call wLog("INTERNAL")
            call wLog("tempXNStep")
            call wLog(tempXNStep)
            call wLog(" ")

            !Border Points

            call wLog("NEIGHBOURS")
            do j = 1, size(MSH%xMaxNeigh, 2)

                if(.not. MSH%considerNeighbour(j)) cycle

                tempXNStep = find_xNStep(MSH%xMinNeigh(:,j), MSH%xMaxNeigh(:,j), MSH%xStep)

                call wLog("tempXNStep")
                call wLog(tempXNStep)

                do i = 1, product(tempXNStep)
                    call get_Permutation(i, MSH%xMaxNeigh(:,j), tempXNStep, &
                        xPoints(:,counterXPoints + i), MSH%xMinNeigh(:,j), &
                        snapExtremes = .true.);
                end do

                MSH%indexNeigh(1, j) = counterXPoints + 1
                counterXPoints = counterXPoints + product(tempXNStep);
                MSH%indexNeigh(2,j) = counterXPoints
            end do

            call wLog(" ")

            RDF%xPoints => xPoints

        else
            call wLog("LINEAR GENERATION")

            allocate(xPoints(MSH%nDim, MSH%xNTotal))

            counterXPoints = 0

            call wLog("MSH%xNInit")
            call wLog(MSH%xNInit)
            call wLog("MSH%xNEnd")
            call wLog(MSH%xNEnd)
            call wLog("MSH%xNTotal")
            call wLog(MSH%xNTotal)

            do i = MSH%xNInit, MSH%xNEnd
                counterXPoints = counterXPoints +1
                call get_Permutation(i, MSH%xMaxGlob, MSH%xNStep, &
                    xPoints(:,counterXPoints), MSH%xMinGlob,  &
                    snapExtremes = .true.);
            end do
            call wLog(" ")

            RDF%xPoints => xPoints

        end if

    end subroutine set_XPoints

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
        integer, dimension (MSH%nDim) :: tempXNStep

        periods(:) = .false.

        if (MSH%meshMod == "unv") then
            RDF%xPoints => coordList
            call wLog("-> defining_UNV_extremes")
            call get_Global_Extremes_Mesh(RDF%xPoints, MSH%xMinGlob, MSH%xMaxGlob, RDF%comm)
            RDF%xNTotal = MSH%xNTotal
            RDF%xMinExt = MSH%xMinGlob
            RDF%xMaxExt = MSH%xMaxGlob
            MSH%xMinExt = MSH%xMinGlob
            MSH%xMaxExt = MSH%xMaxGlob
            MSH%xNTotal = size(RDF%xPoints, 2)
            RDF%xNTotal = MSH%xNTotal
            MSH%xMinBound = MSH%xMinGlob
            MSH%xMaxBound = MSH%xMaxGlob
            call wLog("     MSH%xMinBound = ")
            call wLog(MSH%xMinBound)
            call wLog("     MSH%xMaxBound = ")
            call wLog(MSH%xMaxBound)
            call wLog(" ")

        else if(RDF%independent) then
            call wLog("-> set_procPerDim")
            call set_procPerDim (MSH%nb_procs, MSH%nDim, MSH%procPerDim)
            call wLog("-> MPI_CART_CREATE")
            call MPI_CART_CREATE (MSH%comm, MSH%nDim, MSH%procPerDim, periods, .false., MSH%topComm, code)
            call wLog("-> MPI_CART_COORDS")
            call MPI_CART_COORDS (MSH%topComm, MSH%rang, MSH%nDim, MSH%coords, code)

            call wLog("-> set_neighbours")
            call set_neighbours (MSH)
            call wLog("-> get_NeighbourCriteria")
            call get_NeighbourCriteria (MSH)

            call wLog("-> define_generation_geometry")
            call define_generation_geometry (MSH, RDF)
            call wLog("-> Getting Global Matrix Reference")
            call get_XPoints_globCoords(RDF, MSH)
            call wLog("     RDF%origin = ")
            call wLog(RDF%origin)
            call wLog(" ")

        else
            call wLog("-> define_generation_geometry")
            call define_generation_geometry (MSH, RDF)
            tempXNStep = find_xNStep(MSH%xMinGlob, MSH%xMaxGlob, MSH%xStep)
            call get_Permutation(MSH%xNInit, MSH%xMaxGlob, tempXNStep, &
                                 MSH%xMinBound, MSH%xMinGlob,  &
                                 snapExtremes = .true.)
            call get_Permutation(MSH%xNEnd, MSH%xMaxGlob, tempXNStep, &
                                 MSH%xMaxBound, MSH%xMinGlob,  &
                                 snapExtremes = .true.);
            call wLog("     MSH%xMinBound (Linear)= ")
            call wLog(MSH%xMinBound)
            call wLog("     MSH%xMaxBound (Linear) = ")
            call wLog(MSH%xMaxBound)
            call wLog(" ")
            call wLog("-> Getting Global Matrix Reference")
            call get_XPoints_globCoords(RDF, MSH)
            call wLog("     RDF%origin = ")
            call wLog(RDF%origin)
            call wLog(" ")

        end if

    end subroutine define_topography

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine define_generation_geometry (MSH, RDF)

        implicit none

        !INPUT AND OUTPUT
        type(MESH), intent(inout) :: MSH
        type(RF)  , intent(inout) :: RDF

        !LOCAL
        double precision, dimension(MSH%nDim) :: delta, half, ovlp, nonOvlp, minVol
        double precision, dimension(MSH%nDim) :: locNO, locO
        integer, dimension(MSH%nDim) :: nPointsNO, nPointsO, nPointsTot, nPointsLoc
        integer, dimension(MSH%nDim) :: xNStepGlob
        integer :: sliceSize

        !FFTW
        integer(C_INTPTR_T) :: L, M, N
        integer(C_INTPTR_T) :: local_LastDim
        integer(C_INTPTR_T) :: local_j_offset
        integer(C_INTPTR_T) :: alloc_local

        !GLOBAL MANIPULATIONS

        !Defining MSH%xStep
        call wLog(" ")
        call wLog("     IN  RDF%corrL          = ")
        call wLog(RDF%corrL)
        call wLog("     IN  MSH%pointsPerCorrL = ")
        call wLog(MSH%pointsPerCorrL)

        MSH%xStep = RDF%corrL/dble(MSH%pointsPerCorrL)

        call wLog("     OUT MSH%xStep          = ")
        call wLog(MSH%xStep)
        call wLog(" ")

        !Rounding Global Extremes

        call wLog(" ")
        call wLog("        IN  MSH%xMinGlob = ")
        call wLog(MSH%xMinGlob)
        call wLog("        IN  MSH%xMaxGlob = ")
        call wLog(MSH%xMaxGlob)
        call wLog("            delta        = ")
        call wLog(MSH%xMaxGlob - MSH%xMinGlob)

        call roundToMultiple(MSH%xMinGlob, MSH%xStep, up=.false.)
        call roundToMultiple(MSH%xMaxGlob, MSH%xStep, up=.true.)

        call wLog("        OUT MSH%xMinGlob (Round 1) = ")
        call wLog(MSH%xMinGlob)
        call wLog("        OUT MSH%xMaxGlob (Round 1) = ")
        call wLog(MSH%xMaxGlob)
        call wLog("            delta        (Round 1) = ")
        call wLog(MSH%xMaxGlob - MSH%xMinGlob)
        call wLog(" ")

        !GLOBAL CASE (Not INDEPENDENT)
        if(.not. MSH%independent) then

            !Global extremes
            half  = (MSH%xMaxGlob + MSH%xMinGlob)/2.0D0
            delta = MSH%xMaxGlob - MSH%xMinGlob
            call roundToMultiple(delta, MSH%xStep*dble(MSH%procPerDim)*2.0D0, up=.true.)
            MSH%xMinGlob = half - delta/2.0D0
            MSH%xMaxGlob = half + delta/2.0D0

            call wLog("        OUT MSH%xMinGlob (Round 2) = ")
            call wLog(MSH%xMinGlob)
            call wLog("        OUT MSH%xMaxGlob (Round 2) = ")
            call wLog(MSH%xMaxGlob)
            call wLog("            delta        (Round 2) = ")
            call wLog(MSH%xMaxGlob - MSH%xMinGlob)
            call wLog(" ")


            xNStepGlob = find_xNStep(xMaxExt=(MSH%xMaxGlob - MSH%xMinGlob), xStep=MSH%xStep)

            if(RDF%method == FFT) then

                if(xNStepGlob(MSH%nDim) < MSH%nb_procs) then
                   xNStepGlob(MSH%nDim) = MSH%nb_procs
                   half  = (MSH%xMaxGlob + MSH%xMinGlob)/2.0D0
                   delta(MSH%nDim) = xNStepGlob(MSH%nDim)
                   call roundToMultiple(delta, MSH%xStep*dble(MSH%procPerDim)*2.0D0, up=.true.)
                   MSH%xMinGlob = half - delta/2.0D0
                   MSH%xMaxGlob = half + delta/2.0D0

                   call wLog("        OUT MSH%xMinGlob (Round 3 - For FFT) = ")
                   call wLog(MSH%xMinGlob)
                   call wLog("        OUT MSH%xMaxGlob (Round 3 - For FFT) = ")
                   call wLog(MSH%xMaxGlob)
                   call wLog("            delta        (Round 3 - For FFT) = ")
                   call wLog(MSH%xMaxGlob - MSH%xMinGlob)
                   call wLog(" ")

                   xNStepGlob = find_xNStep(xMaxExt=(MSH%xMaxGlob - MSH%xMinGlob), xStep=MSH%xStep)
                end if

                call wLog("xNStepGlob = ")
                call wLog(xNStepGlob)

                if(xNStepGlob(MSH%nDim) < MSH%nb_procs) stop("ERROR!! When using parallel FFT the last dimension should have at least 1 slice by proc")

                if(RDF%nDim == 2) then
                    L = xNStepGlob(1)
                    M = xNStepGlob(2)
                    alloc_local = fftw_mpi_local_size_2d(M, L, RDF%comm, &
                                                         local_LastDim, local_j_offset) !FOR MPI
                else if(RDF%nDim == 3) then
                    L = xNStepGlob(1)
                    M = xNStepGlob(2)
                    N = xNStepGlob(3)
                    alloc_local = fftw_mpi_local_size_3d(N, M, L, RDF%comm, &
                                                         local_LastDim, local_j_offset) !FOR MPI
                else
                    stop("Inside define_generation_geometry no mesh division for FFT in this dimension")
                end if

                call wLog("local_LastDim = ")
                call wLog(int(local_LastDim))
                call wLog("local_j_offset = ")
                call wLog(int(local_j_offset))

                MSH%xNGlob = product(xNStepGlob)
                MSH%xNInit = local_j_offset + 1
                MSH%xNEnd  = MSH%xNInit + local_LastDim - 1

                sliceSize = 1
                if(MSH%nDim > 1) sliceSize = product(xNStepGlob(1:RDF%nDim -1))
                MSH%xNInit = (MSH%xNInit - 1) * sliceSize + 1
                MSH%xNEnd  =  MSH%xNEnd *sliceSize

            else
                MSH%xNGlob = product(xNStepGlob)

                MSH%xNEnd = (MSH%rang + 1) * MSH%XNGlob/MSH%nb_procs + 1
                MSH%xNInit  = MSH%rang * MSH%XNGlob/MSH%nb_procs + 1

                if(MSH%rang == MSH%nb_procs-1) then
                    MSH%xNEnd = MSH%XNGlob
                else
                    MSH%xNEnd = MSH%xNEnd - 1
                end if
            end if

            MSH%xNStep = find_xNStep(MSH%xMinGlob, MSH%xMaxGlob, MSH%xStep)
            MSH%xNTotal = MSH%xNEnd - MSH%xNInit + 1
            RDF%xNTotal = MSH%xNTotal

            call wLog("        OUT MSH%xNInit = ")
            call wLog(int(MSH%xNInit))
            call wLog("        OUT MSH%xNEnd = ")
            call wLog(int(MSH%xNEnd))
            call wLog("        OUT MSH%xNTotal = ")
            call wLog(int(MSH%xNTotal))
            call wLog("        OUT RDF%xNTotal = ")
            call wLog(int(RDF%xNTotal))
            call wLog("        OUT MSH%xNTotal = ")
            call wLog(int(MSH%xNTotal))
            call wLog("        OUT MSH%xNStep = ")
            call wLog(int(MSH%xNStep))

        !INDEPENDENT CASE (For Localization)
        else if(MSH%independent) then


            !Rounding overlap
            call wLog(" ")
            call wLog("     IN   MSH%overlap =  ")
            call wLog(MSH%overlap)

            if(MSH%independent) then
                MSH%overlap = ceiling(MSH%overlap*RDF%corrL/(2*MSH%xStep)) * 2*MSH%xStep/RDF%corrL
            else
                MSH%overlap = 0
            end if

            call roundToMultiple(MSH%overlap, MSH%xStep, up=.true.)

            call wLog("     OUT  MSH%overlap =  ")
            call wLog(MSH%overlap)
            call wLog(" ")


            !Size of the non-overlapping areas
            ovlp    = MSH%overlap*RDF%corrL*dble(MSH%procPerDim-1)
            nonOvlp = (MSH%xMaxGlob - MSH%xMinGlob) - ovlp
            where(nonOvlp < 0.0D0)
                nonOvlp = 0.0D0
            end where
            call roundToMultiple(nonOvlp, MSH%xStep*dble(MSH%procPerDim), up=.true.)

            call wLog("        OVERLAP TOTAL = ")
            call wLog(ovlp)
            call wLog("        NON-OVERLAP TOTAL = ")
            call wLog(nonOvlp)

            nPointsO   = find_xNStep(xMaxExt=MSH%overlap*RDF%corrL, xStep=MSH%xStep) - 2 !The extremes are considered in the non-overlapping area
            nPointsNO  = find_xNStep(xMaxExt=nonOvlp/dble(MSH%procPerDim), xStep=MSH%xStep)
            where(MSH%procPerDim == 1) nPointsO = 0

            call wLog("        nPointsO = ")
            call wLog(nPointsO)
            call wLog("        nPointsNO = ")
            call wLog(nPointsNO)

            if(any(nPointsO<1 .and. MSH%procPerDim > 1) .or. any(nPointsNO<1)) stop("Step or Overlap not adapted to the points per Correlation asked")

            !Global Extremes

            nPointsTot = (nPointsO * (MSH%procPerDim-1)) + (nPointsNO * MSH%procPerDim)
            call wLog("        nPointsTot = ")
            call wLog(nPointsTot)

            delta = dble(nPointsTot-1) * MSH%xStep
            MSH%xMaxGlob = MSH%xMinGlob + delta

            call wLog("        OUT MSH%xMinGlob (Round 2) = ")
            call wLog(MSH%xMinGlob)
            call wLog("        OUT MSH%xMaxGlob (Round 2) = ")
            call wLog(MSH%xMaxGlob)
            call wLog("            delta        (Round 2) = ")
            call wLog(MSH%xMaxGlob - MSH%xMinGlob)
            call wLog(" ")

            !Local Extremes
            nPointsLoc = nPointsO + 1 + nPointsNO
            where(MSH%coords == 0 .or. MSH%coords == (MSH%procPerDim-1)) nPointsLoc = (nPointsO + 1)/2 + nPointsNO
            delta = dble(nPointsLoc-1)*MSH%xStep

            nPointsTot = 1
            where(MSH%coords /= 0)
                nPointsTot = (nPointsO+1)/2 + nPointsNO  !Adding the processor in 0 (Only a half overlap area)
                nPointsTot = nPointsTot + (nPointsO + nPointsNO)*(MSH%coords-1) !Adding the others processors
            end where

            MSH%xMinExt = MSH%xMinGlob + (dble(nPointsTot-1) * MSH%xStep)
            MSH%xMaxExt = MSH%xMinExt + delta

            !Copy to the internal boundary
            MSH%xMinInt   = MSH%xMinExt
            MSH%xMaxInt   = MSH%xMaxExt

            call wLog("        OUT MSH%xMinExt = ")
            call wLog(MSH%xMinInt)
            call wLog("        OUT MSH%xMaxExt = ")
            call wLog(MSH%xMaxInt)
            call wLog("        OUT MSH%xMinInt = ")
            call wLog(MSH%xMinInt)
            call wLog("        OUT MSH%xMaxInt = ")
            call wLog(MSH%xMaxInt)

            call wLog("-> get_overlap_geometry")
            call get_overlap_geometry (MSH, RDF%corrL)

            MSH%xNStep  = find_xNStep(MSH%xMinBound, MSH%xMaxBound, MSH%xStep)
            MSH%xNTotal = product(MSH%xNStep)
            RDF%xNTotal = MSH%xNTotal

            call wLog("        OUT RDF%xNTotal = ")
            call wLog(int(RDF%xNTotal))
            call wLog("        OUT MSH%xNTotal = ")
            call wLog(int(MSH%xNTotal))
            call wLog("        OUT MSH%xNStep = ")
            call wLog(int(MSH%xNStep))

        end if

    end subroutine define_generation_geometry

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

        call wLog(" ")
        call wLog("-> Redimensioning for Overlap")
        call wLog(" ")

        !Redimensioning the internal part
        call wLog("    Redimensioning the internal and external part")

        call wLog(" IN MSH%xMinInt")
        call wLog(MSH%xMinInt)
        call wLog(" IN MSH%xMaxInt")
        call wLog(MSH%xMaxInt)
        call wLog(" IN MSH%xMinExt")
        call wLog(MSH%xMinExt)
        call wLog(" IN MSH%xMaxExt")
        call wLog(MSH%xMaxExt)

        do neighPos = 1, 2*MSH%nDim

            if(MSH%neigh(neighPos) < 0) cycle

            where(MSH%neighShift(:,neighPos) < 0)
                MSH%xMinInt = MSH%xMinExt + MSH%overlap*corrL/2.0D0
                MSH%xMinExt = MSH%xMinExt - MSH%overlap*corrL/2.0D0
            elsewhere(MSH%neighShift(:,neighPos) > 0)
                MSH%xMaxInt = MSH%xMaxExt - MSH%overlap*corrL/2.0D0
                MSH%xMaxExt = MSH%xMaxExt + MSH%overlap*corrL/2.0D0
            end where

        end do

        call wLog(" OUT MSH%xMinInt")
        call wLog(MSH%xMinInt)
        call wLog(" OUT MSH%xMaxInt")
        call wLog(MSH%xMaxInt)
        call wLog(" OUT MSH%xMinExt")
        call wLog(MSH%xMinExt)
        call wLog(" OUT MSH%xMaxExt")
        call wLog(MSH%xMaxExt)
        call wLog(" ")

        !Dimensioning overlapping area
        call wLog("    Dimensioning neighbours limits")

        do neighPos = 1, size(MSH%neigh)
            if(MSH%neigh(neighPos) < 0) cycle

            where(MSH%neighShift(:,neighPos) > 0)
                MSH%xMaxNeigh(:,neighPos) = MSH%xMaxExt - MSH%overlap*corrL/2.0D0 !- MSH%xStep
                MSH%xMinNeigh(:,neighPos) = MSH%xMaxInt + MSH%xStep
                MSH%xOrNeigh(:,neighPos)  = MSH%xMaxInt
            elsewhere(MSH%neighShift(:,neighPos) < 0)
                MSH%xMaxNeigh(:,neighPos) = MSH%xMinInt - MSH%xStep
                MSH%xMinNeigh(:,neighPos) = MSH%xMinExt + MSH%xStep + MSH%overlap*corrL/2.0D0
                MSH%xOrNeigh(:,neighPos)  = MSH%xMinInt
            elsewhere
                MSH%xMaxNeigh(:,neighPos) = MSH%xMaxInt
                MSH%xMinNeigh(:,neighPos) = MSH%xMinInt
                MSH%xOrNeigh(:,neighPos)  = MSH%xMinInt
            end where

        end do


        !Dimensioning Bounding Box
        MSH%xMinBound = MSH%xMinInt
        MSH%xMaxBound = MSH%xMaxInt

        do i = 1, size(MSH%xMaxNeigh, 2)

            if(.not. MSH%considerNeighbour(i)) cycle

            where (MSH%xMinNeigh(:,i) < MSH%xMinBound) MSH%xMinBound = MSH%xMinNeigh(:,i)
            where (MSH%xMaxNeigh(:,i) > MSH%xMaxBound) MSH%xMaxBound = MSH%xMaxNeigh(:,i)
        end do

        call wLog(" OUT MSH%xMinBound")
        call wLog(MSH%xMinBound)
        call wLog(" OUT MSH%xMaxBound")
        call wLog(MSH%xMaxBound)
        call wLog(" ")

        call show_MESHneigh(MSH, " ", onlyExisting = .true., forLog = .true.)

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
        double precision, dimension(:), allocatable :: xMinInt, xMaxInt;
        integer :: effectComm

        if(present(communicator)) then
            effectComm = communicator
        else
            effectComm = MPI_COMM_WORLD
        end if

        nDim = size(xPoints, 1)


        allocate(xMinInt(nDim))
        allocate(xMaxInt(nDim))

        xMinInt = minval(xPoints, 2)
        xMaxInt = maxval(xPoints, 2)

        do i = 1, nDim
            call MPI_ALLREDUCE (xMinInt(i), xMinGlob(i), 1, MPI_DOUBLE_PRECISION, &
                MPI_MIN, effectComm, code)
            call MPI_ALLREDUCE (xMaxInt(i), xMaxGlob(i), 1, MPI_DOUBLE_PRECISION, &
                MPI_MAX, effectComm, code)
        end do

        deallocate(xMinInt)
        deallocate(xMaxInt)

    end subroutine get_Global_Extremes_Mesh

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function find_xNStep(xMinExt, xMaxExt, xStep) result (xNStep)
        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: xMaxExt, xStep;
        double precision, dimension(:), intent(in), optional :: xMinExt

        !OUTPUT
        integer, dimension(:), allocatable :: xNStep
        integer :: i

        allocate(xNStep(size(xStep)))

        if(present(xMinExt)) then
            xNStep = 1 + nint((xMaxExt-xMinExt)/(xStep));
        else
            !write(*,*) "xMin not present"
            xNStep = 1 + nint(xMaxExt/xStep);
        end if



        do i = 1, size(xStep)
            if(xNStep(i) < 1) then
                write(*,*) "ERROR!!! Inside find_xNStep, xMinExt is greater than xMaxExt"
                write(*,*) "xNStep = ", xNStep

                !write(get_fileId(),*) "ERROR!!! Inside find_xNStep, xMinExt is greater than xMaxExt"
                !write(get_fileId(),*) " xMaxExt = ", xMaxExt
                !write(get_fileId(),*) " xMinExt = ", xMinExt
                !write(get_fileId(),*) " xStep = ", xStep
                !write(get_fileId(),*) " xNStep = ", xNStep
                stop(" ")
            end if
        end do

    end function find_xNStep

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_XPoints_globCoords(RDF, MSH)
        implicit none

        !INPUT AND OUTPUT
        type(RF)  , intent(inout) :: RDF
        type(MESH), intent(inout) :: MSH


        RDF%origin = find_xNStep(MSH%xMinGlob, MSH%xMinBound , MSH%xStep)


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
        MSH%xNStep = find_xNStep(MSH%xMinInt, MSH%xMaxInt, MSH%xStep)
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
    subroutine get_NeighbourCriteria(MSH)

        !INPUT OUTPUT
        type(MESH), intent(inout) :: MSH

        MSH%considerNeighbour = .true.
        where(MSH%neigh < 0) MSH%considerNeighbour = .false.
        !where(minval(MSH%neighShift,1) < 0) MSH%considerNeighbour = .false.

        call wLog(" MSH%considerNeighbour = ")
        call wLog(MSH%considerNeighbour)

    end subroutine get_NeighbourCriteria

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
