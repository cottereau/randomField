module topography_RF

    !use mpi
    use math_RF
    use write_Log_File
    use type_RF
    use type_MESH
    use fftw3

    implicit none

contains

    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    subroutine round_basic_inputs(MSH, xStep, overlap)

        implicit none
        !INPUT
        type(MESH), intent(in) :: MSH

        !OUTPUT
        double precision, dimension(:), intent(inout) :: xStep, overlap

        !Defining MSH%xStep
        call wLog(" ")
        call wLog("     IN  MSH%corrL          = ")
        call wLog(MSH%corrL)
        call wLog("     IN  MSH%pointsPerCorrL = ")
        call wLog(MSH%pointsPerCorrL)

        xStep = MSH%corrL/dble(MSH%pointsPerCorrL)

        call wLog("     OUT xStep          = ")
        call wLog(xStep)
        call wLog(" ")

        !Rounding overlap
        call wLog(" ")
        call wLog("     IN   MSH%overlap =  ")
        call wLog(MSH%overlap)
        call wLog("     IN   MSH%xStep =  ")
        call wLog(MSH%xStep)
        call wLog("     IN   MSH%independent =  ")
        call wLog(MSH%independent)

        if(MSH%independent) then
            overlap = ceiling(overlap*MSH%corrL/(2*MSH%xStep)) * 2*MSH%xStep/MSH%corrL
        else
            overlap = 0.0D0
        end if

        where(MSH%procPerDim == 1) overlap = 0.0D0

        call wLog("     OUT  overlap =  ")
        call wLog(overlap)
        call wLog(" ")

    end subroutine round_basic_inputs
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    subroutine set_global_extremes(MSH,xMaxGlob, xMinGlob, procExtent, procStart)
        !INPUT
        type(MESH), intent(in) :: MSH
        !OUTPUT
        double precision, dimension(:), intent(inout) :: xMinGlob, xMaxGlob
        double precision, dimension(:), intent(out) :: procExtent, procStart
        !LOCAL
        double precision, dimension(MSH%nDim) :: delta, half, minVol
        double precision, dimension(MSH%nDim) :: ovlp_Tot, Novlp_Tot, ovlp, Novlp
        double precision, dimension(MSH%nDim) :: delta_Proc, stepProc
        double precision, dimension(MSH%nDim) :: locNO, locO
        integer, dimension(MSH%nDim) :: nPointsO, nPointsO_Ext, nPointsTot, nPointsLoc
        integer, dimension(MSH%nDim) :: nPointsNO, nPointsNO_int, nPointsNO_ext
        integer, dimension(MSH%nDim) :: pointInit, pointEnd, pointNb
        integer, dimension(MSH%nDim) :: xNStepGlob
        integer :: sliceSize
        integer :: i

        !Rounding Global Extremes
        call wLog(" ")
        call wLog("        IN  MSH%xMinGlob = ")
        call wLog(MSH%xMinGlob)
        call wLog("        IN  MSH%xMaxGlob = ")
        call wLog(MSH%xMaxGlob)
        call wLog("            delta        = ")
        call wLog(MSH%xMaxGlob - MSH%xMinGlob)

        call roundToMultiple(xMinGlob, MSH%xStep, up=.false.)
        call roundToMultiple(xMaxGlob, MSH%xStep, up=.true.)

        call wLog("        OUT MSH%xMinGlob (Round 1) = ")
        call wLog(MSH%xMinGlob)
        call wLog("        OUT MSH%xMaxGlob (Round 1) = ")
        call wLog(MSH%xMaxGlob)
        call wLog("            delta        (Round 1) = ")
        call wLog(MSH%xMaxGlob - MSH%xMinGlob)
        call wLog(" ")

        call wLog("        IN  MSH%procPerDim = ")
        call wLog(MSH%procPerDim)

        !GLOBAL CASE (Not INDEPENDENT)
        if(.not. MSH%independent) then

            !Global extremes
            half  = (MSH%xMaxGlob + MSH%xMinGlob)/2.0D0
            delta = MSH%xMaxGlob - MSH%xMinGlob
            call roundToMultiple(delta, MSH%xStep*dble(MSH%procPerDim), up=.true.)
            xMaxGlob = MSH%xMinGlob + delta

            procExtent = delta/dble(MSH%procPerDim)
            procStart  = procExtent*MSH%coords + MSH%xMinGlob

            call wLog("        OUT MSH%xMinGlob (Round 2) = ")
            call wLog(MSH%xMinGlob)
            call wLog("        OUT MSH%xMaxGlob (Round 2) = ")
            call wLog(MSH%xMaxGlob)
            call wLog("            delta        (Round 2) = ")
            call wLog(MSH%xMaxGlob - MSH%xMinGlob)
            call wLog(" ")

        !INDEPENDENT CASE (For Localization)
        else if(MSH%independent) then


            !Size of the non-overlapping areas
            ovlp      = MSH%overlap*MSH%corrL
            ovlp_Tot  = ovlp*dble(MSH%procPerDim-1) +2.0D0*(ovlp-MSH%xStep) !We consider an overlap -1 in each extreme
            where(MSH%procPerDim == 1) ovlp = 0.0D0
            where(MSH%procPerDim == 1) ovlp_Tot = 0.0D0
            where(ovlp_Tot <= 0.0D0) ovlp_Tot = 0.0D0
            call wLog("        OVERLAP = ")
            call wLog(ovlp)
            call wLog("        OVERLAP TOTAL = ")
            call wLog(ovlp_Tot)
            Novlp_Tot = (MSH%xMaxGlob - MSH%xMinGlob) - ovlp_Tot
            call wLog("        NON-OVERLAP TOTAL (rest) = ")
            call wLog(Novlp_Tot)
            where(Novlp_Tot < 0.0D0)
                Novlp_Tot = 0.0D0
            end where
            call roundToMultiple(Novlp_Tot, MSH%xStep*dble(MSH%procPerDim), up=.true.)
            call wLog("        NON-OVERLAP TOTAL (Round) = ")
            call wLog(Novlp_Tot)
            Novlp = Novlp_Tot/dble(MSH%procPerDim)

            call wLog("        NON-OVERLAP TOTAL = ")
            call wLog(Novlp_Tot)

            !Redefining Global Extremes
            delta = ovlp_Tot + Novlp_Tot
            xMaxGlob = MSH%xMinGlob + delta
            call wLog("        OUT MSH%xMinGlob (Round 2) = ")
            call wLog(MSH%xMinGlob)
            call wLog("        OUT MSH%xMaxGlob (Round 2) = ")
            call wLog(MSH%xMaxGlob)

            !Extent of a proc
            delta_Proc = (Novlp_Tot)/dble(MSH%procPerDim) + 2.0D0*ovlp - 2*MSH%xStep
            where(MSH%procPerDim == 1)  delta_Proc = (Novlp_Tot)/dble(MSH%procPerDim)
            stepProc   =  Novlp + ovlp
            call wLog("        delta_Proc = ")
            call wLog(delta_Proc)
            call wLog("        stepProc = ")
            call wLog(stepProc)

            procExtent = delta_Proc
            procStart  = stepProc*MSH%coords

        end if

    end subroutine set_global_extremes

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_local_bounding_box (MSH, procStart, procExtent, xMinBound, xMaxBound, &
                                       xNStep, xNTotal, origin)

        implicit none

        type(MESH), intent(in) :: MSH
        double precision, dimension(:), intent(in) :: procStart, procExtent
        !INPUT AND OUTPUT
        double precision, dimension(:), intent(out) :: xMinBound, xMaxBound
        integer         , dimension(:), intent(out) :: xNStep, origin
        integer(kind=8) , intent(out) :: xNTotal

        !FFTW
        integer(C_INTPTR_T) :: L, M, N
        integer(C_INTPTR_T) :: local_LastDim
        integer(C_INTPTR_T) :: local_j_offset
        integer(C_INTPTR_T) :: alloc_local
        integer, dimension(MSH%nDim) :: xNStepGlob

        xMinBound = procStart
        xMaxBound = procStart + procExtent

        !GLOBAL FFT
        if(.not. MSH%independent) then

            if(MSH%method == FFT) then

                xNStepGlob = find_xNStep(xMaxExt=(MSH%xMaxGlob - MSH%xMinGlob), xStep=MSH%xStep)
                call wLog("xNStepGlob = ")
                call wLog(xNStepGlob)

                call wLog("FFT Local Division")

                if(xNStepGlob(MSH%nDim) < MSH%nb_procs) stop("ERROR!! When using parallel FFT the last dimension should have at least 1 slice by proc")

                if(MSH%nDim == 2) then
                    L = xNStepGlob(1)
                    M = xNStepGlob(2)
                    alloc_local = fftw_mpi_local_size_2d(M, L, MSH%comm, &
                                                         local_LastDim, local_j_offset) !FOR MPI
                else if(MSH%nDim == 3) then
                    L = xNStepGlob(1)
                    M = xNStepGlob(2)
                    N = xNStepGlob(3)
                    alloc_local = fftw_mpi_local_size_3d(N, M, L, MSH%comm, &
                                                         local_LastDim, local_j_offset) !FOR MPI
                else
                    stop("Inside set_local_extremes no mesh division for FFT in this dimension")
                end if

                call wLog("local_LastDim = ")
                call wLog(int(local_LastDim))
                call wLog("local_j_offset = ")
                call wLog(int(local_j_offset))

                xMinBound = MSH%xMinGlob
                xMinBound(MSH%nDim) = MSH%xMinGlob(MSH%nDim) + dble(local_j_offset)*MSH%xStep(MSH%nDim)
                xMaxBound = MSH%xMaxGlob
                xMaxBound(MSH%nDim) = xMinBound(MSH%nDim) + dble(local_LastDim-1)*MSH%xStep(MSH%nDim)

            else
                !Taking out the repeated point on the positive extreme
                where(MSH%coords /= (MSH%procPerDim - 1)) xMaxBound = xMaxBound - MSH%xStep
            end if
        end if

        call wLog("        OUT xMinBound = ")
        call wLog(xMinBound)
        call wLog("        OUT xMaxBound = ")
        call wLog(xMaxBound)

        xNStep = find_xNStep(xMinBound, xMaxBound, MSH%xStep)
        origin = find_xNStep(MSH%xMinGlob, xMinBound , MSH%xStep)
        xNTotal = product(xNStep)

        call wLog("        OUT xNStep = ")
        call wLog(xNStep)
        call wLog("        OUT xNTotal = ")
        call wLog(int(xNTotal))
        call wLog("        OUT origin = ")
        call wLog(origin)

    end subroutine set_local_bounding_box

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_overlap_geometry (MSH, xMinInt, xMaxInt, xMinExt, xMaxExt, &
                                     xMaxNeigh, xMinNeigh, xOrNeigh)

        implicit none

        !INPUT
        type(MESH), intent(in) :: MSH

        !OUTPUT
        double precision, dimension(:)  , intent(out) :: xMinInt, xMaxInt, xMinExt, xMaxExt
        double precision, dimension(:,:), intent(out) :: xMaxNeigh, xMinNeigh, xOrNeigh

        !LOCAL VARIABLES
        double precision :: ovlpFraction = 1.0D0
        integer :: neighPos



        if(.not. MSH%independent) then
            call wLog("    ")
            call wLog("     GLOBAL GENERATION")
            call wLog(" ")
            call wLog("-> Copy Global Limits")
            call wLog(" ")
            call wLog(" IN MSH%xMinGlob")
            call wLog(MSH%xMinGlob)
            call wLog(" IN MSH%xMaxGlob")
            call wLog(MSH%xMaxGlob)

            xMinExt = MSH%xMinGlob
            xMaxExt = MSH%xMaxGlob
            xMinInt = xMinExt
            xMaxInt = xMaxExt

            call wLog(" OUT xMinInt")
            call wLog(xMinInt)
            call wLog(" OUT xMaxInt")
            call wLog(xMaxInt)
            call wLog(" OUT xMinExt")
            call wLog(xMinExt)
            call wLog(" OUT xMaxExt")
            call wLog(xMaxExt)

        else if(MSH%independent) then
            call wLog("    ")
            call wLog("     LOCAL GENERATION")
            call wLog(" ")
            call wLog("-> Redimensioning for Overlap")
            call wLog(" ")
            call wLog(" IN MSH%xMinBound")
            call wLog(MSH%xMinBound)
            call wLog(" IN MSH%xMaxBound")
            call wLog(MSH%xMaxBound)
            !Finding External Extremes
            call wLog(" ")
            call wLog("    Finding External Extremes")

            xMinExt = MSH%xMinBound
            xMaxExt = MSH%xMaxBound

            do neighPos = 1, 2*MSH%nDim

                if(MSH%neigh(neighPos) < 0) cycle

                where(MSH%neighShift(:,neighPos) < 0)
                    xMinExt = MSH%xMinBound - MSH%xStep
                elsewhere(MSH%neighShift(:,neighPos) > 0)
                    xMaxExt = MSH%xMaxBound + MSH%xStep
                end where

            end do

            call wLog(" OUT xMinExt")
            call wLog(xMinExt)
            call wLog(" OUT xMaxExt")
            call wLog(xMaxExt)


            !Redimensioning the internal part
            call wLog("    Redimensioning the internal part")
            xMinInt = xMinExt
            xMaxInt = xMaxExt

            do neighPos = 1, 2*MSH%nDim

                if(MSH%neigh(neighPos) < 0) cycle

                where(MSH%neighShift(:,neighPos) < 0)
                    xMinInt = xMinExt + MSH%overlap*MSH%corrL
                elsewhere(MSH%neighShift(:,neighPos) > 0)
                    xMaxInt = xMaxExt - MSH%overlap*MSH%corrL
                end where

            end do

            call wLog(" OUT xMinInt")
            call wLog(xMinInt)
            call wLog(" OUT xMaxInt")
            call wLog(xMaxInt)

            !Dimensioning overlapping area
            call wLog("    Dimensioning neighbours limits")

            do neighPos = 1, size(MSH%neigh)
                if(MSH%neigh(neighPos) < 0) cycle

                where(MSH%neighShift(:,neighPos) > 0)
                    !MSH%xMaxNeigh(:,neighPos) = MSH%xMaxExt - MSH%xStep - MSH%overlap*corrL*(1.0D0-ovlpFraction)
                    xMaxNeigh(:,neighPos) = MSH%xMaxBound
                    xMinNeigh(:,neighPos) = MSH%xMaxInt + MSH%xStep
                    xOrNeigh(:,neighPos)  = MSH%xMaxInt
                elsewhere(MSH%neighShift(:,neighPos) < 0)
                    xMaxNeigh(:,neighPos) = MSH%xMinInt - MSH%xStep
                    xMinNeigh(:,neighPos) = MSH%xMinBound
                    xOrNeigh(:,neighPos)  = MSH%xMinInt
                elsewhere
                    xMaxNeigh(:,neighPos) = MSH%xMaxBound
                    xMinNeigh(:,neighPos) = MSH%xMinBound
                    xOrNeigh(:,neighPos)  = MSH%xMinInt
                end where
            end do

            call show_MESHneigh(MSH, " ", onlyExisting = .true., forLog = .true.)
        end if

    end subroutine set_overlap_geometry

!    !---------------------------------------------------------------------------------
!    !---------------------------------------------------------------------------------
!    !---------------------------------------------------------------------------------
!    !---------------------------------------------------------------------------------
!    subroutine define_topography(RDF, MSH, coordList)
!
!        implicit none
!        !INPUT OUTPUT
!        type(RF)   :: RDF
!        type(MESH) :: MSH
!        !INPUT
!        double precision, dimension(:,:), target, intent(in), optional :: coordList
!        !LOCAL
!        integer :: code
!        integer, dimension (MSH%nDim) :: tempXNStep
!
!
!        if (MSH%meshMod == msh_UNV) then
!            RDF%xPoints => coordList
!
!            RDF%xNTotal = MSH%xNTotal
!            RDF%xMinExt = MSH%xMinGlob
!            RDF%xMaxExt = MSH%xMaxGlob
!            MSH%xMinExt = MSH%xMinGlob
!            MSH%xMaxExt = MSH%xMaxGlob
!            MSH%xNTotal = size(RDF%xPoints, 2)
!            RDF%xNTotal = MSH%xNTotal
!            MSH%xMinBound = MSH%xMinGlob
!            MSH%xMaxBound = MSH%xMaxGlob
!            call wLog("     MSH%xMinBound = ")
!            call wLog(MSH%xMinBound)
!            call wLog("     MSH%xMaxBound = ")
!            call wLog(MSH%xMaxBound)
!            call wLog(" ")
!
!        else if(RDF%independent) then
!
!            call wLog("-> define_generation_geometry")
!            call define_generation_geometry (MSH, RDF)
!            call wLog("-> Getting Global Matrix Reference")
!            call get_XPoints_globCoords(RDF, MSH)
!            call wLog("     RDF%origin = ")
!            call wLog(RDF%origin)
!            call wLog(" ")
!
!        else
!            call wLog("-> define_generation_geometry")
!            call define_generation_geometry (MSH, RDF)
!            tempXNStep = find_xNStep(MSH%xMinGlob, MSH%xMaxGlob, MSH%xStep)
!            call get_Permutation(MSH%xNInit, MSH%xMaxGlob, tempXNStep, &
!                                 MSH%xMinBound, MSH%xMinGlob,  &
!                                 snapExtremes = .true.)
!            call get_Permutation(MSH%xNEnd, MSH%xMaxGlob, tempXNStep, &
!                                 MSH%xMaxBound, MSH%xMinGlob,  &
!                                 snapExtremes = .true.);
!            call wLog("     MSH%xMinBound (Linear)= ")
!            call wLog(MSH%xMinBound)
!            call wLog("     MSH%xMaxBound (Linear) = ")
!            call wLog(MSH%xMaxBound)
!            call wLog(" ")
!            call wLog("-> Getting Global Matrix Reference")
!            call get_XPoints_globCoords(RDF, MSH)
!            call wLog("     RDF%origin = ")
!            call wLog(RDF%origin)
!            call wLog(" ")
!
!        end if
!
!    end subroutine define_topography
!
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine define_generation_geometry (MSH, RDF)
!
!        implicit none
!
!        !INPUT AND OUTPUT
!        type(MESH), intent(inout) :: MSH
!        type(RF)  , intent(inout) :: RDF
!
!        !LOCAL
!        double precision, dimension(MSH%nDim) :: delta, half, minVol
!        double precision, dimension(MSH%nDim) :: ovlp_Tot, Novlp_Tot, ovlp, Novlp
!        double precision, dimension(MSH%nDim) :: delta_Proc, stepProc
!        double precision, dimension(MSH%nDim) :: locNO, locO
!        integer, dimension(MSH%nDim) :: nPointsO, nPointsO_Ext, nPointsTot, nPointsLoc
!        integer, dimension(MSH%nDim) :: nPointsNO, nPointsNO_int, nPointsNO_ext
!        integer, dimension(MSH%nDim) :: pointInit, pointEnd, pointNb
!        integer, dimension(MSH%nDim) :: xNStepGlob
!        integer :: sliceSize
!        integer :: i
!
!        !FFTW
!        integer(C_INTPTR_T) :: L, M, N
!        integer(C_INTPTR_T) :: local_LastDim
!        integer(C_INTPTR_T) :: local_j_offset
!        integer(C_INTPTR_T) :: alloc_local
!
!        !GLOBAL MANIPULATIONS
!
!        !Defining MSH%xStep
!        call wLog(" ")
!        call wLog("     IN  RDF%corrL          = ")
!        call wLog(RDF%corrL)
!        call wLog("     IN  MSH%pointsPerCorrL = ")
!        call wLog(MSH%pointsPerCorrL)
!
!        MSH%xStep = RDF%corrL/dble(MSH%pointsPerCorrL)
!
!        call wLog("     OUT MSH%xStep          = ")
!        call wLog(MSH%xStep)
!        call wLog(" ")
!
!        !Rounding Global Extremes
!
!        call wLog(" ")
!        call wLog("        IN  MSH%xMinGlob = ")
!        call wLog(MSH%xMinGlob)
!        call wLog("        IN  MSH%xMaxGlob = ")
!        call wLog(MSH%xMaxGlob)
!        call wLog("            delta        = ")
!        call wLog(MSH%xMaxGlob - MSH%xMinGlob)
!
!        call roundToMultiple(MSH%xMinGlob, MSH%xStep, up=.false.)
!        call roundToMultiple(MSH%xMaxGlob, MSH%xStep, up=.true.)
!
!        call wLog("        OUT MSH%xMinGlob (Round 1) = ")
!        call wLog(MSH%xMinGlob)
!        call wLog("        OUT MSH%xMaxGlob (Round 1) = ")
!        call wLog(MSH%xMaxGlob)
!        call wLog("            delta        (Round 1) = ")
!        call wLog(MSH%xMaxGlob - MSH%xMinGlob)
!        call wLog(" ")
!
!        !GLOBAL CASE (Not INDEPENDENT)
!        if(.not. MSH%independent) then
!
!            !Global extremes
!            half  = (MSH%xMaxGlob + MSH%xMinGlob)/2.0D0
!            delta = MSH%xMaxGlob - MSH%xMinGlob
!            call roundToMultiple(delta, MSH%xStep*dble(MSH%procPerDim)*2.0D0, up=.true.)
!            MSH%xMinGlob = half - delta/2.0D0
!            MSH%xMaxGlob = half + delta/2.0D0
!
!            call wLog("        OUT MSH%xMinGlob (Round 2) = ")
!            call wLog(MSH%xMinGlob)
!            call wLog("        OUT MSH%xMaxGlob (Round 2) = ")
!            call wLog(MSH%xMaxGlob)
!            call wLog("            delta        (Round 2) = ")
!            call wLog(MSH%xMaxGlob - MSH%xMinGlob)
!            call wLog(" ")
!
!
!            xNStepGlob = find_xNStep(xMaxExt=(MSH%xMaxGlob - MSH%xMinGlob), xStep=MSH%xStep)
!
!            if(RDF%method == FFT) then
!
!                if(xNStepGlob(MSH%nDim) < MSH%nb_procs) then
!                   xNStepGlob(MSH%nDim) = MSH%nb_procs
!                   half  = (MSH%xMaxGlob + MSH%xMinGlob)/2.0D0
!                   delta(MSH%nDim) = xNStepGlob(MSH%nDim)
!                   call roundToMultiple(delta, MSH%xStep*dble(MSH%procPerDim)*2.0D0, up=.true.)
!                   MSH%xMinGlob = half - delta/2.0D0
!                   MSH%xMaxGlob = half + delta/2.0D0
!
!                   call wLog("        OUT MSH%xMinGlob (Round 3 - For FFT) = ")
!                   call wLog(MSH%xMinGlob)
!                   call wLog("        OUT MSH%xMaxGlob (Round 3 - For FFT) = ")
!                   call wLog(MSH%xMaxGlob)
!                   call wLog("            delta        (Round 3 - For FFT) = ")
!                   call wLog(MSH%xMaxGlob - MSH%xMinGlob)
!                   call wLog(" ")
!
!                   xNStepGlob = find_xNStep(xMaxExt=(MSH%xMaxGlob - MSH%xMinGlob), xStep=MSH%xStep)
!                end if
!
!                call wLog("xNStepGlob = ")
!                call wLog(xNStepGlob)
!
!                if(xNStepGlob(MSH%nDim) < MSH%nb_procs) stop("ERROR!! When using parallel FFT the last dimension should have at least 1 slice by proc")
!
!                if(RDF%nDim == 2) then
!                    L = xNStepGlob(1)
!                    M = xNStepGlob(2)
!                    alloc_local = fftw_mpi_local_size_2d(M, L, RDF%comm, &
!                                                         local_LastDim, local_j_offset) !FOR MPI
!                else if(RDF%nDim == 3) then
!                    L = xNStepGlob(1)
!                    M = xNStepGlob(2)
!                    N = xNStepGlob(3)
!                    alloc_local = fftw_mpi_local_size_3d(N, M, L, RDF%comm, &
!                                                         local_LastDim, local_j_offset) !FOR MPI
!                else
!                    stop("Inside define_generation_geometry no mesh division for FFT in this dimension")
!                end if
!
!                call wLog("local_LastDim = ")
!                call wLog(int(local_LastDim))
!                call wLog("local_j_offset = ")
!                call wLog(int(local_j_offset))
!
!                MSH%xNGlob = product(xNStepGlob)
!                MSH%xNInit = local_j_offset + 1
!                MSH%xNEnd  = MSH%xNInit + local_LastDim - 1
!
!                sliceSize = 1
!                if(MSH%nDim > 1) sliceSize = product(xNStepGlob(1:RDF%nDim -1))
!                MSH%xNInit = (MSH%xNInit - 1) * sliceSize + 1
!                MSH%xNEnd  =  MSH%xNEnd *sliceSize
!
!            else
!                MSH%xNGlob = product(xNStepGlob)
!
!                MSH%xNEnd = (MSH%rang + 1) * MSH%XNGlob/MSH%nb_procs + 1
!                MSH%xNInit  = MSH%rang * MSH%XNGlob/MSH%nb_procs + 1
!
!                if(MSH%rang == MSH%nb_procs-1) then
!                    MSH%xNEnd = MSH%XNGlob
!                else
!                    MSH%xNEnd = MSH%xNEnd - 1
!                end if
!            end if
!
!            MSH%xNStep = find_xNStep(MSH%xMinGlob, MSH%xMaxGlob, MSH%xStep)
!            MSH%xNTotal = MSH%xNEnd - MSH%xNInit + 1
!            RDF%xNStep  = MSH%xNStep
!            RDF%xNTotal = MSH%xNTotal
!
!            call wLog("        OUT MSH%xNInit = ")
!            call wLog(int(MSH%xNInit))
!            call wLog("        OUT MSH%xNEnd = ")
!            call wLog(int(MSH%xNEnd))
!            call wLog("        OUT MSH%xNTotal = ")
!            call wLog(int(MSH%xNTotal))
!            call wLog("        OUT RDF%xNTotal = ")
!            call wLog(int(RDF%xNTotal))
!            call wLog("        OUT MSH%xNTotal = ")
!            call wLog(int(MSH%xNTotal))
!            call wLog("        OUT MSH%xNStep = ")
!            call wLog(int(MSH%xNStep))
!
!        !INDEPENDENT CASE (For Localization)
!        else if(MSH%independent) then
!
!            !Rounding overlap
!            call wLog(" ")
!            call wLog("     IN   MSH%overlap =  ")
!            call wLog(MSH%overlap)
!
!            if(MSH%independent) then
!                MSH%overlap = ceiling(MSH%overlap*RDF%corrL/(2*MSH%xStep)) * 2*MSH%xStep/RDF%corrL
!            else
!                MSH%overlap = 0
!            end if
!
!            call roundToMultiple(MSH%overlap, MSH%xStep, up=.true.)
!
!            call wLog("     OUT  MSH%overlap =  ")
!            call wLog(MSH%overlap)
!            call wLog(" ")
!
!
!            !Size of the non-overlapping areas
!            ovlp      = MSH%overlap*RDF%corrL
!            ovlp_Tot  = ovlp*dble(MSH%procPerDim-1) +2.0D0*(ovlp-MSH%xStep) !We consider an overlap -1 in each extreme
!            where(MSH%procPerDim == 1) ovlp = 0.0D0
!            where(MSH%procPerDim == 1) ovlp_Tot = 0.0D0
!            where(ovlp_Tot <= 0.0D0) ovlp_Tot = 0.0D0
!            call wLog("        OVERLAP = ")
!            call wLog(ovlp)
!            call wLog("        OVERLAP TOTAL = ")
!            call wLog(ovlp_Tot)
!            Novlp_Tot = (MSH%xMaxGlob - MSH%xMinGlob) - ovlp_Tot
!            call wLog("        NON-OVERLAP TOTAL (rest) = ")
!            call wLog(Novlp_Tot)
!            where(Novlp_Tot < 0.0D0)
!                Novlp_Tot = 0.0D0
!            end where
!            call roundToMultiple(Novlp_Tot, MSH%xStep*dble(MSH%procPerDim), up=.true.)
!            call wLog("        NON-OVERLAP TOTAL (Round) = ")
!            call wLog(Novlp_Tot)
!            Novlp = Novlp_Tot/dble(MSH%procPerDim)
!
!            call wLog("        NON-OVERLAP TOTAL = ")
!            call wLog(Novlp_Tot)
!
!            !Redefining Global Extremes
!            delta = ovlp_Tot + Novlp_Tot
!            MSH%xMaxGlob = MSH%xMinGlob + delta
!            call wLog("        OUT MSH%xMinGlob (Round 2) = ")
!            call wLog(MSH%xMinGlob)
!            call wLog("        OUT MSH%xMaxGlob (Round 2) = ")
!            call wLog(MSH%xMaxGlob)
!
!            !Extent of a proc
!            delta_Proc = (Novlp_Tot)/dble(MSH%procPerDim) + 2.0D0*ovlp - 2*MSH%xStep
!            where(MSH%procPerDim == 1)  delta_Proc = (Novlp_Tot)/dble(MSH%procPerDim)
!            stepProc   =  Novlp + ovlp
!            call wLog("        delta_Proc = ")
!            call wLog(delta_Proc)
!            call wLog("        stepProc = ")
!            call wLog(stepProc)
!
!            !External Min and Max
!            MSH%xMinExt = stepProc*MSH%coords + MSH%xMinGlob
!            MSH%xMaxExt = MSH%xMinExt + delta_Proc
!            call wLog("        OUT MSH%xMinExt (1) = ")
!            call wLog(MSH%xMinExt)
!            call wLog("        OUT MSH%xMaxExt (1) = ")
!            call wLog(MSH%xMaxExt)
!
!            !Expanding where there are neighbours
!            do i=1, MSH%nDim
!                if(any(MSH%neighShift(i, :) == 1))  MSH%xMaxExt(i) = MSH%xMaxExt(i) + MSH%xStep(i)
!                if(any(MSH%neighShift(i, :) == -1)) MSH%xMinExt(i) = MSH%xMinExt(i) - MSH%xStep(i)
!            end do
!
!            call wLog("        OUT MSH%xMinExt = ")
!            call wLog(MSH%xMinExt)
!            call wLog("        OUT MSH%xMaxExt = ")
!            call wLog(MSH%xMaxExt)
!
!
!            call wLog("-> get_overlap_geometry")
!            if (RDF%method == FFT) then
!                call get_overlap_geometry (MSH, RDF%corrL, fullGen=.true.)
!            else
!                call get_overlap_geometry (MSH, RDF%corrL, fullGen=.false.)
!            end if
!
!            MSH%xNStep  = find_xNStep(MSH%xMinBound, MSH%xMaxBound, MSH%xStep)
!            MSH%xNTotal = product(MSH%xNStep)
!            RDF%xNTotal = MSH%xNTotal
!            RDF%xNStep  = MSH%xNStep
!
!            call wLog("        OUT RDF%xNTotal = ")
!            call wLog(int(RDF%xNTotal))
!            call wLog("        OUT MSH%xNTotal = ")
!            call wLog(int(MSH%xNTotal))
!            call wLog("        OUT MSH%xNStep = ")
!            call wLog(int(MSH%xNStep))
!            call wLog("        OUT RDF%xNStep = ")
!            call wLog(int(RDF%xNStep))
!
!        end if
!
!    end subroutine define_generation_geometry

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine get_xPoints_globCoords(MSH)
!        implicit none
!
!        !INPUT AND OUTPUT
!        type(RF)  , intent(inout) :: RDF
!        type(MESH), intent(inout) :: MSH
!
!        call wLog("MSH%xMinGlob = ")
!        call wLog(MSH%xMinGlob)
!        call wLog("MSH%xMinBound = ")
!        call wLog(MSH%xMinBound)
!
!        RDF%origin = find_xNStep(MSH%xMinGlob, MSH%xMinBound , MSH%xStep)
!
!    end subroutine get_xPoints_globCoords

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_communications_topology(MSH)

        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH
        !LOCAL
        integer :: code
        logical, dimension(MSH%nDim) :: periods

        periods(:) = .false.

        call wLog("-> set_procPerDim")
        call set_procPerDim (MSH, MSH%procPerDim)
        call wLog("     MSH%procPerDim")
        call wLog(MSH%procPerDim)
        call wLog("-> MPI_CART_CREATE")
        call MPI_CART_CREATE (MSH%comm, MSH%nDim, MSH%procPerDim, periods, .false., MSH%topComm, code)
        call wLog("-> MPI_CART_COORDS")
        call MPI_CART_COORDS (MSH%topComm, MSH%rang, MSH%nDim, MSH%coords, code)
        call wLog("     MSH%coords")
        call wLog(MSH%coords)

        if(MSH%independent) then
            call wLog("-> set_neighbours")
            call set_neighbours (MSH)
            call wLog("-> get_NeighbourCriteria")
            call get_NeighbourCriteria (MSH)
        end if

    end subroutine set_communications_topology
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

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_procPerDim (MSH, procPerDim)
        implicit none

        !INPUT
        type(MESH), intent(in) :: MSH

        !OUTPUT
        integer, dimension(:) :: procPerDim

        !LOCAL VARIABLES
        integer :: i, j;
        double  precision :: procRootDim, logProc2;

        if(size(procPerDim)/=MSH%nDim) then
            write(*,*) "Error inside 'set_procPerDim', dimensions are not compatible"
            write(*,*) "size(procPerDim) = ", size(procPerDim)
            write(*,*) "MSH%nDim         = ", MSH%nDim
            stop(" ")
        end if

        if(MSH%method == FFT .and. (.not. MSH%independent)) then

            procPerDim(:) = 1
            !procPerDim(MSH%nDim) = MSH%nb_procs
            procPerDim(1) = MSH%nb_procs

        else

            procRootDim = dble(MSH%nb_procs)**(1/dble(MSH%nDim))
            logProc2   = log(dble(MSH%nb_procs))/log(2.0D0)

            if (areEqual(procRootDim, dble(nint(procRootDim)))) then
                call wLog("    Exact Division")
                procPerDim(:) = nint(dble(MSH%nb_procs)**(1.0d0/MSH%nDim))

            else if(areEqual(logProc2, dble(nint(logProc2)))) then
                call wLog("    Power of two")

                procPerDim(:) = 1
                if(MSH%nb_procs /= 1) then
                    do j = 1, nint(logProc2)
                        i = cyclicMod(j, MSH%nDim)
                        procPerDim(i) = procPerDim(i)*2
                    end do
                end if
            else if(.not. MSH%independent) then
                call wLog("    1D Division")
                procPerDim(:) = 1
                procPerDim(MSH%nDim) = MSH%nb_procs
            else
                stop ("ERROR, no mesh division algorithm for this number of procs")
            end if

        end if

        call wLog("    procPerDim = ")
        call wLog(procPerDim)

    end subroutine set_procPerDim

end module topography_RF
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
