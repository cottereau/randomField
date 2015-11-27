module calls_RF

    use displayCarvalhol
    use write_Log_File
    use math_RF
    use constants_RF
    use mpi
    use writeResultFile_RF
    use type_RF
    use localization_RF
    use type_MESH
    use common_variables_RF
    use randomFieldND
    use mesh_RF

    implicit none

    interface createRandomField
       module procedure create_RF_Unstruct_Init
    end interface createRandomField

contains

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine create_RF_Unstruct_noInit (xPoints, corrL, corrMod, Nmc,   &
!                                          randField, method, seedStart,   &
!                                          margiFirst, fieldAvg, fieldVar, &
!                                          comm, rang, nb_procs, calculate, MSH)
!        !INPUT
!        double precision, dimension(1:, 1:), intent(in), target :: xPoints;
!        double precision, dimension(1:)    , intent(in) :: corrL;
!        integer                            , intent(in) :: corrMod;
!        integer                            , intent(in) :: Nmc;
!        integer                            , intent(in) :: method
!        integer                            , intent(in) :: seedStart
!        integer                            , intent(in) :: margiFirst;
!        double precision                   , intent(in) :: fieldAvg
!        double precision                   , intent(in) :: fieldVar;
!        integer                            , intent(in) :: comm, rang, nb_procs
!        logical, dimension(1:), optional   , intent(in) :: calculate
!        type(MESH), intent(inout) :: MSH
!
!        !OUTPUT
!        double precision, dimension(:, :), intent(out), target :: randField;
!
!        !LOCAL
!        type(RF) :: RDF
!
!        write(*,*) "Inside create_RF_Unstruct_noInit"
!
!        !Initializing RF
!        call init_RF(RDF, size(corrL), Nmc, comm, rang, nb_procs)
!        RDF%xPoints   => xPoints
!        RDF%randField => randField
!        RDF%xNTotal    = size(RDF%xPoints, 2)
!        RDF%corrL      = corrL
!        RDF%corrMod    = corrMod
!        RDF%Nmc        = Nmc
!        RDF%method     = method
!        RDF%seedStart  = seedStart
!        RDF%margiFirst = margiFirst
!        RDF%fieldAvg   = fieldAvg
!        RDF%fieldVar   = fieldVar
!        if(present(calculate)) RDF%calculate  = calculate
!
!        call create_RF_Unstruct_Init(RDF, MSH)
!
!    end subroutine create_RF_Unstruct_noInit

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine create_RF_Unstruct_Init (RDF, MSH)
        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        type(MESH), intent(inout) :: MSH
        !LOCAL
        logical, dimension(:), allocatable :: effectCalc;

        if(RDF%rang == 0) write(*,*) "Inside create_RF_Unstruct_Init"

        !Generating standard Gaussian Field
        call gen_Std_Gauss(RDF, MSH)

    end subroutine create_RF_Unstruct_Init

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss (RDF, MSH)
        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        type(MESH), intent(inout) :: MSH

        !LOCAL VARIABLES
        integer :: i;
        integer :: minIndexNeigh, maxIndexNeigh
        !logical, dimension(size(MSH%neigh)) :: considerNeighbour
        !integer, dimension(16) :: testVec
        integer :: partitionType = 1
        integer, dimension(RDF%nDim) :: minPos, maxPos

        !testVec = [(i, i = 1, 16)]

        !Normalization
        call wLog(" ")
        call wLog("->Normalizing Coordinates")
        call wLog(" ")
        do i = 1, RDF%nDim
            RDF%xPoints(i,:)   = RDF%xPoints(i,:)/RDF%corrL(i)
            MSH%xStep(i)       = MSH%xStep(i) /RDF%corrL(i)
            MSH%xMinInt(i)     = MSH%xMinInt(i)/RDF%corrL(i)
            MSH%xMaxInt(i)     = MSH%xMaxInt(i)/RDF%corrL(i)
            MSH%xMinExt(i)     = MSH%xMinExt(i)/RDF%corrL(i)
            MSH%xMaxExt(i)     = MSH%xMaxExt(i)/RDF%corrL(i)
            MSH%xMinGlob(i)    = MSH%xMinGlob(i)/RDF%corrL(i)
            MSH%xMaxGlob(i)    = MSH%xMaxGlob(i)/RDF%corrL(i)
            MSH%xMaxNeigh(i,:) = MSH%xMaxNeigh(i,:)/RDF%corrL(i)
            MSH%xMinNeigh(i,:) = MSH%xMinNeigh(i,:)/RDF%corrL(i)
            MSH%xMaxBound(i)   = MSH%xMaxBound(i)/RDF%corrL(i)
            MSH%xMinBound(i)   = MSH%xMinBound(i)/RDF%corrL(i)
            MSH%xOrNeigh(i,:)  = MSH%xOrNeigh(i,:)/RDF%corrL(i)
        end do

        if(RDF%independent) then
            !RDF%xRange = MSH%xMaxExt - MSH%xMinExt !Delta max in between two wave numbers to avoid periodicity
            RDF%xRange = MSH%xMaxBound - MSH%xMinBound !Delta max in between two wave numbers to avoid periodicity
        else
            RDF%xRange = MSH%xMaxGlob - MSH%xMinGlob !Delta max in between two wave numbers to avoid periodicity
        end if

        !Generating Standard Gaussian Field
        call wLog("")
        call wLog("GENERATING RANDOM FIELDS")
        call wLog("-------------------------------")
        if(RDF%rang == 0) write(*,*)"GENERATING RANDOM FIELDS"
        if(RDF%rang == 0) write(*,*) "-------------------------------"
        call wLog("")

        select case (RDF%method)
            case(ISOTROPIC)
                call wLog(" ISOTROPIC")
                if(RDF%rang == 0) write(*,*)"ISOTROPIC"
                call gen_Std_Gauss_Isotropic(RDF, MSH)
            case(SHINOZUKA)
                call wLog(" SHINOZUKA")
                if(RDF%rang == 0) write(*,*)"SHINOZUKA"
                call gen_Std_Gauss_Shinozuka(RDF, MSH)
            case(RANDOMIZATION)
                call wLog(" RANDOMIZATION")
                if(RDF%rang == 0) write(*,*)"RANDOMIZATION"
                call gen_Std_Gauss_Randomization(RDF, MSH)
            case(FFT)
                call wLog(" FFT")
                if(RDF%rang == 0) write(*,*)"FFT"
                call gen_Std_Gauss_FFT(RDF, MSH)
        end select

        !RDF%randField = 1.0 ! For Tests

        call wLog("minval(RDF%randField,1) =")
        call wLog(minval(RDF%randField,1))
        call wLog("maxval(RDF%randField,1) =")
        call wLog(maxval(RDF%randField,1))

        if(RDF%independent .and. RDF%nb_procs > 1) then
            call wLog("")
            call wLog("GENERATING OVERLAP")
            call wLog("-------------------------------")
            if(RDF%rang == 0) write(*,*)"GENERATING OVERLAP"
            if(RDF%rang == 0) write(*,*) "-------------------------------"
            call wLog("")

            !RDF%randField = 1.0 ! For Tests
            if(RDF%rang == 0) write(*,*) "    ->Applying Weighting Functions"
            call wLog("    ->Applying Weighting Functions on Field")
            call applyWeightingFunctions_OnMatrix(RDF, MSH, partitionType)
            if(RDF%rang == 0) write(*,*) "    ->addNeighboursFields"
            call wLog("    ->addNeighboursFields")
            call addNeighboursFields(RDF, MSH)

!            ! START For Tests
!            minPos = nint((MSH%xMinInt-MSH%xMinExt)/MSH%xStep) + 1
!            maxPos = nint((MSH%xMaxInt-MSH%xMinExt)/MSH%xStep)
!            call wLog("minPos = ")
!            call wLog(minPos)
!            call wLog("maxPos = ")
!            call wLog(maxPos)
!            if(RDF%nDim == 2) RDF%RF_2D(minPos(1):maxPos(1), &
!                                        minPos(2):maxPos(2)) = 0
!            if(RDF%nDim == 3) RDF%RF_3D(minPos(1):maxPos(1), &
!                                        minPos(2):maxPos(2), &
!                                        minPos(3):maxPos(3)) = 0
!            ! END For Tests

        end if

        !Reverting Normalization
        call wLog(" ")
        call wLog("->Reverting Normalization")
        do i = 1, RDF%nDim
            RDF%xPoints(i,:)   = RDF%xPoints(i,:)*RDF%corrL(i)
            RDF%xRange(i)      = RDF%xRange(i)*RDF%corrL(i)
            MSH%xStep(i)       = MSH%xStep(i)*RDF%corrL(i)
            MSH%xMinInt(i)     = MSH%xMinInt(i)*RDF%corrL(i)
            MSH%xMaxInt(i)     = MSH%xMaxInt(i)*RDF%corrL(i)
            MSH%xMinExt(i)     = MSH%xMinExt(i)*RDF%corrL(i)
            MSH%xMaxExt(i)     = MSH%xMaxExt(i)*RDF%corrL(i)
            MSH%xMinGlob(i)    = MSH%xMinGlob(i)*RDF%corrL(i)
            MSH%xMaxGlob(i)    = MSH%xMaxGlob(i)*RDF%corrL(i)
            MSH%xMaxNeigh(i,:) = MSH%xMaxNeigh(i,:)*RDF%corrL(i)
            MSH%xMinNeigh(i,:) = MSH%xMinNeigh(i,:)*RDF%corrL(i)
            MSH%xMaxBound(i)   = MSH%xMaxBound(i)*RDF%corrL(i)
            MSH%xMinBound(i)   = MSH%xMinBound(i)*RDF%corrL(i)
            MSH%xOrNeigh(i,:)  = MSH%xOrNeigh(i,:)*RDF%corrL(i)
        end do

        !RDF%randField = RDF%rang ! For Tests

    end subroutine gen_Std_Gauss

end module calls_RF
