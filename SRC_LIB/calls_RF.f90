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
       module procedure create_RF_Unstruct_noInit,   &
                        create_RF_Unstruct_Init
    end interface createRandomField

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine create_RF_Unstruct_noInit (xPoints, corrL, corrMod, Nmc,   &
                                          randField, method, seedStart,   &
                                          margiFirst, fieldAvg, fieldVar, &
                                          comm, rang, nb_procs, calculate, MSH)
        !INPUT
        double precision, dimension(1:, 1:), intent(in), target :: xPoints;
        double precision, dimension(1:)    , intent(in) :: corrL;
        integer                            , intent(in) :: corrMod;
        integer                            , intent(in) :: Nmc;
        integer                            , intent(in) :: method
        integer                            , intent(in) :: seedStart
        integer                            , intent(in) :: margiFirst;
        double precision                   , intent(in) :: fieldAvg
        double precision                   , intent(in) :: fieldVar;
        integer                            , intent(in) :: comm, rang, nb_procs
        logical, dimension(1:), optional   , intent(in) :: calculate
        type(MESH), intent(inout) :: MSH

        !OUTPUT
        double precision, dimension(:, :), intent(out), target :: randField;

        !LOCAL
        type(RF) :: RDF

        write(*,*) "Inside create_RF_Unstruct_noInit"

        !Initializing RF
        call init_RF(RDF, size(corrL), Nmc, comm, rang, nb_procs)
        RDF%xPoints   => xPoints
        RDF%randField => randField
        RDF%xNTotal    = size(RDF%xPoints, 2)
        RDF%corrL      = corrL
        RDF%corrMod    = corrMod
        RDF%Nmc        = Nmc
        RDF%method     = method
        RDF%seedStart  = seedStart
        RDF%margiFirst = margiFirst
        RDF%fieldAvg   = fieldAvg
        RDF%fieldVar   = fieldVar
        if(present(calculate)) RDF%calculate  = calculate

        call create_RF_Unstruct_Init(RDF, MSH)

    end subroutine create_RF_Unstruct_noInit

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

        !Discovering Global Extremes
        RDF%xMinGlob  = MSH%xMinGlob
        RDF%xMaxGlob  = MSH%xMaxGlob
        RDF%xMinExt = MSH%xMinExt
        RDF%xMaxExt = MSH%xMaxExt

        !Getting Mesh Information
        RDF%xNStep = MSH%xNStep

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

        !testVec = [(i, i = 1, 16)]

        !Normalization
        call wLog(" ")
        call wLog("->Normalizing Coordinates")
        call wLog(" ")
        do i = 1, RDF%nDim
            RDF%xPoints(i,:) = RDF%xPoints(i,:)/RDF%corrL(i)
            RDF%xMinGlob(i)  = RDF%xMinGlob(i)/RDF%corrL(i)
            RDF%xMaxGlob(i)  = RDF%xMaxGlob(i)/RDF%corrL(i)
            RDF%xMinExt(i) = RDF%xMinExt(i)/RDF%corrL(i)
            RDF%xMaxExt(i) = RDF%xMaxExt(i)/RDF%corrL(i)
            RDF%xMinGlob(i)  = RDF%xMinGlob(i)/RDF%corrL(i)
            RDF%xMaxGlob(i)  = RDF%xMaxGlob(i)/RDF%corrL(i)

            MSH%xMaxNeigh(i,:) = MSH%xMaxNeigh(i,:)/RDF%corrL(i)
            MSH%xMinNeigh(i,:) = MSH%xMinNeigh(i,:)/RDF%corrL(i)
            MSH%xMaxBound(i)   = MSH%xMaxBound(i)/RDF%corrL(i)
            MSH%xMinBound(i)   = MSH%xMinBound(i)/RDF%corrL(i)
            MSH%xOrNeigh(i,:)  = MSH%xOrNeigh(i,:)/RDF%corrL(i)
        end do

        !Generating Standard Gaussian Field

        call wLog("")
        call wLog("GENERATING INTERNAL RANDOM FIELD")
        call wLog("-------------------------------")
        if(RDF%rang == 0) write(*,*)"GENERATING INTERNAL RANDOM FIELD"
        if(RDF%rang == 0) write(*,*) "-------------------------------"
        call wLog("")

        select case (RDF%method)
            case(ISOTROPIC)
                call gen_Std_Gauss_Isotropic(RDF)
            case(SHINOZUKA)
                call gen_Std_Gauss_Shinozuka(RDF)
            case(RANDOMIZATION)
                call gen_Std_Gauss_Randomization(RDF)
            case(FFT)
                call gen_Std_Gauss_FFT(RDF)
        end select

        !RDF%randField = 0.0 ! For Tests

        if(RDF%independent) then
            !Communicating borders to neighbours
            call wLog("")
            call wLog("GENERATING BORDER RANDOM FIELDS")
            call wLog("-------------------------------")
            if(RDF%rang == 0) write(*,*)"GENERATING BORDER RANDOM FIELDS"
            if(RDF%rang == 0) write(*,*) "-------------------------------"
            call wLog("")
            call wLog("    ->Discovering neighbours seed")
            call get_neighbours_info(RDF, MSH)
            call wLog("    ->Discovering neighbours index")
            call getNeighIndexRange(MSH, minIndexNeigh, maxIndexNeigh)
            call wLog("    ->Applying Weighting Functions")
            call applyWeightingFunctions(RDF, MSH, minIndexNeigh, maxIndexNeigh, partitionType)
            call wLog("    ->Adding Neighbours Contribution")
            call takeNeighboursContribution(RDF, MSH, minIndexNeigh, maxIndexNeigh, partitionType)
        end if

        !Reverting Normalization
        call wLog(" ")
        call wLog("->Reverting Normalization")
        do i = 1, RDF%nDim
            RDF%xPoints(i,:) = RDF%xPoints(i,:)*RDF%corrL(i)
            RDF%xMinGlob(i)  = RDF%xMinGlob(i)*RDF%corrL(i)
            RDF%xMaxGlob(i)  = RDF%xMaxGlob(i)*RDF%corrL(i)
            RDF%xMinExt(i) = RDF%xMinExt(i)*RDF%corrL(i)
            RDF%xMaxExt(i) = RDF%xMaxExt(i)*RDF%corrL(i)

            MSH%xMaxNeigh(i,:) = MSH%xMaxNeigh(i,:)*RDF%corrL(i)
            MSH%xMinNeigh(i,:) = MSH%xMinNeigh(i,:)*RDF%corrL(i)
            MSH%xMaxBound(i)   = MSH%xMaxBound(i)*RDF%corrL(i)
            MSH%xMinBound(i)   = MSH%xMinBound(i)*RDF%corrL(i)
            MSH%xOrNeigh(i,:)  = MSH%xOrNeigh(i,:)*RDF%corrL(i)
        end do

        RDF%randField = RDF%rang ! For Tests

    end subroutine gen_Std_Gauss

end module calls_RF
