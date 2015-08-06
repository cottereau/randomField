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

    !use spectra_RF

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
        character (len=*)                  , intent(in) :: corrMod;
        integer                            , intent(in) :: Nmc;
        integer                            , intent(in) :: method
        integer                            , intent(in) :: seedStart
        character (len=*)                  , intent(in) :: margiFirst;
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
        RDF%xMinBound = MSH%xMinBound
        RDF%xMaxBound = MSH%xMaxBound

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
        logical, dimension(size(MSH%neigh)) :: considerNeighbour

        !Normalization
        do i = 1, RDF%nDim
            RDF%xPoints(i,:) = RDF%xPoints(i,:)/RDF%corrL(i)
            RDF%xMinGlob(i)  = RDF%xMinGlob(i)/RDF%corrL(i)
            RDF%xMaxGlob(i)  = RDF%xMaxGlob(i)/RDF%corrL(i)
            RDF%xMinBound(i) = RDF%xMinBound(i)/RDF%corrL(i)
            RDF%xMaxBound(i) = RDF%xMaxBound(i)/RDF%corrL(i)
            RDF%xMinGlob(i)  = RDF%xMinGlob(i)/RDF%corrL(i)
            RDF%xMaxGlob(i)  = RDF%xMaxGlob(i)/RDF%corrL(i)
            MSH%xMaxNeigh(i,:) = MSH%xMaxNeigh(i,:)/RDF%corrL(i)
            MSH%xMinNeigh(i,:) = MSH%xMinNeigh(i,:)/RDF%corrL(i)
        end do

        !Generating Standard Gaussian Field

        write(get_fileId(),*) ""
        write(get_fileId(),*) "GENERATING INTERNAL RANDOM FIELD"
        write(get_fileId(),*) "-------------------------------"
        write(get_fileId(),*) ""

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

        if(RDF%independent) then
            !Communicating borders to neighbours
            write(get_fileId(),*) ""
            write(get_fileId(),*) "GENERATING BORDER RANDOM FIELDS"
            write(get_fileId(),*) "-------------------------------"
            write(get_fileId(),*) ""
            write(get_fileId(),*) "->Discovering neighbours seed"
            call get_neighbours_info(RDF, MSH)
            write(get_fileId(),*) "Creating Overlaps"
            call getNeighIndexRange(MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)
            call applyWeightingFunctions(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)
            call takeNeighboursContribution(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)
            call normalizeOverlap(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)

        end if

        !Reverting Normalization
        write(get_fileId(),*) " "
        write(get_fileId(),*) "->Reverting Normalization"
        do i = 1, RDF%nDim
            RDF%xPoints(i,:) = RDF%xPoints(i,:)*RDF%corrL(i)
            RDF%xMinGlob(i)  = RDF%xMinGlob(i)*RDF%corrL(i)
            RDF%xMaxGlob(i)  = RDF%xMaxGlob(i)*RDF%corrL(i)
            RDF%xMinBound(i) = RDF%xMinBound(i)*RDF%corrL(i)
            RDF%xMaxBound(i) = RDF%xMaxBound(i)*RDF%corrL(i)
            MSH%xMaxNeigh(i,:) = MSH%xMaxNeigh(i,:)*RDF%corrL(i)
            MSH%xMinNeigh(i,:) = MSH%xMinNeigh(i,:)*RDF%corrL(i)
        end do

    end subroutine gen_Std_Gauss

end module calls_RF
