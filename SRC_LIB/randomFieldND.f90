module randomFieldND

    use displayCarvalhol
    use spectra_RF
    use math_RF
    use constants_RF
    use mesh_RF
    use mpi
    use write_Log_File
    use type_RF
    use type_MESH
    use common_variables_RF
    use writeResultFile_RF

    implicit none
    include 'fftw3.f'
    !WARNING before this line we have include 'fftw3.f'
    !use blas

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
        double precision, dimension(:), allocatable :: normFactor

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
            allocate(normFactor(minIndexNeigh:maxIndexNeigh))
            call findNormalizationFactors(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour, normFactor)
            call applyWeightingFunctions(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)
            call takeNeighboursContribution(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour, normFactor)

            !RDF%randField(:, 1) = 0
            !RDF%randField(minIndexNeigh:maxIndexNeigh, 1) = normFactor
            !call create_RF_overlap(RDF, MSH)

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

        if(allocated(normFactor)) deallocate(normFactor)

    end subroutine gen_Std_Gauss

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

        implicit none

        !INPUT
        character (len=*), intent(in) :: margiFirst;
        double precision , intent(in) :: fieldAvg, fieldVar;

        !OUTPUT (IN)
        double precision, dimension(1:, 1:), intent(inout) :: randField;

        !LOCAL VARIABLES
        double precision :: normalVar, normalAvg
        integer          :: error, code, i

        select case (margiFirst)
        case("gaussian")
            normalVar = fieldVar
            normalAvg = fieldAvg
        case("lognormal")
            if(fieldAvg <= 0) then
                write(*,*) ""
                write(*,*) "ERROR - when using lognormal fieldAvg should be a positive number"
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            normalVar = log(1 + fieldVar/(fieldAvg**2))
            normalAvg = log(fieldAvg) - normalVar/2
        end select

        randField(:,:) = randField(:,:) * sqrt(normalVar) &
            + normalAvg;

        if (margiFirst == "lognormal") then
            randField(:,:) = exp(randField(:,:))
        end if

    end subroutine multiVariateTransformation

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss_Shinozuka(RDF, randomK_in)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        logical, intent(in), optional ::randomK_in

        !LOCAL
        double precision, dimension(:, :), allocatable :: phiK, kVec;
        double precision, dimension(:)   , allocatable :: dgemm_mult;
        double precision, dimension(:,:) , allocatable :: k_x_phi, kSign;
        double precision :: ampMult
        integer :: testIndex = 65
        integer :: n, i, j, m
        logical :: randomK

        !write(get_fileId(),*) "Inside Shinozuka"

        randomK = .false.
        if(present(randomK_in)) randomK = randomK_in
        call init_random_seed(RDF%seed)

        !write(get_fileId(),*) "Defining kPoints and SkVec"
        call set_kPoints(RDF)
        call set_SkVec(RDF)

        !write(get_fileId(),*) "Calculating Fields"
        allocate(k_x_phi (RDF%xNTotal, 1))
        allocate(kSign (2**(RDF%nDim-1), RDF%nDim));
        allocate(kVec(RDF%nDim, 1))

        call set_kSign(kSign) !Set the sign permutations for kVec

        allocate(phiK (RDF%kNTotal, size(kSign,1)));

        if(randomK) then
            write(get_fileId(),*) "-----Shinozuka, k random-----"
            ampMult = 2.0d0*sqrt(1/(RDF%kNTotal*(2.0d0*PI)**(dble(RDF%nDim))))
        else
            write(get_fileId(),*) "-----Shinozuka, k discrete-----"
            ampMult = 2.0d0*sqrt(product(RDF%kDelta)/((2.0d0*PI)**(dble(RDF%nDim))))
        end if

        write(get_fileId(),*) "     kNStep  = ", RDF%kNStep
        write(get_fileId(),*) "     kNTotal = ", RDF%kNTotal
        write(get_fileId(),*) "     xNTotal = ", size(RDF%xPoints, 2)

        RDF%randField(:,:) = 0.0d0;

        do n = 1, RDF%Nmc
            write(get_fileId(),*) "  --Generating Field Number ", n

            if(.not. RDF%calculate(n)) cycle

            call random_number(phiK(:,:))
            phiK(:,:) = 2.0D0*pi*phiK(:,:)

            write(get_fileId(),*) "     First PhiK = ", phiK(1,1)
            write(get_fileId(),*) "     Last PhiK  = ", phiK(size(phiK,1), size(phiK,2))

            !Loop on k sign
            do m = 1, size(kSign,1)

                !Changing kPoints Sign
                do i = 1, RDF%nDim
                    if(kSign(m, i) == -1) RDF%kPoints(i,:) = -RDF%kPoints(i,:)
                end do

                !Loop on k
                do j = 1, RDF%kNTotal

                    call DGEMM_simple(RDF%xPoints, RDF%kPoints(:,j:j), k_x_phi(:,:), "T", "N") !x*k

                    RDF%randField(:,n) = sqrt(RDF%SkVec(j)) * &
                                         cos(k_x_phi(:,1) + phiK(j, m)) &
                                         + RDF%randField(:,n)

                end do !END Loop on k

                !Reverting kPoints Sign
                do i = 1, RDF%nDim
                    if(kSign(m, i) == -1) RDF%kPoints(i,:) = -RDF%kPoints(i,:)
                end do

            end do !END Loop on k sign

        end do !END Loop on Nmc

        RDF%randField(:,:) = ampMult * RDF%randField(:,:);

        if(allocated(dgemm_mult))       deallocate(dgemm_mult)
        if(allocated(phiK))             deallocate(phiK);
        if(allocated(k_x_phi))          deallocate(k_x_phi)
        if(allocated(kSign))            deallocate(kSign)

    end subroutine gen_Std_Gauss_Shinozuka

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss_Randomization(RDF)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF

        !LOCAL
        logical :: randomK;

        write(get_fileId(),*) "-----Randomization-----"

        randomK = .true.

        call gen_Std_Gauss_Shinozuka(RDF, randomK)

    end subroutine gen_Std_Gauss_Randomization


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine gen_Std_Gauss_Isotropic(RDF)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF

        !LOCAL
        double precision, dimension(:)  , allocatable :: gammaN, phiN, thetaN, psiN;
        double precision, dimension(:)  , allocatable :: rVec;
        logical         , dimension(:)  , allocatable :: effectCalc;
        double precision :: rMax, Sk
        double precision, dimension(1) :: rMaxVec
        integer          :: i, j, k, m;
        integer          :: rNTotal;
        double precision :: step, rDelta;
        double precision, dimension(:), allocatable :: dgemm_mult;

        !Allocating
        allocate(rVec (RDF%nDim));
        allocate(dgemm_mult(RDF%xNTotal))

        write(get_fileId(),*) "-----Inside Isotropic-----"

        !r Definition
        call set_kMaxND(RDF%corrMod, rMaxVec)
        call set_kPoints(RDF)
        rMax = rMaxVec(1)
        rDelta  = maxval(RDF%kDelta(:))/5.0D0 !Delta min in between two wave numbers to avoid periodicity
        rNTotal = ceiling(rMax/rDelta) + 1;

        !write(get_fileId(),*) "rMax = ", rMax
        !write(get_fileId(),*) "rDelta = ", rDelta
        !write(get_fileId(),*) "rNTotal = ", rNTotal
        !write(get_fileId(),*) "RDF%calculate = ", RDF%calculate

        !Generating random field samples
        step      = rMax/dble(rNTotal)
        RDF%randField(:,:) = 0.0D0;

        call init_random_seed(RDF%seed)

        if (RDF%nDim == 2) then
            allocate(psiN   (rNTotal));
            allocate(thetaN (rNTotal));
            allocate(gammaN (rNTotal));
            do k = 1, RDF%Nmc
                if(RDF%calculate(k)) then
                    call random_number(psiN(:))
                    call random_number(thetaN(:))
                    call random_number(gammaN(:))
                    psiN   = 2d0*pi*psiN
                    thetaN = 2d0*pi*psiN
                    gammaN = 2d0*pi*gammaN

                    do j = 1, rNTotal
                        rVec           = [cos(thetaN(j)) * (j-1)*step, &
                                          sin(thetaN(j)) * (j-1)*step]
                        Sk             = get_SpectrumND([(j-1)*step], RDF%corrMod); !Obs, here Sk is a scalar
                        call DGEMM ( "T", "N", RDF%xNTotal, 1, RDF%nDim, &
                            1.0d0, RDF%xPoints, RDF%nDim, rVec, RDF%nDim, 0.0d0, dgemm_mult, RDF%xNTotal)

                        RDF%randField(:,k) = sqrt(Sk*(j-1)*(dble(step**2))) * gammaN(j) &
                                            * cos(                           &
                                            dgemm_mult                &
                                            + psiN(j)                 &
                                            )                          &
                                            + RDF%randField(:,k)
                    end do
                else
                    RDF%randField(:,k) = 0.0
                end if
            end do

        else if (RDF%nDim == 3) then
            !write(*,*) "nDim = 3 !!!"
            !write(*,*) "k = ",k;
            allocate(psiN   (rNTotal));
            allocate(thetaN (rNTotal));
            allocate(phiN   (rNTotal));
            allocate(gammaN (rNTotal));
            do k = 1, RDF%Nmc
                if(RDF%calculate(k)) then
                    !write(*,*) "k = ",k;
                    !write(*,*) "rNTotal = ",rNTotal;
                    call random_number(phiN(:))
                    call random_number(thetaN(:))
                    call random_number(gammaN(:))
                    call random_number(psiN(:))

                    psiN   = 2*pi*psiN
                    thetaN = 2*pi*psiN
                    phiN   = pi*phiN
                    gammaN = sqrt(12.0)*(gammaN -0.5d0)

                    do j = 1, rNTotal
                        !write(*,*) "j = ", j
                        rVec           = [cos(thetaN(j))*sin(phiN(j)) * (j-1)*step, &
                                          sin(thetaN(j))*sin(phiN(j)) * (j-1)*step, &
                                          cos(phiN(j))                * (j-1)*step]
                        Sk             = get_SpectrumND([(j-1)*step], RDF%corrMod);
                        call DGEMM ( "T", "N", RDF%xNTotal, 1, RDF%nDim, &
                                    1.0d0, RDF%xPoints, RDF%nDim, rVec, RDF%nDim, 0.0d0, dgemm_mult, RDF%xNTotal)
                        RDF%randField(:,k) = sqrt(Sk*sin(phiN(j))*step*((j-1)*step)**2) * gammaN(j) &
                                              * cos(                                             &
                                              dgemm_mult                                   &
                                              + psiN(j)                                    &
                                              )                                            &
                                              + RDF%randField(:,k)
                    end do
                else
                    RDF%randField(:,k) = 0.0
                end if
            end do

        else
            write(*,*) "ERROR The number of dimensions is not accepted in this method (Isotropic)";
            write(*,*) "RDF%nDim = ", RDF%nDim;
            stop
        end if

        !if(rang == 0) write(*,*) "Spectra (Sk) cut in: ", Sk

        RDF%randField(:,:) = sqrt((1.0d0)/((2.0d0*pi)**(RDF%nDim)))&
                             * RDF%randField(:,:)

        if(allocated(dgemm_mult))   deallocate(dgemm_mult)
        if(allocated(phiN))         deallocate(phiN);
        if(allocated(psiN))         deallocate(psiN);
        if(allocated(thetaN))       deallocate(thetaN);
        if(allocated(gammaN))       deallocate(gammaN);
        if(allocated(rVec))         deallocate(rVec);

    end subroutine gen_Std_Gauss_Isotropic

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss_FFT(RDF)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF

        !LOCAL
        double precision, dimension(:,:,:), allocatable :: gammak, phik
        double complex, dimension(:,:,:), allocatable :: Dk, Dk_base
        double precision, dimension(:,:,:), allocatable :: realOut
        integer, dimension(3)  :: DkStart
        integer, dimension(3) :: posVec, kS, kE, kSc, kEc;
        integer :: pos, i
        integer*8 plan, planTest2
        double complex,   dimension(2000) :: testV1, testV2

!        !START TEST
!        testV1(:) = 0.0D0
!        testV2(:) = 0.0D0
!        !write(*,*) "size(testV1(686:1314)) = ", size(testV1(686:1314))
!        !write(*,*) "size([(sin(dble(i)/1000), i =0, 628)]) = ", size([(sin(dble(i)/1000), i =0, 628)])
!        testV1(686:1314) = [(sin(dble(i)/100), i =0, 628)]
!
!
!        open (unit = 98 , file = "beforeFFT", action = 'write')
!        open (unit = 99 , file = "afterFFT", action = 'write')
!        open (unit = 100 , file = "afteriFFT", action = 'write')
!
!        write(*,*) "BEFORE FFT testV1"
!        write(98,"(2F15.5)") testV1
!
!
!        call dfftw_plan_dft_1d(planTest2, 2000, &
!                               testV1, testV2, FFTW_FORWARD, FFTW_ESTIMATE)
!        call dfftw_execute(planTest2)
!
!        write(*,*) "AFTER FFT testV2"
!        write(99,"(2F15.5)") testV2
!
!        call dfftw_plan_dft_1d(planTest2, 2000, &
!                               testV2, testV1, FFTW_BACKWARD, FFTW_ESTIMATE)
!        call dfftw_execute(planTest2)
!        write(*,*) "AFTER iFFT testV1"
!        write(100,"(2F15.5)") testV1
!
!        call dfftw_destroy_plan(planTest2)!
!
!        close (98)
!        close (99)
!        close (100)
!        !END TEST

        !write(*,*) "AFTER FFT testV1 = ", testV1

        write(*,*) "Inside FFT"
        write(get_fileId(),*) "Inside FFT"

!        call init_random_seed(RDF%seed)
!

!        call dispCarvalhol(RDF%kMax, "RDF%kMax")
!        call dispCarvalhol(RDF%kNStep, "RDF%kNStep")
!        call dispCarvalhol(RDF%kDelta, "RDF%kDelta")

        write(get_fileId(),*) "Defining kPoints and SkVec"
        call set_kPoints(RDF)
        call set_SkVec(RDF)
!
        allocate(gammak (RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
        allocate(phik   (RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
        allocate(Dk_base(RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
        allocate(Dk     (2*RDF%kNStep(1)-2, 2*RDF%kNStep(2)-2, 2*RDF%kNStep(3)-2))
!        allocate(realOut(RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
!
        call random_number(gammak(:,:,:))
        call random_number(phik(:,:,:))

        Dk_base = gammak*sqrt(RDF%Sk3D)*exp(2*pi*phik);

        write(get_fileId(),*) "Symmetrization"
        !Symmetrization
        kS  = 2
        kE  = RDF%kNStep-1
        kSc = RDF%kNStep + 1;
        kEc = 2*RDF%kNStep -2;

        write(get_fileId(),*) " RDF%kNStep = ", RDF%kNStep
        write(get_fileId(),*) " shape(Dk_base) = ", shape(Dk_base)
        write(get_fileId(),*) " shape(Dk) = ", shape(Dk)
        write(get_fileId(),*) " kS  = ", kS
        write(get_fileId(),*) " kE  = ", kE
        write(get_fileId(),*) " kSc = ", kSc
        write(get_fileId(),*) " kEc = ", kEc

        write(get_fileId(),*) "Copy"
        !Copy
        Dk(1:RDF%kNStep(1), 1:RDF%kNStep(2), 1:RDF%kNStep(3)) = Dk_base(:,:,:)

        Dk(kSc(1):kEc(1)  , 1:RDF%kNStep(2), 1:RDF%kNStep(3)) = Dk_base(kS(1):kE(1),:,:)
        Dk(1:RDF%kNStep(1), kSc(2):kEc(2)  , 1:RDF%kNStep(3)) = Dk_base(:,kS(2):kE(2),:)
        Dk(1:RDF%kNStep(1), 1:RDF%kNStep(2), kSc(3):kEc(3)  ) = Dk_base(:,:,kS(3):kE(3))

        Dk(kSc(1):kEc(1)  , kSc(2):kEc(2)  , 1:RDF%kNStep(3)) = Dk_base(kS(1):kE(1), kS(2):kE(2), :          ) ![-1,-1, 1]
        Dk(1:RDF%kNStep(1), kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = Dk_base(:          , kS(2):kE(2), kS(3):kE(3)) ![ 1,-1,-1]
        Dk(kSc(1):kEc(1)  , 1:RDF%kNStep(2), kSc(3):kEc(3)  ) = Dk_base(kS(1):kE(1), :          , kS(3):kE(3)) ![-1, 1,-1]

        Dk(kSc(1):kEc(1)  , kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = Dk_base(kS(1):kE(1), kS(2):kE(2), kS(3):kE(3)) ![-1,-1,-1]

        write(get_fileId(),*) "Hermitian Conjugate"
        !Hermitian Conjugate
        Dk(kSc(1):kEc(1), :, :) = -conjg(Dk(kEc(1):kSc(1):-1, :, :))
        Dk(:, kSc(2):kEc(2), :) = -conjg(Dk(:, kEc(2):kSc(2):-1, :))
        Dk(:, :, kSc(3):kEc(3)) = -conjg(Dk(:, :, kEc(3):kSc(3):-1))

        write(get_fileId(),*) "FFT Execution"
        !FFT Execution
        !call dfftw_plan_dft_3d(plan, size(Dk,1), size(Dk,2), size(Dk,3), &
        !                       Dk, Dk, FFTW_BACKWARD, FFTW_ESTIMATE)
        !
        !call dfftw_execute(plan)
        !call dfftw_destroy_plan(plan)

        write(get_fileId(),*) "Mapping"
        do pos = 1, RDF%xNTotal
             write(get_fileId(),*) "pos = ", pos
            call find_Permutation(pos, RDF%kNStep, posVec)
            !write(*,*) "pos = ", pos, "posVec = ", posVec
            write(get_fileId(),*) "allocated(RDF%randField) = ", allocated(RDF%randField)
            RDF%randField(pos,1) = real(Dk(posVec(1), posVec(2), posVec(3)))
        end do

!        Dk(1:RDF%kNStep(1), 1:RDF%kNStep(2), 1:RDF%kNStep(3)) =       Dk_base(:             , :             , :             ) ![ 1, 1, 1]

!        Dk(kSc(1):kEc(1)  , 1:RDF%kNStep(2), 1:RDF%kNStep(3)) = conjg(Dk_base(kE(1):kS(1):-1, :             , :             )) ![-1, 1, 1]
!        Dk(1:RDF%kNStep(1), kSc(2):kEc(2)  , 1:RDF%kNStep(3)) = conjg(Dk_base(:             , kE(2):kS(2):-1, :             )) ![ 1,-1, 1]
!        Dk(1:RDF%kNStep(1), 1:RDF%kNStep(2), kSc(3):kEc(3)  ) = conjg(Dk_base(:             , :             , kE(3):kS(3):-1)) ![ 1, 1,-1]
!
!        Dk(kSc(1):kEc(1)  , kSc(2):kEc(2)  , 1:RDF%kNStep(3)) = conjg(Dk_base(kE(1):kS(1):-1, kE(2):kS(2):-1, :             )) ![-1,-1, 1]
!        Dk(1:RDF%kNStep(1), kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = conjg(Dk_base(:             , kE(2):kS(2):-1, kE(3):kS(3):-1)) ![ 1,-1,-1]
!        Dk(kSc(1):kEc(1)  , 1:RDF%kNStep(2), kSc(3):kEc(3)  ) = conjg(Dk_base(kE(1):kS(1):-1, :             , kE(3):kS(3):-1)) ![-1, 1,-1]
!
!        Dk(kSc(1):kEc(1)  , kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = conjg(Dk_base(kE(1):kS(1):-1, kE(2):kS(2):-1, kE(3):kS(3):-1)) ![-1,-1,-1]
!
!        testV1(1:8)  = Dk_base(1:8,1,1)
!        testV1(9:14) = conjg(Dk_base(7:2:-1,1,1))

!        write(*,*) "BEFORE FFT testV1 = ", testV1
!
!        call dfftw_plan_dft_1d(plan, 14, &
!                               testV1, testV1, FFTW_BACKWARD, FFTW_ESTIMATE)
!        call dfftw_execute(plan)
!        call dfftw_destroy_plan(plan)!
!
!        write(*,*) "AFTER FFT testV1 = ", testV1



!        write(*,*) " Corners 1st quadrant"
!        write(*,*) " Dk(1,1,1) = ", Dk(2, 2, 2)
!        write(*,*) " Dk() = ", Dk(kEc(1), 2, 2)
!        write(*,*) " Dk() = ", Dk(2, kEc(2), 2)
!        write(*,*) " Dk() = ", Dk(2, 2, kEc(3))
!        write(*,*) " Dk() = ", Dk(kEc(1), kEc(2), 2)
!        write(*,*) " Dk() = ", Dk(2, kEc(2), kEc(3))
!        write(*,*) " Dk() = ", Dk(kEc(1), 2, kEc(3))
!        write(*,*) " Dk(N,N,N) = ", Dk(kEc(1), kEc(2), kEc(3))

        !write(*,*) "Dk 111    = ", Dk(1:5, 1:5, 1:5)
        !write(*,*) "Dk -1-1-1 = ",Dk(kSc(1)-5:kSc(1)  , kSc(2)-5:kSc(2)  , kSc(3)-5:kSc(3)  )

!        write(*,*) "BEFORE Dk(1:10,1:10,1)", Dk(1:10,1:10,1)
!
!        call dfftw_plan_dft_3d(plan, RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3), &
!                               Dk_base, Dk_base, FFTW_BACKWARD, FFTW_ESTIMATE)
!        call fftw_plan_dft_c2r_3d(plan, kEc(1), kEc(2), kEc(3), &
!                                  Dk, realOut, FFTW_ESTIMATE)


!        call dfftw_plan_dft_3d(plan, kEc(1), kEc(2), kEc(3), &
!                               Dk, Dk, FFTW_BACKWARD, FFTW_ESTIMATE)
!        call dfftw_execute(plan)
!        call dfftw_destroy_plan(plan)
!
!        realOut(:,:,:) = real(Dk(1:RDF%kNStep(1), 1:RDF%kNStep(2), 1:RDF%kNStep(3)))
!
!        !write(*,*) "AFTER Dk(1:10,1:10,1)", Dk(1:10,1:10,1)
!        write(*,*) "AFTER realOut(1:10,1:10,1)", realOut(1:10,1:10,1)
!
!        !Mapping generated field back to coordinates
!        DkStart = [1, 1, 1]
!
!        do pos = 1, RDF%xNTotal
!            call find_Permutation(pos, RDF%kNStep, posVec)
!            !write(*,*) "pos = ", pos, "posVec = ", posVec
!            RDF%randField(pos,1) = realOut(posVec(1), posVec(2), posVec(3))
!        end do
!
        deallocate(gammak)
        deallocate(phik)
        deallocate(Dk)
        deallocate(Dk_base)
        !deallocate(realOut)

    end subroutine gen_Std_Gauss_FFT


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_neighbours_info (RDF, MSH)
        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF

        !INPUT
        type(MESH), intent(in) :: MSH

        !LOCAL
        integer :: neighPos, stage, direction
        integer, allocatable, dimension(:) :: tempSeed
        integer :: request, code, tag
        integer, dimension(MPI_STATUS_SIZE) :: status

        do neighPos = 1, size(MSH%neigh)
            if(MSH%neigh(neighPos) < 0) cycle !Check if this neighbour exists
            call calculate_random_seed(tempSeed, RDF%seedStart+MSH%neigh(neighPos))
            RDF%neighSeed(:,neighPos) = tempSeed
        end do

        do stage = 1, 2 !Sending and then receiving

            do neighPos = 1, size(MSH%neigh)

                if(MSH%neigh(neighPos) < 0) cycle !Check if this neighbour exists

                if(stage == 1) then
                    tag = findTag(MSH, neighPos, neighPos, send = .true.)
                    call MPI_ISEND (RDF%xMaxBound(:)-RDF%xMinBound(:), RDF%nDim, MPI_DOUBLE_PRECISION, &
                        MSH%neigh(neighPos), tag, RDF%comm, request, code)
                else if(stage == 2) then
                    tag = findTag(MSH, neighPos, neighPos, send = .false.)
                    call MPI_RECV (RDF%neighRange(:,neighPos), RDF%nDim, MPI_DOUBLE_PRECISION, &
                        MSH%neigh(neighPos), tag, RDF%comm, status, code)
                end if

            end do

        end do

        deallocate(tempSeed)

    end subroutine get_neighbours_info

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine getNeighIndexRange(MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)
        !INPUT
        type(MESH), intent(in) :: MSH
        !OUTPUT
        integer, intent(out) :: minIndexNeigh, maxIndexNeigh
        logical, dimension(:), intent(out) :: considerNeighbour

        considerNeighbour = .true.
        where(MSH%neigh < 0) considerNeighbour = .false.
        where(minval(MSH%neighShift,1) < 0) considerNeighbour = .false.

        !Global Min and Max positions
        minIndexNeigh = minval(pack(MSH%indexNeigh(1,:), considerNeighbour))
        maxIndexNeigh = maxval(pack(MSH%indexNeigh(2,:), considerNeighbour))

    end subroutine getNeighIndexRange

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine takeNeighboursContribution(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour, normFactor)
        !INPUT
        type(MESH), intent(in) :: MSH
        integer, intent(in) :: minIndexNeigh, maxIndexNeigh
        logical, dimension(:), intent(in) ::considerNeighbour
        double precision, dimension(minIndexNeigh:maxIndexNeigh), intent(in) :: normFactor
        !OUTPUT
        type(RF), intent(inout) :: RDF
        !LOCAL
        double precision, dimension(MSH%nDim) :: neighOrCorner
        logical :: sndRcv
        integer :: i, direction, neighPos, minPos, maxPos
        double precision, dimension(minIndexNeigh:maxIndexNeigh) :: power
        type(RF) :: tmpRDF


        !Taking the contributions from neighbours------------------------------------------------------
        write(get_fileId(), *) "   Taking the contributions from neighbours "
        call init_RF(tmpRDF, RDF%nDim, RDF%Nmc, RDF%comm, RDF%rang, RDF%nb_procs)
        call copy_RF_properties(RDF, tmpRDF)

        do direction = 1, size(MSH%neigh)

            if(.not. considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction

            write(get_fileId(), *) " "
            write(get_fileId(), *) "  DIRECTION      = ", direction
            write(get_fileId(), *) "  Neighbour Rank = ", MSH%neigh(direction)

            !write(get_fileId(),*) "Direction = ", direction

            !Preparing the xPoints of a given direction
            call copy_RF_xPoints(MSH, RDF, tmpRDF, tmpRDF%xPoints_Local, direction)
            call allocate_randField(tmpRDF, tmpRDF%randField_Local)
            minPos = MSH%indexNeigh(1,direction)
            maxPos = MSH%indexNeigh(2,direction)

            do neighPos = 1, size(MSH%neigh)

                if(.not. considerNeighbour(neighPos)) cycle !Don't consider Neighbours in this direction

                sndRcv = .true.

                !Checking if this neighbour should be taking into account in this part of the field
                do i = 1, MSH%nDim
                    if(       (MSH%neighShift(i, neighPos) /= 0) &
                        .and. (MSH%neighShift(i, neighPos) /= MSH%neighShift(i, direction))) then
                        sndRcv = .false.
                        exit
                    end if
                end do

                !From this point we know that we want the contribution of this neighbour in this direction
                if (sndRcv) then

                    tmpRDF%xMinBound = 0.0D0
                    tmpRDF%xMaxBound = RDF%neighRange(:,neighPos)
                    tmpRDF%seed      = RDF%neighSeed(:,neighPos)

                    !Generating Standard Gaussian Field
                    select case (tmpRDF%method)
                        case(SHINOZUKA)
                            call gen_Std_Gauss_Shinozuka(tmpRDF)
                        case(ISOTROPIC)
                            call gen_Std_Gauss_Isotropic(tmpRDF)
                        case(RANDOMIZATION)
                            call gen_Std_Gauss_Randomization(tmpRDF)
                    end select

                    !Finding origin for Shape Function
                    neighOrCorner = originCorner + MSH%overlap*MSH%neighShift(:, neighPos)

                    !Shape Function multiplication and sum of the contribution
                    power(minPos:maxPos) = 0
                    do i = 1, MSH%nDim
                        if(MSH%neighShift(i, direction) == 0) cycle

                        power(minPos:maxPos) = ((tmpRDF%xPoints(i, :) - neighOrCorner(i))/MSH%overlap)**2 &
                                                   + power(minPos:maxPos)
                    end do

                    power(minPos:maxPos) = sqrt(power(minPos:maxPos))
                    where (power > 1.0D0) power = 1.0

                    RDF%randField(minPos:maxPos,1) = RDF%randField(minPos:maxPos,1) + &
                                                     (tmpRDF%randField(:,1) * sqrt((cos((PI)*(power(minPos:maxPos)))+1.0D0)/2.0D0))
                end if
            end do !Neighbours
        end do !Directions

        !Normalizing Field---------------------------------------------------------------
        RDF%randField(minIndexNeigh:maxIndexNeigh,1) = RDF%randField(minIndexNeigh:maxIndexNeigh,1)/normFactor(:)

        call finalize_RF(tmpRDF)

    end subroutine takeNeighboursContribution

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine applyWeightingFunctions(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)
        !INPUT
        type(RF), intent(in) :: RDF
        type(MESH), intent(in) :: MSH
        integer, intent(in) :: minIndexNeigh, maxIndexNeigh
        logical, dimension(:), intent(in) ::considerNeighbour

        !LOCAL
        integer :: direction, minPos, maxPos, i
        double precision, dimension(minIndexNeigh:maxIndexNeigh) :: powerWF
        double precision, dimension(MSH%nDim) :: originCorner

        !Modify extremes of local Random Field-------------------------------------------------------

        !Loop on directions
        do direction = 1, size(MSH%neigh)

            if(.not. considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction

            !Positions in temp Vector
            minPos = MSH%indexNeigh(1,direction)
            maxPos = MSH%indexNeigh(2,direction)

            !Finding origin
            originCorner = MSH%xMinNeigh(:, direction)
            where(MSH%neighShift(:, direction) == -1) originCorner = MSH%xMaxNeigh(:, direction)
            !write(get_fileId(), *) "  originCorner = ", originCorner

            !Shape Function multiplication (Obs: redefinition of power, should come after Normalization values calculation)
            powerWF(minPos:maxPos) = 0.0D0
            do i = 1, MSH%nDim
                if(MSH%neighShift(i, direction) == 0) cycle

                powerWF(minPos:maxPos) = ((RDF%xPoints(i, minPos:maxPos) - originCorner(i))/MSH%overlap(i))**2 &
                                         + powerWF(minPos:maxPos)
            end do

            powerWF(minPos:maxPos) = sqrt(powerWF(minPos:maxPos))
            where (powerWF(minPos:maxPos) > 1.0D0) powerWF(minPos:maxPos) = 1.0D0

        end do !Direction

        RDF%randField(minIndexNeigh:maxIndexNeigh,1) = RDF%randField(minIndexNeigh:maxIndexNeigh,1) &
                                                 * sqrt((cos(PI*(powerWF))+1.0D0)/2.0D0)

    end subroutine applyWeightingFunctions

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine findNormalizationFactors(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour, normFactor)

        implicit none

        !INPUT
        type(RF), intent(in) :: RDF
        type(MESH), intent(in) :: MSH
        integer, intent(in) :: minIndexNeigh, maxIndexNeigh
        logical, dimension(:), intent(in) ::considerNeighbour

        !OUTPUT
        double precision, dimension(minIndexNeigh:maxIndexNeigh), intent(out) :: normFactor

        !LOCAL
        integer :: i, j, n
        double precision, dimension(MSH%nDim, 2**MSH%nDim) :: originList
        double precision, dimension(minIndexNeigh:maxIndexNeigh) :: power
        integer, dimension(MSH%nDim) :: tempXNStep
        double precision :: distance
        integer :: minPos, maxPos
        integer :: direction, zeroVal

        !Find normalization values-------------------------------------------------------
        normFactor(:) = 0
        power(:) = 0

        !Loop on directions
        write(get_fileId(), *) "   Finding normalization Factors "
        do direction = 1, size(MSH%neigh)

            if(.not. considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction

            !Positions in temp Vector
            minPos = MSH%indexNeigh(1,direction)
            maxPos = MSH%indexNeigh(2,direction)



            !Finding all vertices
            tempXNStep = 2
            do i = 1, 2**MSH%nDim
                call get_Permutation(i, MSH%xMaxNeigh(:,direction), tempXNStep, originList(:, i), &
                                     MSH%xMinNeigh(:,direction), snapExtremes = .true.);
            end do

            !Normalization values
            normFactor(minPos:maxPos) = 0 !TO DELETE
            zeroVal = MSH%nDim - (sum((MSH%neighShift(:, direction))**2))

            !We sweep all vertices in this direction
            do j = 1, size(originList, 2)
                power(minPos:maxPos) = 0
                do i = 1, MSH%nDim
                    if(MSH%neighShift(i, direction) == 0) cycle

                    power(minPos:maxPos) = ((RDF%xPoints(i, minPos:maxPos) - originList(i, j))/MSH%overlap(i))**2 &
                                               + power(minPos:maxPos);
                end do

                power(minPos:maxPos) = sqrt(power(minPos:maxPos))
                where (power(minPos:maxPos) > 1.0D0) power(minPos:maxPos) = 1.0D0
                normFactor(minPos:maxPos) = sqrt((cos((PI)*(power(minPos:maxPos)))+1.0D0)/2.0D0) + normFactor(minPos:maxPos)


            end do !Neighbours

            normFactor(minPos:maxPos) = normFactor(minPos:maxPos)/(2.0**zeroVal)

        end do !Direction

    end subroutine findNormalizationFactors

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine create_RF_overlap(RDF, MSH)
        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF

        !INPUT
        type(MESH), intent(in) :: MSH

!        !LOCAL
!        integer :: i, j, n
!        double precision, dimension(:, :), allocatable :: originList
!        double precision, dimension(:), allocatable :: normFactor, power, powerSF
!        double precision, dimension(:), allocatable :: originCorner, neighOrCorner
!        integer :: minPos, maxPos, neighPos, minPosGlob, maxPosGlob
!        integer, dimension(:), allocatable :: tempXNStep, shift, counterDir
!        integer :: shiftScal
!        integer :: direction, zeroVal, nVertices
!        logical, dimension(:,:), allocatable :: first
!        logical ::sndRcv
!        logical, dimension(size(MSH%neigh)) :: considerNeighbour
!        type(RF) :: tmpRDF
!
!        allocate(originList(MSH%nDim, 2**MSH%nDim))
!        allocate(originCorner(MSH%nDim))
!        allocate(neighOrCorner(MSH%nDim))
!        allocate(tempXNStep(MSH%nDim))
!        allocate(shift(MSH%nDim))
!        allocate(first(-1:1, MSH%nDim))
!
!        considerNeighbour = .true.
!        do direction = 1, size(MSH%neigh)
!            if(minval(MSH%neighShift(:,direction)) == -1) considerNeighbour(direction) = .false.
!            if(MSH%neigh(direction) < 0) considerNeighbour(direction) = .false.
!        end do
!
!        minPosGlob = minval(pack(MSH%indexNeigh(1,:), considerNeighbour))
!        maxPosGlob = maxval(pack(MSH%indexNeigh(2,:), considerNeighbour))
!        allocate(normFactor(minPosGlob:maxPosGlob)) !Obs: first index doesn't start in "1"
!        allocate(power(minPosGlob:maxPosGlob)) !Obs: first index doesn't start in "1"
!        allocate(powerSF(minPosGlob:maxPosGlob))
!        allocate(counterDir(MSH%nDim))
!
!        !Modify extremes of local Random Field and find normalization values-------------------------------------------------------
!        normFactor(:) = 0
!        power(:) = 0
!
!        !Loop on directions
!        do direction = 1, size(MSH%neigh)
!
!            if(MSH%neigh(direction) <0) cycle !Don't consider Neighbours in this direction
!
!            !Positions in temp Vector
!            minPos = MSH%indexNeigh(1,direction)
!            maxPos = MSH%indexNeigh(2,direction)
!
!            !Finding origin
!            originCorner = MSH%xMinNeigh(:, direction)
!            do i = 1, MSH%nDim
!                if(MSH%neighShift(i, direction) == -1) then
!                    originCorner(i) = MSH%xMaxNeigh(i, direction)
!                end if
!            end do
!
!            !write(get_fileId(), *) "  originCorner = ", originCorner
!
!            !Finding all vertices
!            tempXNStep = 2
!            do i = 1, 2**MSH%nDim
!                call get_Permutation(i, MSH%xMaxNeigh(:,direction), tempXNStep, originList(:, i), &
!                                     MSH%xMinNeigh(:,direction), snapExtremes = .true.);
!            end do
!
!            !Normalization values
!            first = .true.
!            normFactor(minPos:maxPos) = 0
!            power(minPos:maxPos) = 0
!            zeroVal = MSH%nDim - (sum((MSH%neighShift(:, direction))**2))
!
!
!            !We sweep al vertices in this direction
!            do j = 1, size(originList, 2)
!
!                power(minPos:maxPos) = 0
!
!                do i = 1, MSH%nDim
!
!                    if(MSH%neighShift(i, direction) == 0) cycle
!
!                    power(minPos:maxPos) = (RDF%xPoints(i, minPos:maxPos) - originList(i, j))**2 &
!                                               + power(minPos:maxPos);
!
!                end do
!
!                where (sqrt(power(minPos:maxPos)) > MSH%overlap) power(minPos:maxPos) = MSH%overlap**2
!                power(minPos:maxPos) = sqrt(power(minPos:maxPos)/(MSH%overlap**2))
!                normFactor(minPos:maxPos) = sqrt((cos((PI)*(power(minPos:maxPos)))+1.0D0)/2.0D0) + normFactor(minPos:maxPos)
!
!
!            end do !Neighbours
!
!            normFactor(minPos:maxPos) = normFactor(minPos:maxPos)/(2.0**zeroVal)
!
!            !Shape Function multiplication (Obs: redefinition of power, should come after Normalization values calculation)
!            powerSF(minPos:maxPos) = 0
!            do i = 1, MSH%nDim
!                if(MSH%neighShift(i, direction) == 0) cycle
!
!                powerSF(minPos:maxPos) = (RDF%xPoints(i, minPos:maxPos) - originCorner(i))**2 &
!                                       + powerSF(minPos:maxPos)
!            end do
!
!            where (sqrt(powerSF(minPos:maxPos)) > MSH%overlap) powerSF(minPos:maxPos) = MSH%overlap**2
!            powerSF(minPos:maxPos) = sqrt(powerSF(minPos:maxPos)/(MSH%overlap**2))
!
!        end do !Direction

!        RDF%randField(minPosGlob:maxPosGlob,1) = RDF%randField(minPosGlob:maxPosGlob,1) &
!                                                 * sqrt((cos(PI*(powerSF))+1.0D0)/2.0D0)

!        RDF%randField(minPosGlob:maxPosGlob,1) = powerSF


!        !Taking the contributions from neighbours------------------------------------------------------
!        write(get_fileId(), *) "   Taking the contributions from neighbours "
!        call init_RF(tmpRDF, RDF%nDim, RDF%Nmc, RDF%comm, RDF%rang, RDF%nb_procs)
!        call copy_RF_properties(RDF, tmpRDF)
!
!        do direction = 1, size(MSH%neigh)
!
!            if(.not. considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction
!
!            write(get_fileId(), *) " "
!            write(get_fileId(), *) "  DIRECTION      = ", direction
!            write(get_fileId(), *) "  Neighbour Rank = ", MSH%neigh(direction)
!
!            !write(get_fileId(),*) "Direction = ", direction
!
!            !Finding origin
!            originCorner = MSH%xMinNeigh(:, direction)
!            do i = 1, MSH%nDim
!                if(MSH%neighShift(i, direction) == -1) then
!                    originCorner(i) = MSH%xMaxNeigh(i, direction)
!                end if
!            end do
!
!            !write(get_fileId(),*) "originCorner = ", originCorner
!
!            !Preparing the xPoints of a given direction
!            call copy_RF_xPoints(MSH, RDF, tmpRDF, tmpRDF%xPoints_Local, direction)
!            call allocate_randField(tmpRDF, tmpRDF%randField_Local)
!            minPos = MSH%indexNeigh(1,direction)
!            maxPos = MSH%indexNeigh(2,direction)
!
!            do neighPos = 1, size(MSH%neigh)
!
!                if(.not. considerNeighbour(neighPos)) cycle !Don't consider Neighbours in this direction
!
!                sndRcv = .true.
!
!                !Checking if this neighbour should be taking into account in this part of the field
!                do i = 1, MSH%nDim
!                    if(       (MSH%neighShift(i, neighPos) /= 0) &
!                        .and. (MSH%neighShift(i, neighPos) /= MSH%neighShift(i, direction))) then
!                        sndRcv = .false.
!                        exit
!                    end if
!                end do
!
!                !From this point we know that we want the contribution of this neighbour in this direction
!                if (sndRcv) then
!
!                    !write(get_fileId(),*) "sndRcv = ", sndRcv
!
!                    tmpRDF%xMinBound = 0.0D0
!                    tmpRDF%xMaxBound = RDF%neighRange(:,neighPos)
!                    tmpRDF%seed      = RDF%neighSeed(:,neighPos)
!
!                    !Generating Standard Gaussian Field
!                    !write(get_fileId(),*) "Generating Standard Gaussian Field"
!                    select case (tmpRDF%method)
!                        case(SHINOZUKA)
!                            write(get_fileId(),*) "SHINOZUKA"
!                            call gen_Std_Gauss_Shinozuka(tmpRDF)
!                        case(ISOTROPIC)
!                            call gen_Std_Gauss_Isotropic(tmpRDF)
!                        case(RANDOMIZATION)
!                            call gen_Std_Gauss_Randomization(tmpRDF)
!                    end select
!
!                    !Finding origin for Shape Function
!                    neighOrCorner = originCorner + MSH%overlap*MSH%neighShift(:, neighPos)
!
!                    !Shape Function multiplication and sum of the contribution
!                    power(minPos:maxPos) = 0
!                    do i = 1, MSH%nDim
!                        if(MSH%neighShift(i, direction) == 0) cycle
!
!                        power(minPos:maxPos) = (tmpRDF%xPoints(i, :) - neighOrCorner(i))**2 &
!                                                   + power(minPos:maxPos)
!                    end do
!
!                    where (sqrt(power(minPos:maxPos)) > MSH%overlap) power(minPos:maxPos) = MSH%overlap**2
!                    power(minPos:maxPos) = power(minPos:maxPos)/(MSH%overlap**2)
!
!                    RDF%randField(minPos:maxPos,1) = RDF%randField(minPos:maxPos,1) + &
!                                                     (tmpRDF%randField(:,1) * sqrt((cos((PI)*(power(minPos:maxPos)))+1.0D0)/2.0D0))
!                end if
!            end do !Neighbours
!        end do !Directions
!
!        !Normalizing Field---------------------------------------------------------------
!
!        RDF%randField(minPosGlob:maxPosGlob,1) = RDF%randField(minPosGlob:maxPosGlob,1)/normFactor(:)

    end subroutine create_RF_overlap
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine reorderRandomFieldStruct(RDF, MSH)
        implicit none
        
        type(RF)  , intent(inout) :: RDF
        type(MESH), intent(in)    :: MSH

        !LOCAL
        integer, dimension(MSH%nDim) :: ind3D, offset
        double precision, dimension(MSH%nDim) :: orig, coordValue, tempCoordValue, xStep
        double precision, dimension(RDF%Nmc)  :: tempRFValue, RFValue
        integer :: ind1D
        integer :: i

        orig  = (dble(RDF%origin-1)*MSH%xStep + MSH%xMinGlob)
        xStep = MSH%xStep

        offset(1) = 1
        do i = 2, MSH%nDim
          offset(i) = product(MSH%xNStep(1:i-1))
        end do

        !write(get_fileId(),*) "  "
        !write(get_fileId(),*) " orig        = ", orig
        !write(get_fileId(),*) " offset      = ", offset
        !write(get_fileId(),*) " RDF%origin  = ", RDF%origin
        !write(get_fileId(),*) " xStep       = ", xStep
        !write(get_fileId(),*) " MSH%xNTotal = ", MSH%xNTotal
        !write(get_fileId(),*) " MSH%xNStep  = ", MSH%xNStep
        !write(get_fileId(),*) "  "

        do i = 1, size(RDF%randField,1)
          coordValue = RDF%xPoints(:,i)
          RFValue    = RDF%randField(i,:)
          ind3D = nint((coordValue-orig)/xStep)
          ind1D = sum(ind3D*offset) + 1

          !The point is not where it was supposed to be
          do while (ind1D /= i)
            !Finding index
            ind3D = nint((coordValue-orig)/xStep)
            ind1D = sum(ind3D*offset) + 1
            !Saving temp data
            tempRFValue    = RDF%randField(ind1D,:)
            tempCoordValue = RDF%xPoints(:, ind1D)
            !Replacing data
            RDF%randField(ind1D,:) = RFvalue
            RDF%xPoints(:, ind1D)  = coordValue
            !Going to the next coordinate (the one that was in the index we took)
            RFValue     = tempRFValue
            coordValue  = tempCoordValue
          end do

          
       end do

    end subroutine reorderRandomFieldStruct
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine copy_RF_properties(orRDF, destRDF)
        implicit none

        !INPUT AND OUTPUT
        type(RF), intent(in)   :: orRDF
        type(RF) ::destRDF

        destRDF%nDim     = orRDF%nDim
        destRDF%Nmc      = orRDF%Nmc
        destRDF%comm     = orRDF%comm
        destRDF%rang     = orRDF%rang
        destRDF%nb_procs = orRDF%nb_procs
        destRDF%corrL    = orRDF%corrL
        destRDF%corrMod  = orRDF%corrMod
        destRDF%kMax     = orRDF%kMax
        destRDF%xMinGlob = orRDF%xMinGlob
        destRDF%xMaxGlob = orRDF%xMaxGlob
        destRDF%calculate   = orRDF%calculate
        destRDF%method      = orRDF%method
        destRDF%corrL       = orRDF%corrL
        destRDF%kMax        = orRDF%kMax
        destRDF%independent =  orRDF%independent

    end subroutine copy_RF_properties

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine copy_RF_xPoints(MSH, orRDF, destRDF, destXPoints, direction)
        implicit none

        !INPUT AND OUTPUT
        type(MESH), intent(in) :: MSH
        type(RF), intent(in)   :: orRDF
        type(RF) ::destRDF
        integer, intent(in) :: direction
        double precision, dimension(:, :), allocatable, intent(out), target :: destXPoints;

        !LOCAL
        integer :: i, j, totalSize

        totalSize = MSH%indexNeigh(2,direction) - MSH%indexNeigh(1,direction) + 1

        if(allocated(destXPoints)) then
            if(size(destXPoints) /= totalSize) then
                nullify(destRDF%xPoints)
                deallocate(destXPoints)
            end if
        end if

        if(.not.allocated(destXPoints)) allocate(destXPoints(MSH%nDim, 1:totalSize))

        destRDF%xPoints => destXPoints

        destRDF%xPoints(:,:) = orRDF%xPoints(:, MSH%indexNeigh(1,direction):MSH%indexNeigh(2,direction))
        destRDF%xNTotal = totalSize

    end subroutine copy_RF_xPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function findTag(MSH, neighPos, direction, send) result(tag)

        implicit none
        !INPUT
        type(MESH), intent(in) :: MSH
        integer, intent(in) :: neighPos, direction
        logical, intent(in) :: send
        integer :: sinal

        !OUTPUT
        integer :: tag

        !LOCAL
        double precision :: i


        if(MSH%neigh(neighPos) < 0 .or. MSH%neigh(direction) < 0) then
            write(*,*) "Inside findTag , Invalid Neighbour"
            tag = -1
        else
            tag = 0

            if(send) then
                do i = 1, MSH%nDim
                    tag = tag + MSH%neighShift(i, direction)*3**(MSH%nDim - i)
                end do
            else
                do i = 1, MSH%nDim
                    sinal = 1
                    if(MSH%neighShift(i, neighPos) /= 0) sinal = -1
                    tag = tag + sinal*MSH%neighShift(i, direction)*3**(MSH%nDim - i)
                end do
            end if

            tag = tag + (3**(MSH%nDim)-1)/2

        end if

    end function findTag

end module randomFieldND
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
