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



contains
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
        double complex, dimension(:,:,:), allocatable :: Dk_base
        !double complex, dimension(:,:,:), allocatable :: Dk
        double precision, dimension(:,:,:), allocatable :: realOut
        integer, dimension(3)  :: DkStart
        integer, dimension(3) :: posVec, kS, kE, kSc, kEc;
        integer :: pos, i
        integer*8 plan, planTest2
        double complex,   dimension(2000) :: testV1, testV2
        double precision :: sizesProd



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

        write(get_fileId(),*) "Inside FFT"

!        call init_random_seed(RDF%seed)
!

!        call dispCarvalhol(RDF%kMax, "RDF%kMax")
!        call dispCarvalhol(RDF%kNStep, "RDF%kNStep")
!        call dispCarvalhol(RDF%kDelta, "RDF%kDelta")

        write(get_fileId(),*) "Defining kPoints and SkVec"
        call set_kPoints(RDF)
        call set_SkVec(RDF)

        allocate(gammak (RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
        allocate(phik   (RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
        allocate(Dk_base(RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
        !allocate(Dk     (0:2*RDF%kNStep(1)-2-1, 0:2*RDF%kNStep(2)-2-1, 0:2*RDF%kNStep(3)-2-1))
        allocate(realOut(RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))

        call random_number(gammak(:,:,:))
        call random_number(phik(:,:,:))

        !Dk_base = product(xMaxGlob - xMinGlob) * gammak*sqrt(RDF%Sk3D)*exp(2*pi*phik);
        if(RDF%independent) then
            sizesProd = product(RDF%xMaxBound - RDF%xMinBound)
        else
            sizesProd = product(RDF%xMaxGlob - RDF%xMinGlob)
        end if

        Dk_base(:,:,:) = sqrt(sizesProd*RDF%Sk3D);


!        i = 0
!        write(get_fileId(),*) " Dk_base 3x ", i
!        write(get_fileId(),*) Dk_base (i,i,i)
        i = 1
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)
        i = 2
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)
        i = 3
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)
        i = 4
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)
        i = 5
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)


!        write(get_fileId(),*) "Symmetrization"
!        !Symmetrization
!        kS  = 2
!        kE  = RDF%kNStep-1
!        kSc = RDF%kNStep + 1;
!        kEc = 2*RDF%kNStep -2;
!
!        write(get_fileId(),*) " RDF%kNStep = ", RDF%kNStep
!        write(get_fileId(),*) " shape(Dk_base) = ", shape(Dk_base)
!        write(get_fileId(),*) " shape(Dk) = ", shape(Dk)
!        write(get_fileId(),*) " kS  = ", kS
!        write(get_fileId(),*) " kE  = ", kE
!        write(get_fileId(),*) " kSc = ", kSc
!        write(get_fileId(),*) " kEc = ", kEc
!
!        write(get_fileId(),*) "Copy"
!        !Copy
!        Dk(1:RDF%kNStep(1), 1:RDF%kNStep(2), 1:RDF%kNStep(3)) = Dk_base(:,:,:)
!
!        Dk(kSc(1):kEc(1)  , 1:RDF%kNStep(2), 1:RDF%kNStep(3)) = Dk_base(kS(1):kE(1),:,:)
!        Dk(1:RDF%kNStep(1), kSc(2):kEc(2)  , 1:RDF%kNStep(3)) = Dk_base(:,kS(2):kE(2),:)
!        Dk(1:RDF%kNStep(1), 1:RDF%kNStep(2), kSc(3):kEc(3)  ) = Dk_base(:,:,kS(3):kE(3))
!
!        Dk(kSc(1):kEc(1)  , kSc(2):kEc(2)  , 1:RDF%kNStep(3)) = Dk_base(kS(1):kE(1), kS(2):kE(2), :          ) ![-1,-1, 1]
!        Dk(1:RDF%kNStep(1), kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = Dk_base(:          , kS(2):kE(2), kS(3):kE(3)) ![ 1,-1,-1]
!        Dk(kSc(1):kEc(1)  , 1:RDF%kNStep(2), kSc(3):kEc(3)  ) = Dk_base(kS(1):kE(1), :          , kS(3):kE(3)) ![-1, 1,-1]
!
!        Dk(kSc(1):kEc(1)  , kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = Dk_base(kS(1):kE(1), kS(2):kE(2), kS(3):kE(3)) ![-1,-1,-1]
!
!        write(get_fileId(),*) "Hermitian Conjugate"
!        !Hermitian Conjugate
!        Dk(kSc(1):kEc(1), :, :) = -conjg(Dk(kEc(1):kSc(1):-1, :, :))
!        Dk(:, kSc(2):kEc(2), :) = -conjg(Dk(:, kEc(2):kSc(2):-1, :))
!        Dk(:, :, kSc(3):kEc(3)) = -conjg(Dk(:, :, kEc(3):kSc(3):-1))

        write(get_fileId(),*) " "
        write(get_fileId(),*) "iFFT Execution"
        !FFT Execution
        call dfftw_plan_dft_c2r(plan, RDF%nDim, RDF%kNStep, &
                               Dk_base, realOut, FFTW_ESTIMATE);
        call dfftw_execute(plan)
        call dfftw_destroy_plan(plan)

        !realOut = realOut*sqrt((2*PI)**RDF%nDim/sizesProd)/RDF%kNTotal
        realOut = realOut/RDF%kNTotal

        !write(*,*) " realOut = ", realOut

        write(get_fileId(),*) "AFTER iFFT"
        i = 1
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)
        i = 2
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)
        i = 3
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)
        i = 4
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)
        i = 5
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)

        write(get_fileId(),*) " "
        write(get_fileId(),*) "FFT Execution"
        call dfftw_plan_dft_r2c(plan, RDF%nDim, RDF%kNStep, &
                                 realOut, Dk_base, FFTW_ESTIMATE);
        Dk_base = Dk_base*RDF%kNTotal
        call dfftw_execute(plan)
        call dfftw_destroy_plan(plan)

        write(get_fileId(),*) "AFTER FFT"

        i = 1
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)
        i = 2
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)
        i = 3
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)
        i = 4
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)
        i = 5
        write(get_fileId(),*) " Dk_base 3x ", i
        write(get_fileId(),*) Dk_base (i,i,i)

        realOut = realOut*sqrt((2*PI)**RDF%nDim/sizesProd)

        write(get_fileId(),*) " "
        write(get_fileId(),*) "AFTER real mutiplication"
        i = 1
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)
        i = 2
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)
        i = 3
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)
        i = 4
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)
        i = 5
        write(get_fileId(),*) " realOut 3x ", i
        write(get_fileId(),*) realOut (i,i,i)


        write(get_fileId(),*) "Mapping"
        do pos = 1, RDF%xNTotal
            !write(get_fileId(),*) "pos = ", pos
            call find_Permutation(pos, RDF%kNStep, posVec)
            !posVec = posVec -1
            !write(*,*) "pos = ", pos, "posVec = ", posVec
            !write(get_fileId(),*) "allocated(RDF%randField) = ", allocated(RDF%randField)
            RDF%randField(pos,1) = (realOut(posVec(1), posVec(2), posVec(3)))
        end do

        deallocate(gammak)
        deallocate(phik)
        !deallocate(Dk)
        deallocate(Dk_base)
        deallocate(realOut)

    end subroutine gen_Std_Gauss_FFT

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
