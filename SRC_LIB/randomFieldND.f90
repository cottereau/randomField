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
    use fftw3
    !use blas
    !use, intrinsic :: iso_c_binding
    implicit none




contains
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

        implicit none

        !INPUT
        integer          , intent(in) :: margiFirst;
        double precision , intent(in) :: fieldAvg, fieldVar;

        !OUTPUT (IN)
        double precision, dimension(1:, 1:), intent(inout) :: randField;

        !LOCAL VARIABLES
        double precision :: normalVar, normalAvg
        integer          :: error, code, i

        select case (margiFirst)
        case(fom_GAUSSIAN)
            normalVar = fieldVar
            normalAvg = fieldAvg
        case(fom_LOGNORMAL)
            if(fieldAvg <= 0.0D0) then
                write(*,*) ""
                write(*,*) "ERROR - when using lognormal fieldAvg should be a positive number"
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            normalVar = log(1 + fieldVar/(fieldAvg**2))
            normalAvg = log(fieldAvg) - normalVar/2
        end select

        randField(:,:) = randField(:,:) * sqrt(normalVar) &
            + normalAvg;

        if (margiFirst == fom_LOGNORMAL) then
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

        !!write(get_fileId(),*) "Inside Shinozuka"

        randomK = .false.
        if(present(randomK_in)) randomK = randomK_in
        call init_random_seed(RDF%seed)

        !!write(get_fileId(),*) "Defining kPoints and SkVec"
        call set_kPoints(RDF)
        call set_SkVec(RDF)

        !!write(get_fileId(),*) "Calculating Fields"
        allocate(k_x_phi (RDF%xNTotal, 1))
        allocate(kSign (2**(RDF%nDim-1), RDF%nDim));
        allocate(kVec(RDF%nDim, 1))

        call set_kSign(kSign) !Set the sign permutations for kVec

        allocate(phiK (RDF%kNTotal, size(kSign,1)));

        if(randomK) then
            !write(get_fileId(),*) "-----Shinozuka, k random-----"
            ampMult = 2.0d0*sqrt(1/(RDF%kNTotal*(2.0d0*PI)**(dble(RDF%nDim))))
        else
            !write(get_fileId(),*) "-----Shinozuka, k discrete-----"
            ampMult = 2.0d0*sqrt(product(RDF%kDelta)/((2.0d0*PI)**(dble(RDF%nDim))))
        end if

        !write(get_fileId(),*) "     kNStep  = ", RDF%kNStep
        !write(get_fileId(),*) "     kNTotal = ", RDF%kNTotal
        !write(get_fileId(),*) "     xNTotal = ", size(RDF%xPoints, 2)

        RDF%randField(:,:) = 0.0d0;

        do n = 1, RDF%Nmc
            !write(get_fileId(),*) "  --Generating Field Number ", n

            if(.not. RDF%calculate(n)) cycle

            call random_number(phiK(:,:))
            phiK(:,:) = 2.0D0*pi*phiK(:,:)

            !write(get_fileId(),*) "     First PhiK = ", phiK(1,1)
            !write(get_fileId(),*) "     Last PhiK  = ", phiK(size(phiK,1), size(phiK,2))

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

        !write(get_fileId(),*) "-----Randomization-----"

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

        !write(get_fileId(),*) "-----Inside Isotropic-----"

        !r Definition
        call set_kMaxND(RDF%corrMod, rMaxVec)
        call set_kPoints(RDF)
        rMax = rMaxVec(1)
        rDelta  = maxval(RDF%kDelta(:))/5.0D0 !Delta min in between two wave numbers to avoid periodicity
        rNTotal = ceiling(rMax/rDelta) + 1;

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

        !RDF%randField = 1.0 ! For Tests
        !RDF%randField = RDF%rang ! For Tests

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
        !INPUT
        type(RF), intent(inout) :: RDF
        !LOCAL
        integer(C_INTPTR_T) :: L, M, N
        integer(C_INTPTR_T) :: local_LastDim
        integer(C_INTPTR_T) :: local_j_offset
        real(C_DOUBLE), pointer :: data_real_2D(:,:), data_real_3D(:,:,:)
        type(C_PTR) :: cdata, plan
        integer(C_INTPTR_T) :: alloc_local
        integer :: sliceSize
        integer :: kNLocal
        integer, dimension(RDF%nDim) :: kNStepLocal
        double precision, dimension(:), allocatable :: gammaK, phiK
        integer :: i, j, k, ind
        double precision :: trashNumber

        call wLog("gen_Std_Gauss_FFT_init")
        call wLog(" ")

        L = RDF%xNStep(1)
        if(RDF%nDim >= 2) M = RDF%xNStep(2)
        if(RDF%nDim >= 3) N = RDF%xNStep(2)

        if(RDF%independent) then
            call wLog("    LOCAL")

            if(RDF%nDim == 2) then
                alloc_local = L*M
                cdata = fftw_alloc_real(alloc_local)
                call c_f_pointer(cdata, data_real_2D, [L, M])

            else if(RDF%nDim == 3) then
                alloc_local = L*M*N
                cdata = fftw_alloc_real(alloc_local)
                call c_f_pointer(cdata, data_real_3D, [L, M, N])

            else
                stop("Inside gen_Std_Gauss_FFT dimension not yet implemented for this generation method")

            end if

            call wLog("Defining kPoints and SkVec")
            call set_kPoints(RDF)
            call set_SkVec(RDF)

            !call wLog("SkVec")
            !do i = 1, size(RDF%SkVec)
            !    call wLog(RDF%SkVec(i))
            !end do

            !call gen_Std_Gauss_FFT_step2_monoproc(RDF, RDF%SkVec)

            allocate(gammaK(RDF%kNTotal))
            allocate(phik(RDF%kNTotal))

            call wLog("shape(gammaK)")
            call wLog(shape(gammaK))
            call wLog("shape(phiK)")
            call wLog(shape(phiK))
            call wLog("shape(RDF%SkVec)")
            call wLog(shape(RDF%SkVec))

            call random_number(gammaK(:))
            call random_number(phiK(:))

            call wLog("Not calculating Field")

            gammaK       = gammaK -0.5
            RDF%SkVec(:) = gammak*sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);
            !RDF%SkVec(:) =  sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);

            if(allocated(gammaK)) deallocate(gammaK)
            if(allocated(phik))   deallocate(phik)

            if(RDF%nDim == 2) then
                plan = fftw_plan_r2r(RDF%nDim, [int(M), int(L)], data_real_2D, data_real_2D, &
                                        [FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
                data_real_2D(:,:) = reshape(RDF%SkVec, [L, M])
                call fftw_execute_r2r(plan, data_real_2D, data_real_2D)
                RDF%randField(:,1) = pack(data_real_2D(1:RDF%xNStep(1), 1:RDF%xNStep(2)), .true.)
                call fftw_destroy_plan(plan)

            else if(RDF%nDim == 3) then
                plan = fftw_plan_r2r(RDF%nDim, [int(N), int(M), int(L)], data_real_3D, data_real_3D, &
                                        [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
                data_real_3D(:,:,:) = reshape(RDF%SkVec, [L, M,N])
                call fftw_execute_r2r(plan, data_real_3D, data_real_3D)
                RDF%randField(:,1) = pack(data_real_3D(1:RDF%xNStep(1), 1:RDF%xNStep(2), 1:RDF%xNStep(3)), .true.)
                call fftw_destroy_plan(plan)

            else
                stop("No FFT method implemented only for 2D and 3D cases")
            end if

            !RDF%randField(:,1) = RDF%SkVec(:) !FOR TESTS
            !RDF%randField(:,1) = RDF%kPoints(:) !FOR TESTS
            !RDF%randField(:,1) = RDF%rang!FOR TESTS

        else
            call fftw_mpi_init()
            call wLog("    GLOBAL")

            if(RDF%nDim == 2) then
                L = RDF%xNStep(1)
                M = RDF%xNStep(2)
                alloc_local = fftw_mpi_local_size_2d(M, L, RDF%comm, &
                                                     local_LastDim, local_j_offset) !FOR MPI
                cdata = fftw_alloc_real(alloc_local)
                call c_f_pointer(cdata, data_real_2D, [L, local_LastDim])

            else if(RDF%nDim == 3) then
                L = RDF%xNStep(1)
                M = RDF%xNStep(2)
                N = RDF%xNStep(3)
                alloc_local = fftw_mpi_local_size_3d(N, M, L, RDF%comm, &
                                                     local_LastDim, local_j_offset) !FOR MPI
                cdata = fftw_alloc_real(alloc_local)
                call c_f_pointer(cdata, data_real_3D, [L, M, local_LastDim])

            else
                stop("Inside gen_Std_Gauss_FFT dimension not yet implemented for this generation method")

            end if

            !Defining kInit and kEnd
            RDF%kNInit = local_j_offset + 1
            RDF%kNEnd  = RDF%kNInit + local_LastDim - 1
            call wLog("local_j_offset")
            call wLog(int(local_j_offset))
            call wLog("local_LastDim")
            call wLog(int(local_LastDim))

            call wLog("RDF%kNInit")
            call wLog(int(RDF%kNInit))
            call wLog("RDF%kNEnd")
            call wLog(int(RDF%kNEnd))

            sliceSize = 1
            if(RDF%nDim > 1) sliceSize = product(RDF%xNStep(1:RDF%nDim -1))
            RDF%kNInit = (RDF%kNInit - 1) * sliceSize + 1
            RDF%kNEnd  = RDF%kNEnd *sliceSize

            !RDF%origin  = [1, int(local_j_offset) + 1]
            !RDF%kExtent = [L , local_M]
            RDF%origin = 1
            RDF%origin(RDF%nDim) = int(local_j_offset) + 1
            RDF%kExtent = RDF%xNStep
            RDF%kExtent(RDF%nDim) = int(local_LastDim)

            call set_kPoints(RDF)
            call set_SkVec(RDF)

            !call wLog("RDF%kPoints")
            !call wLog(RDF%kPoints)
            !call wLog("RDF%SkVec")
            !call wLog(RDF%SkVec)



            kNLocal = size(RDF%kPoints,2)
            !kNStepLocal = RDF%kNStep
            !kNStepLocal(RDF%nDim) = kNLocal/sliceSize

            call wLog("RDF%origin")
            call wLog(RDF%origin)
            call wLog("RDF%kExtent")
            call wLog(RDF%kExtent)
            call wLog("kNLocal")
            call wLog(kNLocal)
            allocate(gammaK(kNLocal))
            allocate(phik(kNLocal))

            !Putting away the random numbers from others k (that are in others procs)
            do i = 1, RDF%kNInit-1
                call random_number(trashNumber)
            end do
            call random_number(gammaK(:))
            do i = RDF%kNEnd+1, product(RDF%kNStep)
                call random_number(trashNumber)
            end do
            do i = 1, RDF%kNInit-1
                call random_number(trashNumber)
            end do
            call random_number(phiK(:))

            call wLog("shape(gammaK)")
            call wLog(shape(gammaK))
            call wLog("shape(phiK)")
            call wLog(shape(phiK))
            call wLog("shape(RDF%SkVec)")
            call wLog(shape(RDF%SkVec))

            gammaK       = gammaK -0.5
            RDF%SkVec(:) =  gammak*sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);
            !RDF%SkVec(:) =  sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);

            if(allocated(gammaK)) deallocate(gammaK)
            if(allocated(phik))   deallocate(phik)

            !cdata = fftw_alloc_real(alloc_local)
            !call c_f_pointer(cdata, data_real_2D, [L,local_M])

            if(RDF%nDim == 2) then
                plan = fftw_mpi_plan_r2r(RDF%nDim, [M, L], data_real_2D, data_real_2D, &
                                         RDF%comm, [FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
                data_real_2D(:,:) = reshape(RDF%SkVec, [L, local_LastDim])
                call fftw_mpi_execute_r2r(plan, data_real_2D, data_real_2D)
                RDF%randField(:,1) = pack(data_real_2D(1:L,1:local_LastDim), .true.)

            else if(RDF%nDim == 3) then
                plan = fftw_mpi_plan_r2r(RDF%nDim, [N, M, L], data_real_3D, data_real_3D, &
                                         RDF%comm, [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
                data_real_3D(:,:,:) = reshape(RDF%SkVec, [L, M, local_LastDim])
                call fftw_mpi_execute_r2r(plan, data_real_3D, data_real_3D)
                RDF%randField(:,1) = pack(data_real_3D(1:L, 1:M, 1:local_LastDim), .true.)
            end if

            !call wLog("L")
            !call wLog(int(L))
            !call wLog("local_LastDim")
            !call wLog(int(local_LastDim))
            call wLog("shape(RDF%SkVec)")
            call wLog(shape(RDF%SkVec))
            call wLog("shape(RDF%randField)")
            call wLog(shape(RDF%randField))


            call fftw_destroy_plan(plan)
            call fftw_free(cdata)

            call wLog("L")
            call wLog(int(L))
            call wLog("M")
            call wLog(int(M))
            call wLog("local_LastDim")
            call wLog(int(local_LastDim))
            call wLog("local_j_offset")
            call wLog(int(local_j_offset))
            call wLog("RDF%kNInit")
            call wLog(int(RDF%kNInit))
            call wLog("RDF%kNEnd")
            call wLog(int(RDF%kNEnd))


            !call gen_Std_Gauss_FFT_step2_multiproc(RDF, RDF%SkVec)
        end if

        if(allocated(gammaK)) deallocate(gammaK)
        if(allocated(phik)) deallocate(phik)


        !call gen_Std_Gauss_FFT_step2_monoproc(RDF, RDF%SkVec)
        !call gen_Std_Gauss_FFT_step2_multiproc(RDF, RDF%SkVec)

    end subroutine gen_Std_Gauss_FFT


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss_FFT_step2_monoproc(RDF, SkVec)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        double precision, dimension(RDF%kNTotal), target :: SkVec
        !double precision, dimension(RDF%xNTotal), target :: randFieldVec

        !START EXEMPLE
        integer(C_INTPTR_T) :: L, M, N
        !integer(4), parameter :: Lb = 3
        !integer(4), parameter :: Mb = 4
        type(C_PTR) :: plan, cdata_complex, cdata_real
        complex(C_DOUBLE_COMPLEX), pointer :: data_complex(:,:)
        real(C_DOUBLE), pointer :: data_real_2D(:,:), data_real_3D(:,:,:)
        integer(C_INTPTR_T) :: i, j, alloc_local, local_M, local_j_offset
        double precision, dimension(RDF%kNTotal) :: gammaK, phiK
        double precision, dimension(12), target :: vecTest=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]


        call random_number(gammaK(:))
        call random_number(phiK(:))
        gammaK   = gammaK -0.5
        SkVec(:) =  gammak*sqrt(SkVec)*cos(2.0D0*PI*phik);

        call wLog("Inside gen_Std_Gauss_FFT")

        !get local data size and allocate (note dimension reversal)
        !alloc_local = fftw_mpi_local_size_2d(M, L, MPI_COMM_WORLD, &
        !                                     local_M, local_j_offset) !FOR MPI
        !cdata_complex = fftw_alloc_complex(alloc_local)

        cdata_real = c_loc(SkVec)
        !cdata_real= fftw_alloc_real(L*M)


        if(RDF%nDim == 1) then
            stop("No FFT method implemented for 1D cases")
            L = RDF%kNStep(1)

        else if(RDF%nDim == 2) then

            L = RDF%kNStep(1)
            M = RDF%kNStep(2)
            call c_f_pointer(cdata_real, data_real_2D, [L,M])
            !plan = fftw_plan_r2r_2d(int(M), int(L), data_real_2D, data_real_2D, FFTW_REDFT01, &
            !                        FFTW_REDFT01, FFTW_ESTIMATE)
            plan = fftw_plan_r2r(RDF%nDim, [int(M), int(L)], data_real_2D, data_real_2D, &
                                    [FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
            !FFTW_BACKWARD
            call fftw_execute_r2r(plan, data_real_2D, data_real_2D)
            RDF%randField(:,1) = pack(data_real_2D(1:RDF%xNStep(1), 1:RDF%xNStep(2)), .true.)
            call fftw_destroy_plan(plan)

        else if(RDF%nDim == 3) then

            L = RDF%kNStep(1)
            M = RDF%kNStep(2)
            N = RDF%kNStep(3)
            call c_f_pointer(cdata_real, data_real_3D, [L,M,N])
            !plan = fftw_plan_r2r_3d(int(N), int(M), int(L), data_real_3D, data_real_3D, FFTW_REDFT01, &
            !                        FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
            plan = fftw_plan_r2r(RDF%nDim, [int(N), int(M), int(L)], data_real_3D, data_real_3D, &
                                    [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
            call fftw_execute_r2r(plan, data_real_3D, data_real_3D)
            RDF%randField(:,1) = pack(data_real_3D(1:RDF%xNStep(1), 1:RDF%xNStep(2), 1:RDF%xNStep(3)), .true.)
            call fftw_destroy_plan(plan)

        else
            stop("No FFT method implemented only for 2D and 3D cases")
        end if


    end subroutine gen_Std_Gauss_FFT_step2_monoproc

    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    subroutine allocate_randField(RDF, randField)
        !INPUT AND OUTPUT
        type(RF)   :: RDF
        double precision, dimension(:,:), allocatable, target :: randField

        if(allocated(randField)) then
            if(.not.(size(randField,1) == RDF%xNTotal .and. size(randField,1) == RDF%Nmc)) then
                nullify(RDF%randField)
                deallocate(randField)
            end if
        end if

        if(.not.allocated(randField)) then
            allocate(randField(RDF%xNTotal, RDF%Nmc))
            RDF%randField => randField
            call wLog("Inside allocate_randField")
            call wLog("RDF%xNStep = ")
            call wLog(RDF%xNStep)
            if(RDF%nDim == 2) RDF%RF_2D(1:RDF%xNStep(1),1:RDF%xNStep(2)) => randField
            if(RDF%nDim == 3) RDF%RF_3D(1:RDF%xNStep(1),1:RDF%xNStep(2),1:RDF%xNStep(3)) => randField
        end if

    end subroutine allocate_randField


!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine gen_Std_Gauss_FFT(RDF)
!
!        implicit none
!
!        !INPUT OUTPUT
!        type(RF), intent(inout) :: RDF
!
!        !write(get_fileId(),*) "Inside FFT"
!
!        !write(get_fileId(),*) "Defining kPoints and SkVec"
!        call set_kPoints(RDF)
!        call set_SkVec(RDF)
!
!        !write(get_fileId(),*) "Step 2"
!        !if(RDF%independent) then
!            call gen_Std_Gauss_FFT_step2_indep(RDF, RDF%SkVec, RDF%randField(:,1))
!        !else
!        !    call gen_Std_Gauss_FFT_step2_glob(RDF, RDF%SkVec, RDF%randField(:,1))
!        !end if
!    end subroutine gen_Std_Gauss_FFT

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine gen_Std_Gauss_FFT_step2_glob(RDF, SkVec, randFieldVec)
!
!        implicit none
!
!        !INPUT OUTPUT
!        type(RF), intent(inout) :: RDF
!        double precision, dimension(:), target :: SkVec
!        double precision, dimension(:), target :: randFieldVec
!
!        !POINTERS TO INPUT
!        double precision, dimension(RDF%kNTotal) :: gammaK, phiK
!        double precision, dimension(:,:)  , pointer :: SkVec_2D
!        double precision, dimension(:,:,:), pointer :: SkVec_3D
!        double precision, dimension(:,:)  , pointer :: RF_2D
!        double precision, dimension(:,:,:)  , pointer :: RF_3D
!
!        !LOCAL
!        integer, dimension(RDF%nDim) :: M;
!        double complex  , dimension(:,:)  , pointer :: SkSym_2D
!        double complex  , dimension(:,:,:), pointer :: SkSym_3D
!        double complex  , dimension(product(2*RDF%xNStep-2)), target :: SkSym
!        double precision, dimension(product(2*RDF%xNStep-2)), target :: RFSym
!        double precision, dimension(:,:)  , pointer :: RFSym_2D
!        double precision, dimension(:,:,:)  , pointer :: RFSym_3D
!        integer, dimension(RDF%nDim) :: kS, kE, kSc, kEc, kCore;
!        integer :: pos, i,j
!        integer*8 plan, planTest1, planTest2
!        double precision :: sizesProd

        !START EXEMPLE
!        integer(C_INTPTR_T), parameter :: L = 3
!        integer(C_INTPTR_T), parameter :: M = 4
!        type(C_PTR) :: plan, cdata
!        complex(C_DOUBLE_COMPLEX), pointer :: data(:,:)
!        integer(C_INTPTR_T) :: i, j, alloc_local, local_M, local_j_offset

!        call wLog("Inside FFT Global")

!        !get local data size and allocate (note dimension reversal)
!        alloc_local = fftw_mpi_local_size_2d(M, L, MPI_COMM_WORLD, &
!            local_M, local_j_offset)
!        cdata = fftw_alloc_complex(alloc_local)
!        call c_f_pointer(cdata, data, [L,local_M])
!
!        !   create MPI plan for in-place forward DFT (note dimension reversal)
!        plan = fftw_mpi_plan_dft_2d(M, L, data, data, MPI_COMM_WORLD, &
!            FFTW_FORWARD, FFTW_ESTIMATE)
!
!        ! initialize data to some function my_function(i,j)
!        do j = 1, local_M
!            do i = 1, L
!                data(i, j) = i+ j*L
!            end do
!        end do
!
!        ! compute transform (as many times as desired)
!        call fftw_mpi_execute_dft(plan, data, data)
!
!        call fftw_destroy_plan(plan)
!        call fftw_free(cdata)
!
!        !END EXEMPLE
!
!
!        M = 2*RDF%xNStep-2
!
!        if(RDF%nDim == 2) then
!            SkVec_2D(1:RDF%kNStep(1),1:RDF%kNStep(2)) => SkVec
!            RF_2D(1:RDF%xNStep(1),1:RDF%xNStep(2))    => randFieldVec
!            SkSym_2D(1:M(1),1:M(2)) => SkSym
!            RFSym_2D(1:M(1),1:M(2)) => RFSym
!        else if(RDF%nDim == 3) then
!            SkVec_3D(1:RDF%kNStep(1),1:RDF%kNStep(2),1:RDF%kNStep(3)) => SkVec
!            RF_3D(1:RDF%xNStep(1),1:RDF%xNStep(2),1:RDF%xNStep(3))    => randFieldVec
!            SkSym_3D(1:M(1),1:M(2),1:M(3)) => SkSym
!            RFSym_3D(1:M(1),1:M(2),1:M(3)) => RFSym
!        end if
!
!        !START IMPLEMENTATION 2D---------------------------------------
!
!        !VISUALIZATION of SkVec_2D BEFORE
!        !RF_2D = SkVec_2D
!        call wLog("Creating Field by FFT")
!
!        call random_number(gammaK(:))
!        call random_number(phiK(:))
!        gammaK = gammaK -0.5
!
!        SkVec(:) =  gammak*sqrt(SkVec)*cos(2.0D0*PI*phik);
!
!        kCore =  RDF%kNStep
!        kS  = 2
!        kE  = kCore - 1
!        kSc = kCore + 1
!        kEc = M
!
!        !Core
!        SkSym_2D(1:kCore(1), 1:kCore(2)) = SkVec_2D(1:kCore(1), 1:kCore(2))
!        !Sides
!        SkSym_2D(kSc(1):kEc(1)  , 1:kCore(2)) = SkVec_2D(kS(1):kE(1),:)
!        SkSym_2D(1:kCore(1), kSc(2):kEc(2)  ) = SkVec_2D(:,kS(2):kE(2))
!        !Diagonal
!        SkSym_2D(kSc(1):kEc(1)  , kSc(2):kEc(2)  ) = SkVec_2D(kS(1):kE(1), kS(2):kE(2))
!        !Hermitian Conjugate
!        SkSym_2D(kSc(1):kEc(1), :) = -conjg(SkSym_2D(kEc(1):kSc(1):-1, :))
!        SkSym_2D(:, kSc(2):kEc(2)) = -conjg(SkSym_2D(:, kEc(2):kSc(2):-1))
!
!
!        call dfftw_plan_dft_2d(planTest2, size(SkSym_2D,1), size(SkSym_2D,2), SkSym_2D, SkSym_2D, &
!                               FFTW_BACKWARD, FFTW_ESTIMATE)
!        call dfftw_execute_dft(planTest2, SkSym_2D, SkSym_2D)
!        call dfftw_destroy_plan(planTest2)
!        !write(get_fileId(),*) "After FFT"
!
!        !write(get_fileId(),*) "shape(RF_2D)    = ", shape(RF_2D)
!        !write(get_fileId(),*) "shape(SkSym_2D) = ", shape(SkSym_2D)
!
!        !VISUALIZATION of SkSym_2D
!        RF_2D = real(SkSym_2D(1:kCore(1),1:kCore(2)))/dble(product(shape(SkSym_2D)))
!
!        if(associated(RF_2D))    nullify(RF_2D)
!        if(associated(RF_3D))    nullify(RF_3D)
!        if(associated(RFSym_2D)) nullify(RFSym_2D)
!        if(associated(RFSym_3D)) nullify(RFSym_3D)
!        if(associated(SkVec_2D)) nullify(SkVec_2D)
!        if(associated(SkVec_3D)) nullify(SkVec_3D)
!        if(associated(SkSym_2D)) nullify(SkSym_2D)
!        if(associated(SkSym_3D)) nullify(SkSym_3D)
!
!    end subroutine gen_Std_Gauss_FFT_step2_glob

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
!    subroutine gen_Std_Gauss_FFT_step2_indep(RDF, SkVec, randFieldVec)
!
!        implicit none
!
!        !INPUT OUTPUT
!        type(RF), intent(inout) :: RDF
!        double precision, dimension(:), target :: SkVec
!        double precision, dimension(:), target :: randFieldVec
!
!        !POINTERS TO INPUT
!        double precision, dimension(RDF%kNTotal) :: gammaK, phiK
!
!        SkVec(:) =  gammak*sqrt(SkVec)*cos(2.0D0*PI*phik);
!        double precision, dimension(:,:)  , pointer :: SkVec_2D
!        double precision, dimension(:,:,:), pointer :: SkVec_3D
!        double precision, dimension(:,:)  , pointer :: RF_2D
!        double precision, dimension(:,:,:)  , pointer :: RF_3D
!
!        !LOCAL
!        !integer, dimension(RDF%nDim) :: M;
!        double complex  , dimension(:,:)  , pointer :: SkSym_2D
!        double complex  , dimension(:,:,:), pointer :: SkSym_3D
!        double complex  , dimension(product(2*RDF%xNStep-2)), target :: SkSym
!        double precision, dimension(product(2*RDF%xNStep-2)), target :: RFSym
!        double precision, dimension(:,:)  , pointer :: RFSym_2D
!        double precision, dimension(:,:,:)  , pointer :: RFSym_3D
!        integer, dimension(RDF%nDim) :: kS, kE, kSc, kEc, kCore;
!        integer :: pos, i,j
!        integer*8 plan, planTest1, planTest2
!        double precision :: sizesProd
!
!        !TESTS
!        integer, parameter :: N = 2000, N_h = 1000
!        double complex,   dimension(N) :: testV1, testV2
!        double complex  , dimension(N) :: testComplxSym
!        double precision, dimension(N) :: testReal
!        double complex  , dimension(N,N) :: testComplxSym2D
!        double complex  , dimension(N,N) :: testComplxSym2D_B
!        double precision, dimension(N,N) :: testReal2D
!        double complex,   dimension(N,N) :: arr
!
!        !START C EXEMPLE
!        integer(C_INTPTR_T) :: L
!        integer(C_INTPTR_T) :: M
!        integer(C_INTPTR_T) :: N
!        type(C_PTR) :: plan, cdata
!        complex(C_DOUBLE_COMPLEX), pointer :: data(:,:)
!        integer(C_INTPTR_T) :: i, j, alloc_local, local_M, local_j_offset
!
!
!        !write(get_fileId(),*) "RDF%xNStep   = ", RDF%xNStep
!        !write(get_fileId(),*) "shape(SkVec) = ", shape(SkVec)
!        !write(get_fileId(),*) "RDF%kNStep   = ", RDF%kNStep
!
!        !M = 2*RDF%xNStep-2
!
!        if(RDF%nDim == 2) then
!            SkVec_2D(1:RDF%kNStep(1),1:RDF%kNStep(2)) => SkVec
!            RF_2D(1:RDF%xNStep(1),1:RDF%xNStep(2))    => randFieldVec
!            L = RDF%kNStep(1)
!            M = RDF%kNStep(2)
!            !SkSym_2D(1:M(1),1:M(2)) => SkSym
!            !RFSym_2D(1:M(1),1:M(2)) => RFSym
!        else if(RDF%nDim == 3) then
!            SkVec_3D(1:RDF%kNStep(1),1:RDF%kNStep(2),1:RDF%kNStep(3)) => SkVec
!            RF_3D(1:RDF%xNStep(1),1:RDF%xNStep(2),1:RDF%xNStep(3))    => randFieldVec
!            !SkSym_3D(1:M(1),1:M(2),1:M(3)) => SkSym
!            !RFSym_3D(1:M(1),1:M(2),1:M(3)) => RFSym
!        end if

        !START TEST 1D---------------------------------------------------
!        testV1(:) = 0.0D0
!        testV2(:) = 0.0D0
!        testV1(686:1314) = [(sin(dble(i)/100), i =0, 628)]
!        testReal(:) = real(testV1(1:N))
!
!
!        open (unit = 18 , file = "FFT1D/beforeFFT", action = 'write')
!        open (unit = 19 , file = "FFT1D/afterFFT", action = 'write')
!        open (unit = 20 , file = "FFT1D/afteriFFT", action = 'write')
!        open (unit = 21 , file = "FFT1D/beforeFFT_2", action = 'write')
!        open (unit = 22 , file = "FFT1D/afterFFT_2", action = 'write')
!        open (unit = 23 , file = "FFT1D/afteriFFT_2", action = 'write')
!
!        write(*,*) "BEFORE FFT testV1"
!        write(18,"(2F15.5)") testV1
!        write(21,"(F15.5)") testReal
!
!        call dfftw_plan_dft_1d(planTest2, N, &
!                               testV1, testV2, FFTW_FORWARD, FFTW_ESTIMATE)
!        call dfftw_execute(planTest2)
!
!        testComplxSym = 0.0D0
!        testComplxSym(1:N_h) = testV2(1:N_h)
!
!        write(*,*) "AFTER FFT testV2"
!        write(19,"(2F15.5)") testV2
!        write(22,"(2F15.5)") testComplxSym
!
!
!        call dfftw_plan_dft_1d(planTest2, N, &
!                               testV2, testV1, FFTW_BACKWARD, FFTW_ESTIMATE)
!        call dfftw_execute(planTest2)
!
!        call dfftw_plan_dft_c2r(planTest1, 1, [N], &
!                                testComplxSym, testReal, FFTW_ESTIMATE);
!        call dfftw_execute(planTest1)
!
!        write(*,*) "AFTER iFFT testV1"
!        write(20,"(2F15.5)") testV1
!        write(23,"(F15.5)") testReal
!
!        call dfftw_destroy_plan(planTest2)
!        call dfftw_destroy_plan(planTest1)
!
!        close (18)
!        close (19)
!        close (20)
!        close (21)
!        close (22)
!        close (23)
        !END TEST 1D---------------------------------------------------------------

        !START TEST 2D---------------------------------------------------

!        testV1(:) = 0.0D0
!        testV2(:) = 0.0D0
!        testV1(686:1314) = [(sin(dble(i)/100), i =0, 628)]
!        testReal(:) = real(testV1(1:N))
!
!        do j = 1, RDF%xNStep(2)
!            do i = 1, RDF%xNStep(1)
!                !RF_2D(i,j) =  testReal(i)*testReal(j)
!                testReal2D(i,j) =  testReal(i)*testReal(j)
!            end do
!        end do
!
!        arr = testReal2D
!
!        !RF_2D = abs(arr)
!
!        call dfftw_plan_dft_2d(planTest2, N,N, arr,arr, &
!                                   FFTW_FORWARD, FFTW_ESTIMATE)
!        call dfftw_execute_dft(planTest2, arr, arr)
!        call dfftw_destroy_plan(planTest2)
!
!        !RF_2D = abs(arr)
!
!        !Making Symmetric
!        !write(get_fileId(),*) "Symmetrization"
!        !Symmetrization
!        kS  = 2
!        kE  = N/2-1
!        kSc = N/2 +2;
!        kEc = N;
!
!        !Core
!        testComplxSym2D_B(1:N/2+1,1:N/2+1) = arr(1:N/2+1,1:N/2+1)
!        !Sides
!        testComplxSym2D_B(1:N/2+1,kSc(2):kEc(2)) = arr(1:N/2+1,kS(2):kE(2))
!        testComplxSym2D_B(kSc(1):kEc(1),1:N/2+1) = arr(kS(1):kE(1),1:N/2+1)
!        !Corner
!        testComplxSym2D_B(kSc(1):kEc(1),kSc(2):kEc(2)) = arr(kS(1):kE(1),kS(2):kE(2))
!        !Hermitian
!        testComplxSym2D_B(:,kSc(2):kEc(2)) = -conjg(testComplxSym2D_B(:,kEc(2):kSc(2):-1))
!        testComplxSym2D_B(kSc(1):kEc(1),:) = -conjg(testComplxSym2D_B(kEc(1):kSc(1):-1,:))
!
!        !RF_2D = abs(testComplxSym2D_B)
!
!        !call dfftw_plan_dft_2d(planTest2, N,N, arr,arr, &
!        !                       FFTW_BACKWARD, FFTW_ESTIMATE)
!        !call dfftw_execute_dft(planTest2, arr, arr)
!        !call dfftw_destroy_plan(planTest2)
!        !RF_2D = abs(arr)
!
!        call dfftw_plan_dft_2d(planTest2, N,N, testComplxSym2D_B, testComplxSym2D_B, &
!                               FFTW_BACKWARD, FFTW_ESTIMATE)
!        call dfftw_execute_dft(planTest2, testComplxSym2D_B, testComplxSym2D_B)
!        call dfftw_destroy_plan(planTest2)
!        RF_2D = real(testComplxSym2D_B)/(N*N)
!
!            !RF_2D = testReal2D!Information To Plot on Paraview
!
!            !call dfftw_plan_dft_r2c(planTest1, RDF%nDim, RDF%xNStep, &
!            !                        testReal2D, testComplxSym2D, FFTW_ESTIMATE);
!            !call dfftw_execute(planTest1)
!            !call dfftw_destroy_plan(planTest1)
!
!            !testComplxSym2D = testComplxSym2D/dble(product(RDF%xNStep))
!
!            !RF_2D = abs(testComplxSym2D)

    !        !Remaking symmetric (only 1/4)
    !        !Core
    !        testComplxSym2D_B = 0.0D0
    !        testComplxSym2D_B(1:N/2+1, 1:N/2+1)= testComplxSym2D(1:N/2+1, 1:N/2+1)
    !
    !
    !        !call dfftw_plan_dft_c2r(planTest1, RDF%nDim, RDF%xNStep, &
    !        !                        testComplxSym2D_B, testReal2D, FFTW_ESTIMATE);
    !        !call dfftw_execute(planTest1)
    !
    !        call dfftw_plan_dft_c2r(planTest1, RDF%nDim, RDF%xNStep, &
    !                                testComplxSym2D, testReal2D, FFTW_ESTIMATE);
    !        call dfftw_execute(planTest1)
    !
    !        call dfftw_destroy_plan(planTest1)
    !
    !        write(*,*) "Calculating RF_2D"
    !        RF_2D = abs(testComplxSym2D)
            !RF_2D = real(testComplxSym2D_B(1::2, 1::2))
            !RF_2D = real(testComplxSym2D(1::2, 1::2)-testComplxSym2D_B(1::2, 1::2))
            !RF_2D = abs(testComplxSym2D)!Information To Plot on Paraview
            !RF_2D = testReal2D!Information To Plot on Paraview

        !END TEST 2D---------------------------------------------------

        !START IMPLEMENTATION 2D---------------------------------------

!        !VISUALIZATION of SkVec_2D BEFORE
!        !RF_2D = SkVec_2D
!        call wLog("Creating Field by FFT")
!
!        call random_number(gammaK(:))
!        call random_number(phiK(:))
!        gammaK = gammaK -0.5
!
!        SkVec(:) =  gammak*sqrt(SkVec)*cos(2.0D0*PI*phik);
!
!
!        !VISUALIZATION of SkVec_2D AFTER
!        !RF_2D = SkVec_2D
!
!        !write(get_fileId(),*) "Symmetrization"
!        kCore =  RDF%kNStep
!        kS  = 2
!        kE  = kCore - 1
!        kSc = kCore + 1
!        kEc = M;
!
!        !Core
!        SkSym_2D(1:kCore(1), 1:kCore(2)) = SkVec_2D(1:kCore(1), 1:kCore(2))
!        !Sides
!        SkSym_2D(kSc(1):kEc(1)  , 1:kCore(2)) = SkVec_2D(kS(1):kE(1),:)
!        SkSym_2D(1:kCore(1), kSc(2):kEc(2)  ) = SkVec_2D(:,kS(2):kE(2))
!        !Diagonal
!        SkSym_2D(kSc(1):kEc(1)  , kSc(2):kEc(2)  ) = SkVec_2D(kS(1):kE(1), kS(2):kE(2))
!        !Hermitian Conjugate
!        SkSym_2D(kSc(1):kEc(1), :) = -conjg(SkSym_2D(kEc(1):kSc(1):-1, :))
!        SkSym_2D(:, kSc(2):kEc(2)) = -conjg(SkSym_2D(:, kEc(2):kSc(2):-1))
!
!        !VISUALIZATION of SkSym_2D
!        !RF_2D = SkSym_2D
!
!        !write(get_fileId(),*) "Before FFT"
!        call dfftw_plan_dft_2d(planTest2, size(SkSym_2D,1), size(SkSym_2D,2), SkSym_2D, SkSym_2D, &
!                               FFTW_BACKWARD, FFTW_ESTIMATE)
!        call dfftw_execute_dft(planTest2, SkSym_2D, SkSym_2D)
!        call dfftw_destroy_plan(planTest2)
!
!        !VISUALIZATION of SkSym_2D
!        RF_2D = real(SkSym_2D(1:kCore(1),1:kCore(2)))/dble(product(shape(SkSym_2D)))
!
!        !VISUALIZATION of SkSym (1:2)
!        !RF_2D = 5.0
!        !RF_2D(1,1)   = real(SkSym_2D(1,1))
!        !RF_2D(2:,1)  = real(SkSym_2D(2:M(1)-1:2,1))
!        !RF_2D(1,2:)  = real(SkSym_2D(1,2:M(2)-1:2))
!
!        !RF_2D(2:,2:) = real(SkSym_2D(2:M(1)-1:2, 2:M(2)-1:2))
!
!        !RF_2D(kCore(1),kCore(2))   = real(SkSym_2D(kCore(1),kCore(2)))
!        !RF_2D(2:,kCore(1))  = real(SkSym_2D(2:M(1)-1:2,kCore(1)))
!        !RF_2D(kCore(1),2:)  = real(SkSym_2D(kCore(1),2:M(2)-1:2))
!
!
!        !VISUALIZATION of the result
!        !RF_2D = aimag(SkSym_2D(1:RDF%kNStep(1), 1:RDF%kNStep(2)))/dble(product(shape(SkSym_2D)))
!        !RF_2D = real(SkSym_2D(1:RDF%kNStep(1), 1:RDF%kNStep(2)))/dble(product(shape(SkSym_2D)))

        !END IMPLEMENTATION 2D--------------------------------

!        !START 2D SECOND IMPLEMENTATION---------------------
!        call wLog("Creating Field by FFT (Imlpementation 2)")
!
!        call random_number(gammaK(:))
!        call random_number(phiK(:))
!        gammaK = gammaK -0.5
!
!        SkVec(:) =  gammak*sqrt(SkVec)*cos(2.0D0*PI*phik);
!
!
!        !call dfftw_plan_dft_r2r_2d(planTest2, size(SkVec_2D,1), size(SkVec_2D,2), SkVec_2D, SkVec_2D, &
!        !                       FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
!        call dfftw_plan_dft_r2r(planTest2, RDF%nDim, shape(SkVec_2D), SkVec_2D, SkVec_2D, &
!                                [FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
!        call dfftw_execute_dft(planTest2, SkVec_2D, SkVec_2D)
!        call dfftw_destroy_plan(planTest2)
!
!        RF_2D = SkVec_2D


        !END 2D SECOND IMPLEMENTATION-----------------------

        !START TEST 3D---------------------------------------------------

!        !write(get_fileId(),*) "Symmetrization"
!        !kS  = 2
!        !kE  = RDF%kNStep-1
!        !kSc = RDF%kNStep + 1;
!        !kEc = 2*RDF%kNStep -2;
!
!        !write(get_fileId(),*) " RDF%kNStep = ", RDF%kNStep
!        !write(get_fileId(),*) " RDF%xNStep = ", RDF%xNStep
!        !write(get_fileId(),*) " shape(SkVec_3D) = ", shape(SkVec_3D)
!        !write(get_fileId(),*) " shape(SkSym_3D) = ", shape(SkSym_3D)
!        !write(get_fileId(),*) " shape(RF_3D)    = ", shape(RF_3D)
!        do i = 1, 5
!        !write(get_fileId(),*) " i = ", i
!        !write(get_fileId(),*) " SkVec_3D(i,i,i) = ", SkVec_3D(i,i,i)
!        end do
!        i = 1
!        !write(get_fileId(),*) " i = ", i
!        !write(get_fileId(),*) " SkVec_3D(i,i,i) = ", SkVec_3D(i,i,i)
!        !write(get_fileId(),*) " SkSym_3D(i,i,i) = ", SkSym_3D(i,i,i)
!        !!write(get_fileId(),*) " kS  = ", kS
!        !!write(get_fileId(),*) " kE  = ", kE
!        !!write(get_fileId(),*) " kSc = ", kSc
!        !!write(get_fileId(),*) " kEc = ", kEc
!
!        SkVec_3D = sqrt(SkVec_3D)
!
!        !write(get_fileId(),*) "Copy"
!        SkSym_3D = 0.0D0
!        SkSym_3D(1:RDF%kNStep(1), 1:RDF%kNStep(2), 1:RDF%kNStep(3)) = SkVec_3D(:,:,:)
!
!        !!write(get_fileId(),*) " SkSym_3D = ", SkSym_3D
!
!        !SkSym_3D(kSc(1):kEc(1)  , 1:RDF%kNStep(2), 1:RDF%kNStep(3)) = SkVec_3D(kS(1):kE(1),:,:)
!        !SkSym_3D(1:RDF%kNStep(1), kSc(2):kEc(2)  , 1:RDF%kNStep(3)) = SkVec_3D(:,kS(2):kE(2),:)
!        !SkSym_3D(1:RDF%kNStep(1), 1:RDF%kNStep(2), kSc(3):kEc(3)  ) = SkVec_3D(:,:,kS(3):kE(3))
!
!        !SkSym_3D(kSc(1):kEc(1)  , kSc(2):kEc(2)  , 1:RDF%kNStep(3)) = SkVec_3D(kS(1):kE(1), kS(2):kE(2), :          ) ![-1,-1, 1]
!        !SkSym_3D(1:RDF%kNStep(1), kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = SkVec_3D(:          , kS(2):kE(2), kS(3):kE(3)) ![ 1,-1,-1]
!        !SkSym_3D(kSc(1):kEc(1)  , 1:RDF%kNStep(2), kSc(3):kEc(3)  ) = SkVec_3D(kS(1):kE(1), :          , kS(3):kE(3)) ![-1, 1,-1]
!
!        !SkSym_3D(kSc(1):kEc(1)  , kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = SkVec_3D(kS(1):kE(1), kS(2):kE(2), kS(3):kE(3)) ![-1,-1,-1]
!
!        !!write(get_fileId(),*) "Hermitian Conjugate"
!        !SkSym_3D(kSc(1):kEc(1), :, :) = -conjg(SkSym_3D(kEc(1):kSc(1):-1, :, :))
!        !SkSym_3D(:, kSc(2):kEc(2), :) = -conjg(SkSym_3D(:, kEc(2):kSc(2):-1, :))
!        !SkSym_3D(:, :, kSc(3):kEc(3)) = -conjg(SkSym_3D(:, :, kEc(3):kSc(3):-1))
!
!        do i = 1, 5
!        !write(get_fileId(),*) " i = ", i
!        !write(get_fileId(),*) " SkSym_3D(i,i,i) = ", SkSym_3D(i,i,i)
!        end do
!
!        do i = 37, 40
!        !write(get_fileId(),*) " i = ", i
!        !write(get_fileId(),*) " SkSym_3D(i,i,i) = ", SkSym_3D(i,i,i)
!        end do
!
!        !write(get_fileId(),*) "Executing FFT"
!        call dfftw_plan_dft_c2r(planTest1, RDF%nDim, shape(SkSym_3D), &
!                                SkSym_3D, RFSym_3D, FFTW_ESTIMATE);
!        call dfftw_execute(planTest1)
!        call dfftw_destroy_plan(planTest1)
!
!        do i = 1, 5
!        !write(get_fileId(),*) " i = ", i
!        !write(get_fileId(),*) " RFSym_3D(i,i,i) = ", RFSym_3D(i,i,i)
!        end do
!
!        do i = 37, 40
!        !write(get_fileId(),*) " i = ", i
!        !write(get_fileId(),*) " RFSym_3D(i,i,i) = ", RFSym_3D(i,i,i)
!        end do
!
!        RF_3D = RFSym_3D(1:RDF%xNStep(1),1:RDF%xNStep(2),1:RDF%xNStep(3)) ! To Plot on Paraview
!
!        !VISUALIZATION of SkSym
!        RF_3D = 5.0
!        RF_3D(1,1,1) = real(SkSym_3D(1,1,1))
!        RF_3D(1,1,2:) = real(SkSym_3D(1,1,2::2))
!        RF_3D(2:,1,1) = real(SkSym_3D(2::2,1,1))
!        RF_3D(1,2:,1) = real(SkSym_3D(1,2::2,1))
!        RF_3D(2:,2:,1 ) = real(SkSym_3D(2::2, 2::2, 1))
!        RF_3D(1 ,2:,2:) = real(SkSym_3D(1, 2::2, 2::2))
!        RF_3D(2:,1 ,2:) = real(SkSym_3D(2::2, 1, 2::2))
!        RF_3D(2:,2:,2:) = real(SkSym_3D(2::2, 2::2, 2::2))

        !END TEST 3D---------------------------------------------------





        !write(*,*) "AFTER FFT testV1 = ", testV1



!        call init_random_seed(RDF%seed)
!

!        call dispCarvalhol(RDF%kMax, "RDF%kMax")
!        call dispCarvalhol(RDF%kNStep, "RDF%kNStep")
!        call dispCarvalhol(RDF%kDelta, "RDF%kDelta")

!        !write(get_fileId(),*) "Defining kPoints and SkVec"
!        call set_kPoints(RDF)
!        call set_SkVec(RDF)
!
!        allocate(gammak (RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
!        allocate(phik   (RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
!        allocate(Dk_base(RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
!        !allocate(Dk     (0:2*RDF%kNStep(1)-2-1, 0:2*RDF%kNStep(2)-2-1, 0:2*RDF%kNStep(3)-2-1))
!        allocate(realOut(RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
!
!        call random_number(gammak(:,:,:))
!        call random_number(phik(:,:,:))
!
!        !Dk_base = product(xMaxGlob - xMinGlob) * gammak*sqrt(RDF%Sk3D)*exp(2*pi*phik);
!        if(RDF%independent) then
!            sizesProd = product(RDF%xMaxExt - RDF%xMinExt)
!        else
!            sizesProd = product(RDF%xMaxGlob - RDF%xMinGlob)
!        end if
!
!        Dk_base(:,:,:) = sqrt(sizesProd*SkVec_3D);








!        i = 1
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)
!        i = 2
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)
!        i = 3
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)
!        i = 4
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)
!        i = 5
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)


!        !write(get_fileId(),*) "Symmetrization"
!        !Symmetrization
!        kS  = 2
!        kE  = RDF%kNStep-1
!        kSc = RDF%kNStep + 1;
!        kEc = 2*RDF%kNStep -2;
!
!        !write(get_fileId(),*) " RDF%kNStep = ", RDF%kNStep
!        !write(get_fileId(),*) " shape(Dk_base) = ", shape(Dk_base)
!        !write(get_fileId(),*) " shape(Dk) = ", shape(Dk)
!        !write(get_fileId(),*) " kS  = ", kS
!        !write(get_fileId(),*) " kE  = ", kE
!        !write(get_fileId(),*) " kSc = ", kSc
!        !write(get_fileId(),*) " kEc = ", kEc
!
!        !write(get_fileId(),*) "Copy"
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
!        !write(get_fileId(),*) "Hermitian Conjugate"
!        !Hermitian Conjugate
!        Dk(kSc(1):kEc(1), :, :) = -conjg(Dk(kEc(1):kSc(1):-1, :, :))
!        Dk(:, kSc(2):kEc(2), :) = -conjg(Dk(:, kEc(2):kSc(2):-1, :))
!        Dk(:, :, kSc(3):kEc(3)) = -conjg(Dk(:, :, kEc(3):kSc(3):-1))









!        !write(get_fileId(),*) " "
!        !write(get_fileId(),*) "iFFT Execution"
!        !FFT Execution
!        call dfftw_plan_dft_c2r(plan, RDF%nDim, RDF%kNStep, &
!                               Dk_base, realOut, FFTW_ESTIMATE);
!        call dfftw_execute(plan)
!        call dfftw_destroy_plan(plan)
!
!        !realOut = realOut*sqrt((2*PI)**RDF%nDim/sizesProd)/RDF%kNTotal
!        realOut = realOut/RDF%kNTotal
!
!        !write(*,*) " realOut = ", realOut
!
!        !write(get_fileId(),*) "AFTER iFFT"
!        i = 1
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!        i = 2
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!        i = 3
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!        i = 4
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!        i = 5
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!
!        !write(get_fileId(),*) " "
!        !write(get_fileId(),*) "FFT Execution"
!        call dfftw_plan_dft_r2c(plan, RDF%nDim, RDF%kNStep, &
!                                 realOut, Dk_base, FFTW_ESTIMATE);
!        Dk_base = Dk_base*RDF%kNTotal
!        call dfftw_execute(plan)
!        call dfftw_destroy_plan(plan)
!
!        !write(get_fileId(),*) "AFTER FFT"
!
!        i = 1
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)
!        i = 2
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)
!        i = 3
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)
!        i = 4
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)
!        i = 5
!        !write(get_fileId(),*) " Dk_base 3x ", i
!        !write(get_fileId(),*) Dk_base (i,i,i)
!
!        realOut = realOut*sqrt((2*PI)**RDF%nDim/sizesProd)
!
!        !write(get_fileId(),*) " "
!        !write(get_fileId(),*) "AFTER real mutiplication"
!        i = 1
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!        i = 2
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!        i = 3
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!        i = 4
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!        i = 5
!        !write(get_fileId(),*) " realOut 3x ", i
!        !write(get_fileId(),*) realOut (i,i,i)
!
!
!        !write(get_fileId(),*) "Mapping"
!        do pos = 1, RDF%xNTotal
!            !!write(get_fileId(),*) "pos = ", pos
!            call find_Permutation(pos, RDF%kNStep, posVec)
!            !posVec = posVec -1
!            !write(*,*) "pos = ", pos, "posVec = ", posVec
!            !!write(get_fileId(),*) "allocated(RDF%randField) = ", allocated(RDF%randField)
!            RDF%randField(pos,1) = (realOut(posVec(1), posVec(2), posVec(3)))
!        end do
!
!        !deallocate(Dk)
!        !if(allocated(Dk_base)) deallocate(Dk_base)
!        !if(allocated(realOut)) deallocate(realOut)
!
!        if(associated(RF_2D))    nullify(RF_2D)
!        if(associated(RF_3D))    nullify(RF_3D)
!        if(associated(RFSym_2D)) nullify(RFSym_2D)
!        if(associated(RFSym_3D)) nullify(RFSym_3D)
!        if(associated(SkVec_2D)) nullify(SkVec_2D)
!        if(associated(SkVec_3D)) nullify(SkVec_3D)
!        if(associated(SkSym_2D)) nullify(SkSym_2D)
!        if(associated(SkSym_3D)) nullify(SkSym_3D)
!
!    end subroutine gen_Std_Gauss_FFT_step2_indep

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
