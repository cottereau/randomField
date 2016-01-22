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
    subroutine gen_Std_Gauss_Shinozuka(RDF, MSH, randomK_in)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        logical, intent(in), optional ::randomK_in
        !INPUT
        type(MESH), intent(in) :: MSH

        !LOCAL
        double precision, dimension(:, :), allocatable :: phiK, kVec;
        double precision, dimension(:)   , allocatable :: dgemm_mult;
        double precision, dimension(:,:) , allocatable :: k_x_phi, kSign;
        double precision :: ampMult
        integer :: testIndex = 65
        integer :: n, i, j, m
        logical :: randomK
        integer(kind=8) :: xNTotal

        !!write(get_fileId(),*) "Inside Shinozuka"

        xNTotal = product(MSH%xNStep)
        randomK = .false.
        if(present(randomK_in)) randomK = randomK_in
        call init_random_seed(RDF%seed)

        !!write(get_fileId(),*) "Defining kPoints and SkVec"
        call set_kPoints(RDF, MSH%xStep)
        call set_SkVec(RDF)

        !!write(get_fileId(),*) "Calculating Fields"
        allocate(k_x_phi (xNTotal, 1))
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
                    if(kSign(m, i) == -1) RDF%kPoints(i,:) = - RDF%kPoints(i,:)
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
    subroutine gen_Std_Gauss_Randomization(RDF, MSH)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        !INPUT
        type(MESH), intent(in) :: MSH

        !LOCAL
        logical :: randomK;

        !write(get_fileId(),*) "-----Randomization-----"

        randomK = .true.

        call gen_Std_Gauss_Shinozuka(RDF, MSH, randomK)

    end subroutine gen_Std_Gauss_Randomization


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss_Isotropic(RDF, MSH)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        !INPUT
        type(MESH), intent(in) :: MSH

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

!        !Allocating
!        allocate(rVec (RDF%nDim));
!        allocate(dgemm_mult(RDF%xNTotal))
!
!        !write(get_fileId(),*) "-----Inside Isotropic-----"
!
!        !r Definition
!        call set_kMaxND(RDF%corrMod, rMaxVec)
!        call set_kPoints(RDF)
!        rMax = rMaxVec(1)
!        rDelta  = maxval(RDF%kDelta(:))/5.0D0 !Delta min in between two wave numbers to avoid periodicity
!        rNTotal = ceiling(rMax/rDelta) + 1;
!
!        !Generating random field samples
!        step      = rMax/dble(rNTotal)
!        RDF%randField(:,:) = 0.0D0;
!
!        call init_random_seed(RDF%seed)
!
!        if (RDF%nDim == 2) then
!            allocate(psiN   (rNTotal));
!            allocate(thetaN (rNTotal));
!            allocate(gammaN (rNTotal));
!            do k = 1, RDF%Nmc
!                if(RDF%calculate(k)) then
!                    call random_number(psiN(:))
!                    call random_number(thetaN(:))
!                    call random_number(gammaN(:))
!                    psiN   = 2d0*pi*psiN
!                    thetaN = 2d0*pi*psiN
!                    gammaN = 2d0*pi*gammaN
!
!                    do j = 1, rNTotal
!                        rVec           = [cos(thetaN(j)) * (j-1)*step, &
!                                          sin(thetaN(j)) * (j-1)*step]
!                        Sk             = get_SpectrumND([(j-1)*step], RDF%corrMod); !Obs, here Sk is a scalar
!                        call DGEMM ( "T", "N", RDF%xNTotal, 1, RDF%nDim, &
!                            1.0d0, RDF%xPoints, RDF%nDim, rVec, RDF%nDim, 0.0d0, dgemm_mult, RDF%xNTotal)
!
!                        RDF%randField(:,k) = sqrt(Sk*(j-1)*(dble(step**2))) * gammaN(j) &
!                                            * cos(                           &
!                                            dgemm_mult                &
!                                            + psiN(j)                 &
!                                            )                          &
!                                            + RDF%randField(:,k)
!                    end do
!                else
!                    RDF%randField(:,k) = 0.0
!                end if
!            end do
!
!        else if (RDF%nDim == 3) then
!            !write(*,*) "nDim = 3 !!!"
!            !write(*,*) "k = ",k;
!            allocate(psiN   (rNTotal));
!            allocate(thetaN (rNTotal));
!            allocate(phiN   (rNTotal));
!            allocate(gammaN (rNTotal));
!            do k = 1, RDF%Nmc
!                if(RDF%calculate(k)) then
!                    !write(*,*) "k = ",k;
!                    !write(*,*) "rNTotal = ",rNTotal;
!                    call random_number(phiN(:))
!                    call random_number(thetaN(:))
!                    call random_number(gammaN(:))
!                    call random_number(psiN(:))
!
!                    psiN   = 2*pi*psiN
!                    thetaN = 2*pi*psiN
!                    phiN   = pi*phiN
!                    gammaN = sqrt(12.0)*(gammaN -0.5d0)
!
!                    do j = 1, rNTotal
!                        !write(*,*) "j = ", j
!                        rVec           = [cos(thetaN(j))*sin(phiN(j)) * (j-1)*step, &
!                                          sin(thetaN(j))*sin(phiN(j)) * (j-1)*step, &
!                                          cos(phiN(j))                * (j-1)*step]
!                        Sk             = get_SpectrumND([(j-1)*step], RDF%corrMod);
!                        call DGEMM ( "T", "N", RDF%xNTotal, 1, RDF%nDim, &
!                                    1.0d0, RDF%xPoints, RDF%nDim, rVec, RDF%nDim, 0.0d0, dgemm_mult, RDF%xNTotal)
!                        RDF%randField(:,k) = sqrt(Sk*sin(phiN(j))*step*((j-1)*step)**2) * gammaN(j) &
!                                              * cos(                                             &
!                                              dgemm_mult                                   &
!                                              + psiN(j)                                    &
!                                              )                                            &
!                                              + RDF%randField(:,k)
!                    end do
!                else
!                    RDF%randField(:,k) = 0.0
!                end if
!            end do
!
!        else
!            write(*,*) "ERROR The number of dimensions is not accepted in this method (Isotropic)";
!            write(*,*) "RDF%nDim = ", RDF%nDim;
!            stop
!        end if
!
!        !if(rang == 0) write(*,*) "Spectra (Sk) cut in: ", Sk
!
!        RDF%randField(:,:) = sqrt((1.0d0)/((2.0d0*pi)**(RDF%nDim)))&
!                             * RDF%randField(:,:)
!
!        !RDF%randField = 1.0 ! For Tests
!        !RDF%randField = RDF%rang ! For Tests

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
    subroutine gen_Std_Gauss_FFT(RDF, MSH)

        implicit none
        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        !INPUT
        type(MESH), intent(in) :: MSH
        !LOCAL
        integer(C_INTPTR_T) :: L, M, N
        integer(C_INTPTR_T) :: local_LastDim
        integer(C_INTPTR_T) :: local_LD_offset
        real(C_DOUBLE), pointer :: data_real_2D(:,:), data_real_3D(:,:,:)
        type(C_PTR) :: cdata, plan
        integer(C_INTPTR_T) :: alloc_local
        integer :: sliceSize
        integer :: kNLocal
        integer, dimension(RDF%nDim) :: kNStepLocal
        double precision, dimension(:), allocatable :: gammaK, phiK
        integer :: i, j, k, ind
        double precision :: trashNumber
        integer, dimension(RDF%nDim) :: xNStepGlob
        double precision :: ampMult

        call wLog("gen_Std_Gauss_FFT_init")
        call wLog(" ")


        L = MSH%xNStep(1)
        if(RDF%nDim >= 2) M = MSH%xNStep(2)
        if(RDF%nDim >= 3) N = MSH%xNStep(3)

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

            call wLog("Defining kPoints")
            call set_kPoints(RDF, MSH%xStep)
            call wLog("Defining SkVec")
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

            gammaK       = gammaK -0.5
            !RDF%SkVec(:) = gammak*sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);
            !RDF%SkVec(:) =  sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);

            if(allocated(gammaK)) deallocate(gammaK)
            if(allocated(phik))   deallocate(phik)

            RDF%randField(:,1) = RDF%SkVec

            if(RDF%nDim == 2) then
                plan = fftw_plan_r2r(RDF%nDim, [int(M), int(L)], data_real_2D, data_real_2D, &
                                        [FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
                data_real_2D(:,:) = reshape(RDF%SkVec, [L, M])
                call fftw_execute_r2r(plan, data_real_2D, data_real_2D)
                RDF%randField(:,1) = pack(data_real_2D(1:MSH%xNStep(1), 1:MSH%xNStep(2)), .true.)
                call fftw_destroy_plan(plan)

            else if(RDF%nDim == 3) then
                plan = fftw_plan_r2r(RDF%nDim, [int(N), int(M), int(L)], data_real_3D, data_real_3D, &
                                        [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
                data_real_3D(:,:,:) = reshape(RDF%SkVec, [L, M, N])
                call fftw_execute_r2r(plan, data_real_3D, data_real_3D)
                RDF%randField(:,1) = pack(data_real_3D(1:MSH%xNStep(1), 1:MSH%xNStep(2), 1:MSH%xNStep(3)), .true.)
                call fftw_destroy_plan(plan)

            else
                stop("No FFT method implemented only for 2D and 3D cases")
            end if

            !RDF%randField(:,1) = RDF%SkVec(:) !FOR TESTS
            !RDF%randField(:,1) = RDF%kPoints(:) !FOR TESTS
            !RDF%randField(:,1) = RDF%rang!FOR TESTS

        else
            call wLog("    GLOBAL")
            call fftw_mpi_init()

            xNStepGlob = find_xNStep(MSH%xMinGlob, MSH%xMaxGlob, MSH%xStep)
            call wLog("MSH%xMinGlob = ")
            call wLog(MSH%xMinGlob)
            call wLog("MSH%xMaxGlob = ")
            call wLog(MSH%xMaxGlob)
            call wLog("MSH%xStep = ")
            call wLog(MSH%xStep)
            call wLog("xNStepGlob = ")
            call wLog(xNStepGlob)

            if(RDF%nDim == 2) then
                L = xNStepGlob(1)
                M = xNStepGlob(2)
                alloc_local = fftw_mpi_local_size_2d(M, L, RDF%comm, &
                                                     local_LastDim, local_LD_offset) !FOR MPI
                cdata = fftw_alloc_real(alloc_local)
                call c_f_pointer(cdata, data_real_2D, [L, local_LastDim])

            else if(RDF%nDim == 3) then
                L = xNStepGlob(1)
                M = xNStepGlob(2)
                N = xNStepGlob(3)
                alloc_local = fftw_mpi_local_size_3d(N, M, L, RDF%comm, &
                                                     local_LastDim, local_LD_offset) !FOR MPI
                cdata = fftw_alloc_real(alloc_local)
                call c_f_pointer(cdata, data_real_3D, [L, M, local_LastDim])

            else
                stop("Inside gen_Std_Gauss_FFT dimension not yet implemented for this generation method")

            end if

            !Defining kInit and kEnd
            RDF%kNInit = local_LD_offset + 1
            RDF%kNEnd  = RDF%kNInit + local_LastDim - 1
            call wLog("local_LD_offset")
            call wLog(local_LD_offset)
            call wLog("local_LastDim")
            call wLog(local_LastDim)

            call wLog("RDF%kNInit")
            call wLog(RDF%kNInit)
            call wLog("RDF%kNEnd")
            call wLog(RDF%kNEnd)

            sliceSize = 1
            if(RDF%nDim > 1) sliceSize = product(MSH%xNStep(1:RDF%nDim -1))
            RDF%kNInit = (RDF%kNInit - 1) * sliceSize + 1
            RDF%kNEnd  = RDF%kNEnd *sliceSize

            call wLog("RDF%kNInit")
            call wLog(RDF%kNInit)
            call wLog("RDF%kNEnd")
            call wLog(RDF%kNEnd)

            !RDF%origin  = [1, int(local_LD_offset) + 1]
            !RDF%kExtent = [L , local_M]
            !RDF%origin = 1
            !RDF%origin(RDF%nDim) = int(local_LD_offset) + 1
            RDF%kExtent = MSH%xNStep
            RDF%kExtent(RDF%nDim) = int(local_LastDim)

            call wLog("MSH%origin")
            call wLog(MSH%origin)
            call wLog("MSH%origin (IDEAL)")
            if(MSH%nDim==2) call wLog([1, int(local_LD_offset) + 1])
            if(MSH%nDim==3) call wLog([1, 1, int(local_LD_offset) + 1])
            call wLog("RDF%kExtent")
            call wLog(RDF%kExtent)
            call wLog("RDF%kExtent (IDEAL)")
            if(MSH%nDim==2) call wLog([int(L), int(local_LastDim)])
            if(MSH%nDim==3) call wLog([int(L), int(M), int(local_LastDim)])


            call set_kPoints(RDF, MSH%xStep)
            call set_SkVec(RDF)

            !call wLog("RDF%kPoints")
            !call DispCarvalhol(RDF%kPoints, unit_in = RDF%log_ID)
            !call wLog("RDF%SkVec")
            !call DispCarvalhol(RDF%SkVec, unit_in = RDF%log_ID)

            kNLocal = size(RDF%kPoints,2)
            !kNStepLocal = RDF%kNStep
            !kNStepLocal(RDF%nDim) = kNLocal/sliceSize

            call wLog("MSH%origin")
            call wLog(MSH%origin)
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
!            !RDF%SkVec(:) =  sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);

            if(allocated(gammaK)) deallocate(gammaK)
            if(allocated(phik))   deallocate(phik)

            !cdata = fftw_alloc_real(alloc_local)
            !call c_f_pointer(cdata, data_real_2D, [L,local_M])

            ampMult = 2.0d0*sqrt(product(MSH%xStep)/((2.0d0)**(dble(RDF%nDim))))

            if(RDF%nDim == 2) then
                plan = fftw_mpi_plan_r2r(RDF%nDim, [M, L], data_real_2D, data_real_2D, &
                                         RDF%comm, [FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
                data_real_2D(:,:) = reshape(RDF%SkVec, [L, local_LastDim])
                call fftw_mpi_execute_r2r(plan, data_real_2D, data_real_2D)
                data_real_2D = data_real_2D*sqrt(product(MSH%xStep))
                RDF%randField(:,1) = pack(data_real_2D(1:L,1:local_LastDim), .true.)

            else if(RDF%nDim == 3) then
                plan = fftw_mpi_plan_r2r(RDF%nDim, [N, M, L], data_real_3D, data_real_3D, &
                                         RDF%comm, [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
                data_real_3D(:,:,:) = reshape(RDF%SkVec, [L, M, local_LastDim])
                call fftw_mpi_execute_r2r(plan, data_real_3D, data_real_3D)
                !data_real_3D = data_real_3D*(2.0D0**((RDF%nDim-1))/2.0D0))*sqrt(product(MSH%xStep))
                data_real_3D = data_real_3D*(2.0D0)*sqrt(product(MSH%xStep))
                RDF%randField(:,1) = pack(data_real_3D(1:L, 1:M, 1:local_LastDim), .true.)
            end if

            call wLog("shape(RDF%SkVec)")
            call wLog(shape(RDF%SkVec))
            call wLog("shape(RDF%randField)")
            call wLog(shape(RDF%randField))


            call fftw_destroy_plan(plan)
            call fftw_free(cdata)

            call wLog("L")
            call wLog(L)
            call wLog("M")
            call wLog(M)
            call wLog("local_LastDim")
            call wLog(local_LastDim)
            call wLog("local_LD_offset")
            call wLog(local_LD_offset)
            call wLog("RDF%kNInit")
            call wLog(RDF%kNInit)
            call wLog("RDF%kNEnd")
            call wLog(RDF%kNEnd)
        end if

        if(allocated(gammaK)) deallocate(gammaK)
        if(allocated(phik)) deallocate(phik)

    end subroutine gen_Std_Gauss_FFT

    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    subroutine allocate_randField(RDF, xNstep, randField)
        implicit none
        !INPUT AND OUTPUT
        type(RF)   :: RDF
        double precision, dimension(:,:), allocatable, target :: randField
        !INPUT
        integer, dimension(:), intent(in) :: xNstep
        !LOCAL
        integer(kind=8) :: xNTotal

        xNTotal = product(xNStep)

        if(allocated(randField)) then
            if(.not.(size(randField,1) == xNTotal .and. size(randField,1) == RDF%Nmc)) then
                nullify(RDF%randField)
                deallocate(randField)
            end if
        end if

        if(.not.allocated(randField)) then
            allocate(randField(xNTotal, RDF%Nmc))
            RDF%randField => randField
            call wLog(" Inside allocate_randField")
            call wLog("     IN xNStep = ")
            call wLog(xNStep)
            if(RDF%nDim == 2) RDF%RF_2D(1:xNStep(1),1:xNStep(2)) => randField
            if(RDF%nDim == 3) RDF%RF_3D(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => randField
        end if

    end subroutine allocate_randField

    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    subroutine normalize_randField(RDF, xNTotal, randField)
        implicit none
        !INPUT OUTPUT
        double precision, dimension(:,:), intent(inout) :: randField
        !INPUT
        integer(kind=8), intent(in) :: xNTotal
        type(RF), intent(in) :: RDF
        !LOCAL
        double precision, dimension(RDF%Nmc) :: sumRF, sumRFsquare
        double precision, dimension(RDF%Nmc) :: totalSumRF, totalSumRFsquare;
        double precision, dimension(RDF%Nmc) :: evntAvg, evntStdDev;
        !integer, dimension(:), allocatable :: xNTotal_Vec, deplacement
        integer :: code
        integer :: i


        call wLog("Calculating Average and stdVar")

        !Calculating
        call wLog("sumRF(1) = ")
        call wLog(sumRF(1))
        call wLog("sumRFsquare(1) = ")
        call wLog(sumRFsquare(1))

        !Total Number of Points
        call wLog("xNTotal = ")
        call wLog(xNTotal)

        !Average Correction
        sumRF(:)       = sum( RDF%randField    , dim = 1)
        call MPI_ALLREDUCE (sumRF,totalSumRF,RDF%Nmc,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,RDF%comm,code)
        evntAvg      = totalSumRF/xNTotal;
        call wLog("Initial Average = ")
        call wLog(evntAvg)

        do i = 1, RDF%Nmc
            randField(:,i) = randField(:,i) - evntAvg(i)
        end do

        !Verifying Average
        sumRF(:)       = sum( RDF%randField    , dim = 1)
        call MPI_ALLREDUCE (sumRF,totalSumRF,RDF%Nmc,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,RDF%comm,code)
        evntAvg      = totalSumRF/xNTotal;
        call wLog("Final Average = ")
        call wLog(evntAvg)


        !Standard Deviation Correction
        sumRFsquare(:) = sum((RDF%randField)**2, dim = 1)
        call MPI_ALLREDUCE (sumRFsquare,totalSumRFsquare,RDF%Nmc,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,RDF%comm,code)
        evntStdDev   = sqrt(totalSumRFsquare/dble(xNTotal)) !Mean of random field is supposed to be 0
        call wLog("Initial StdDev = ")
        call wLog(evntStdDev)
        do i = 1, RDF%Nmc
            randField(:,i) = randField(:,i)/evntStdDev(i)
        end do

        !Verifying Standard Deviation
        sumRFsquare(:) = sum((RDF%randField)**2, dim = 1)
        call MPI_ALLREDUCE (sumRFsquare,totalSumRFsquare,RDF%Nmc,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,RDF%comm,code)
        evntStdDev   = sqrt(totalSumRFsquare/dble(xNTotal)) !Mean of random field is supposed to be 0
        call wLog("Final StdDev = ")
        call wLog(evntStdDev)

    end subroutine normalize_randField


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
