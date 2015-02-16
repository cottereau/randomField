module test_func_RF

    use randomFieldND
    use statistics_RF
    use writeResultFile_RF
    use readFile_RF
    use displayCarvalhol
    use mesh_RF
    use constants_RF
    use mpi
    use write_Log_File
    use common_variables_RF

    implicit none

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine func_test_001_speed_xStep(xMin, xMax,                    &
                                        corrL, corrMod,                 &
                                        margiFirst, fieldAvg, fieldVar, &
                                        Nmc, method,                    &
                                        step_initial, step_nIter,       &
                                        step_mult, step_add,            &
                                        xMinGlob, xMaxGlob)

        implicit none

        !INPUTS
        double precision, dimension(:), intent(inout) :: xMin, xMax;
        double precision, dimension(:), intent(in) :: corrL;
        character (len=*)             , intent(in) :: corrMod, margiFirst;
        double precision              , intent(in) :: fieldAvg, fieldVar;
        integer                       , intent(in) :: Nmc, method;
        integer                       , intent(in) :: step_nIter
        double precision, dimension(:), intent(in) :: step_mult, step_add, step_initial
        double precision, dimension(:), intent(in) :: xMinGlob, xMaxGlob;
        integer, parameter :: seedStart = 3

        !LOCAL
        integer :: nDim;
        integer :: i;
        integer :: comm, rang, code;
        double precision, dimension(:, :), allocatable :: randField, xPoints;
        double precision :: t1, t2, t3, t4;
        double precision  , dimension(:) , allocatable :: avg_Gauss, stdDev_Gauss, &
                                                          avg_Trans, stdDev_Trans
        double precision  , dimension(:,:) , allocatable :: avg_Gauss_evnt, stdDev_Gauss_evnt, &
                                                            avg_Trans_evnt, stdDev_Trans_evnt
        double precision, dimension(:), allocatable :: xStep;
        double precision :: multInd
        integer, dimension(:), allocatable :: seed
        double precision :: testRand
        character (len=200) :: step_result_path

        step_result_path = string_vec_join([results_path,"/",results_folder_name,"/",xStep_folder_name])


        nDim = size(xMax)
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK (comm ,rang,code)

        allocate(xStep(nDim))
        allocate(avg_Gauss(step_nIter), stdDev_Gauss(step_nIter))
        allocate(avg_Trans(step_nIter), stdDev_Trans(step_nIter))
        allocate(avg_Gauss_evnt(step_nIter, Nmc), stdDev_Gauss_evnt(step_nIter, Nmc))
        allocate(avg_Trans_evnt(step_nIter, Nmc), stdDev_Trans_evnt(step_nIter, Nmc))

        if(rang == 0) write(*,*) ""
        if(rang == 0) write(*,*) "----------------------------------------------------------------------------";
        if(rang == 0) write(*,*) "----------------------------------------------------------------------------";
        if(rang == 0) write(*,*) "RUNNING SPEED UP TEST VARYING xStep-----------------------------------------";
        if(rang == 0) write(*,*) "----------------------------------------------------------------------------";
        if(rang == 0) write(*,*) "----------------------------------------------------------------------------";

        do i = 1, step_nIter

            if(rang == 0) write(*,*) ""
            if(rang == 0) write(*,*) "ITERATION xStep",       i;

            write(get_fileId(),*) ""
            write(get_fileId(),*) "------------------------------------------------------------------------";
            write(get_fileId(),*) "ITERATION xStep",       i, "--------------------------------------------";
            write(get_fileId(),*) "------------------------------------------------------------------------";

            !Putting everything to the initial condition
            if (allocated(xPoints))   deallocate(xPoints)
            if (allocated(randField)) deallocate(randField)
            call calculate_random_seed(seed, seedStart)

            xStep = step_initial         &
                    *(step_mult)**(i-1) &
                    + step_add*(i-1)

            write(get_fileId(),*) ""
            write(get_fileId(),*) "  xStep = ", xStep;
            write(get_fileId(),*) ""

            call cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
            call set_XPoints(xMin, xMax, xStep, xPoints)
            call restoreBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
            allocate(randField(size(xPoints,2), Nmc))

            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) RAW")
            call dispCarvalhol(transpose(xPoints(:, :)), "transpose(xPoints(:, :)) BEG")
            !call dispCarvalhol(transpose(xPoints(:, size(xPoints,2)-20:)), "transpose(xPoints(:, size(xPoints,2)-20:)) END")

            call write_generation_spec(xMinGlob, xMaxGlob, xStep,              &
                                       corrL, corrMod,                 &
                                       margiFirst, fieldAvg, fieldVar, &
                                       Nmc, method, seed,              &
                                       rang, step_result_path, stringNumb_join("iteration_",i))

            write(get_fileId(),*) "-> Creating Standard Gaussian Random field (unstructured)";
            t1 = MPI_Wtime();
            randField = 0.0
            call create_Std_Gaussian_Field_Unstruct (xPoints, corrL, corrMod, Nmc,  &
                                                      randField, method, seed)
            t2 = MPI_Wtime();

            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            call write_Mono_XMF_h5(xPoints, randField, "gauss_", rang, step_result_path, &
                                            MPI_COMM_WORLD, ["it_", "_proc_"], [i, rang], i)


            write(get_fileId(),*) "-> Calculating Statistics Before Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                  &
                                                  avg_Gauss(i), stdDev_Gauss(i),              &
                                                  comm,                                       &
                                                  avg_Gauss_evnt(i,:), stdDev_Gauss_evnt(i,:))


            write(get_fileId(),*) "-> Transforming Random field (unstructured)";
            t3 = MPI_Wtime();
            call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)
            t4 = MPI_Wtime();

            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            call write_Mono_XMF_h5(xPoints, randField, "trans_", rang, step_result_path, &
                                            MPI_COMM_WORLD, ["it_", "_proc_"], [i, rang], i)

            write(get_fileId(),*) "-> Calculating Statistics After Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                   &
                                                  avg_Trans(i), stdDev_Trans(i),               &
                                                  comm,                                        &
                                                  avg_Trans_evnt(i,:), stdDev_Trans_evnt(i,:))
            call dispCarvalhol(randField(:, :), "randField(:, :) AFTER ALL")
        end do

        !Writing
        if(rang == 0) then
            write(*,*) "-----CONVERGENCE TEST CHANGING xStep (Speed Up)-----"
            write(*,*) "     By Event "
            call DispCarvalhol(avg_Gauss_evnt   , "avg_Gauss_evnt")
            call DispCarvalhol(stdDev_Gauss_evnt, "stdDev_Gauss_evnt")
            write(*,*) ""
            write(*,*) "     GLOBAL "
            write(*,*) ""
            write(*,'(A4, A10, A2, 4A20)') "Iter", "multInd", " )", "Avg-Gauss", "StdDev-Gauss", &
                "Avg-Trans", "StdDev-Trans"
            do i = 1, step_nIter
                multInd = product((step_mult)**(i-1) + step_add*(i-1))
                write(*,'(I4, F10.5, A2, 4F20.8)') i, multInd , " )",              &
                    avg_Gauss(i), stdDev_Gauss(i), &
                    avg_Trans(i), stdDev_Trans(i)
            end do
            write(*,'(A16, 4F20.8)') "Reference      )", &
                0.0d0, 1.0d0,  &
                fieldAvg, fieldVar
        end if

        if (allocated(xStep))     deallocate(xStep)
        if (allocated(xPoints))   deallocate(xPoints)
        if (allocated(randField)) deallocate(randField)
        if (allocated(seed))      deallocate(seed)
        deallocate(avg_Gauss, stdDev_Gauss)
        deallocate(avg_Trans, stdDev_Trans)
        deallocate(avg_Gauss_evnt, stdDev_Gauss_evnt)
        deallocate(avg_Trans_evnt, stdDev_Trans_evnt)

    end subroutine func_test_001_speed_xStep

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine func_test_002_speed_Nmc(xMin, xMax,                      &
                                      corrL, corrMod,                  &
                                      margiFirst, fieldAvg, fieldVar,  &
                                      xStep, method,                   &
                                      nmc_initial, nmc_nIter,          &
                                      nmc_mult, nmc_add,            &
                                      xMinGlob, xMaxGlob)

        implicit none

        !INPUTS
        double precision, dimension(:), intent(in) :: xMin, xMax;
        double precision, dimension(:), intent(in) :: corrL;
        character (len=*)             , intent(in) :: corrMod, margiFirst;
        double precision              , intent(in) :: fieldAvg, fieldVar;
        double precision, dimension(:), intent(in) :: xStep;
        integer                       , intent(in) :: method;
        integer                       , intent(in) :: nmc_initial, nmc_nIter
        double precision              , intent(in) :: nmc_mult, nmc_add
        double precision, dimension(:), intent(in), optional :: xMinGlob, xMaxGlob
        integer :: seedStart = 3

        !LOCAL
        integer :: nDim, Nmc;
        integer :: i, j, fileId, nColumns;
        integer :: code, rang, error, nb_procs;
        character (len=40) :: doubleFmt;
        double precision, dimension(:, :), allocatable :: randField, xPoints;
        double precision :: t1, t2, t3, t4;
        double precision  , dimension(:) , allocatable :: avg_Gauss, stdDev_Gauss, &
                                                          avg_Trans, stdDev_Trans
        double precision :: multInd
        integer :: comm
        integer, dimension(:), allocatable :: seed
        double precision :: testRand
        integer :: Nmc_iter
        character (len=200) :: Nmc_result_path

        Nmc_result_path = string_vec_join([results_path,"/",results_folder_name,"/",Nmc_folder_name])

        nDim = size(xMax)
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK (comm ,rang,code)

        allocate(avg_Gauss(nmc_nIter), stdDev_Gauss(nmc_nIter))
        allocate(avg_Trans(nmc_nIter), stdDev_Trans(nmc_nIter))

        write(get_fileId(),*) "-> Creating xPoints";
        call set_XPoints(xMin, xMax, xStep, xPoints)

        do i = 1, nmc_nIter

            if(rang == 0) write(*,*) ""
            if(rang == 0) write(*,*) "-------------------------------------------------------------------------------";
            if(rang == 0) write(*,*) "ITERATION Nmc Speed",       i, "-----------------------------------------------";
            if(rang == 0) write(*,*) "-------------------------------------------------------------------------------";

            !Putting everything to the initial condition
            if (allocated(randField)) deallocate(randField)
            call calculate_random_seed(seed, seedStart)

            Nmc = nmc_initial         &
                    *(nmc_mult)**(i-1) &
                    + nmc_add*(i-1)

            write(get_fileId(),*) ""
            write(get_fileId(),*) "  Nmc = ", Nmc;
            write(get_fileId(),*) ""

            allocate(randField(size(xPoints,2), Nmc))

            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) RAW")
            !call dispCarvalhol(transpose(xPoints(:, 1:20)), "transpose(xPoints(:, 1:20)) BEG")
            !call dispCarvalhol(transpose(xPoints(:, size(xPoints,2)-20:)), "transpose(xPoints(:, size(xPoints,2)-20:)) END")

            call write_generation_spec(xMinGlob, xMaxGlob, xStep,              &
                                       corrL, corrMod,                 &
                                       margiFirst, fieldAvg, fieldVar, &
                                       Nmc, method, seed,              &
                                       rang, Nmc_result_path, stringNumb_join("iteration_",i))

            t1 = MPI_Wtime();

            write(get_fileId(),*) "-> Creating Standard Gaussian Random field (unstructured)";
            randField = 0.0
            call create_Std_Gaussian_Field_Unstruct (xPoints, corrL, corrMod, Nmc,  &
                                                      randField, method, seed)
            t2 = MPI_Wtime();

            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            call write_Mono_XMF_h5(xPoints, randField, "gauss_", rang, Nmc_result_path, &
                                            MPI_COMM_WORLD, ["it_", "_proc_"], [i, rang], i)

            !call write_MatlabTable(randField, "visu/"//trim(stringNumb_join("Nmc_Iteration_",i)))
            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) GAUSSIAN")


            write(get_fileId(),*) "-> Calculating Statistics Before Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                  &
                                                  avg_Gauss(i), stdDev_Gauss(i),              &
                                                  comm)

            t3 = MPI_Wtime();

            write(get_fileId(),*) "-> Transforming Random field (unstructured)";
            call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

            t4 = MPI_Wtime();

            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            call write_Mono_XMF_h5(xPoints, randField, "trans_", rang, nmc_result_path, &
                                            MPI_COMM_WORLD, ["it_", "_proc_"], [i, rang], i)
            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) TRANSFORMED")

            write(get_fileId(),*) "-> Calculating Statistics After Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                   &
                                                  avg_Trans(i), stdDev_Trans(i),               &
                                                  comm)

        end do

       if(rang == 0) then
           write(*,*) "-----CONVERGENCE TEST CHANGING Nmc (Speed)-----"
           write(*,*) ""
           write(*,*) "     GLOBAL "
           write(*,*) ""
           write(*,'(A4, A10, A2, 4A20)') "Iter", "multInd", " )", "Avg-Gauss", "StdDev-Gauss", &
                                                              "Avg-Trans", "StdDev-Trans"
           do i = 1, nmc_nIter
               Nmc_iter = nmc_initial *(nmc_mult)**(i-1) + nmc_add*(i-1)
               write(*,'(I4, I10, A2, 4F20.8)') i, Nmc_iter , " )",             &
                                                  avg_Gauss(i), stdDev_Gauss(i), &
                                                  avg_Trans(i), stdDev_Trans(i)
           end do
           write(*,'(A16, 4F20.8)') "Reference      )", &
                                    0.0d0, 1.0d0,      &
                                    fieldAvg, fieldVar
       end if

        if (allocated(xPoints))   deallocate(xPoints)
        if (allocated(randField)) deallocate(randField)
        if (allocated(seed))      deallocate(seed)
        deallocate(avg_Gauss, stdDev_Gauss)
        deallocate(avg_Trans, stdDev_Trans)

    end subroutine func_test_002_speed_Nmc

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine func_test_003_conv_xStep(xMin, xMax,                    &
                                        corrL, corrMod,                 &
                                        margiFirst, fieldAvg, fieldVar, &
                                        Nmc, method,                    &
                                        step_initial, step_nIter,       &
                                        step_div,                       &
                                        xMinGlob, xMaxGlob)

        implicit none

        !INPUTS
        double precision, dimension(:), intent(inout) :: xMin, xMax;
        double precision, dimension(:), intent(in) :: corrL;
        character (len=*)             , intent(in) :: corrMod, margiFirst;
        double precision              , intent(in) :: fieldAvg, fieldVar;
        integer                       , intent(in) :: Nmc, method;
        integer                       , intent(in) :: step_nIter
        double precision, dimension(:), intent(in) :: step_initial
        integer         , dimension(:), intent(in) :: step_div
        double precision, dimension(:), intent(in) :: xMinGlob, xMaxGlob;
        integer, parameter :: seedStart = 3

        !LOCAL
        integer :: nDim;
        integer :: i, j, k;
        integer :: comm, rang, code;
        double precision, dimension(:, :), allocatable :: randField, xPoints;
        double precision :: t1, t2, t3, t4;
        double precision  , dimension(:) , allocatable :: avg_Gauss, stdDev_Gauss, &
                                                          avg_Trans, stdDev_Trans
        double precision  , dimension(:,:) , allocatable :: avg_Gauss_evnt, stdDev_Gauss_evnt, &
                                                            avg_Trans_evnt, stdDev_Trans_evnt
        double precision, dimension(:), allocatable :: xStep, xStep_Iter;
        double precision :: multInd
        integer, dimension(:), allocatable :: seed
        double precision :: testRand
        integer, dimension(:), allocatable :: step_integer, xNStep, xNStep_Iter
        logical, dimension(:,:,:), allocatable :: xPoints_mask, randomField_mask
        integer :: beg, step, end, prod_iter
        integer :: nPlanes, iter, cont
        character (len=200) :: step_result_path

        step_result_path = string_vec_join([results_path,"/",results_folder_name,"/",xStep_folder_name])

        nDim = size(xMax)
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK (comm ,rang,code)

        allocate(xStep(nDim))
        allocate(xStep_Iter(nDim))
        allocate(xNStep(nDim))
        allocate(xNStep_Iter(nDim))
        allocate(step_integer(nDim))
        allocate(avg_Gauss(step_nIter), stdDev_Gauss(step_nIter))
        allocate(avg_Trans(step_nIter), stdDev_Trans(step_nIter))
        allocate(avg_Gauss_evnt(step_nIter, Nmc), stdDev_Gauss_evnt(step_nIter, Nmc))
        allocate(avg_Trans_evnt(step_nIter, Nmc), stdDev_Trans_evnt(step_nIter, Nmc))

        !Test
        do i = 1, size(step_div)
            if(step_div(i) < 0) then
                write (*,*) "step_div = ", step_div
                stop "ERROR, In func_test_003_conv_xStep  step_div should be greater than 0"
            end if
        end do
        !write(*,*) "step_integer = ", step_integer

        xStep = step_initial                 &
                *(1.0d0/dble(step_div))**(step_nIter-1)
        !write(*,*) "xStep = ", xStep

        !Putting everything to the initial condition
        if (allocated(xPoints))   deallocate(xPoints)
        if (allocated(randField)) deallocate(randField)
        call calculate_random_seed(seed, seedStart)

!        if(rang == 0) write(*,*) ""
!        if(rang == 0) write(*,*) "  xStep    = ", xStep;
!        if(rang == 0) write(*,*) "  step_div = ", step_div;
!        if(rang == 0) write(*,*) ""

!        !TEST
!        call set_XPoints(xMin, xMax, step_initial, xPoints)
!        write(*,*) "INIT shape(xPoints) = ", shape(xPoints)
!        if (allocated(xPoints))   deallocate(xPoints)
!        !TEST

        write(get_fileId(),*) "-> Creating xPoints";
        call cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
        call set_XPoints(xMin, xMax, xStep, xPoints)
        call restoreBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
!        call dispCarvalhol(transpose(xPoints), "transpose(xPoints)")
        allocate(xPoints_mask(nDim, size(xPoints,2), step_nIter))
        allocate(randomField_mask(size(xPoints,2), Nmc, step_nIter))
        xNStep = find_xNStep(xMin, xMax, xStep);
        !write(*,*) "END shape(xPoints_mask) = ", shape(xPoints_mask)
        !write(*,*) "END shape(randomField_mask) = ", shape(randomField_mask)

        !Building xPoints_mask
        xPoints_mask(:,:,:) = .true.
        randomField_mask(:,:,:) = .true.

        do iter = 1, step_nIter
            xStep_Iter   = step_initial/(dble(step_div)**(dble(iter)-1))
            !xNStep_Iter  = find_xNStep(xMin, xMax, xStep_Iter);
            !step_integer = nint(xStep_Iter/xStep)
            !prod_iter    = product(step_integer)
            !write(*,*) "step_integer = ", step_integer
            !write(*,*) "prod_iter    = ", prod_iter

            do i = 1, size(xPoints,2)
                !write(*,*) "i=",i
                if(randomField_mask(i, 1, iter)) then
                    do k = 1, nDim
                        !write(*,*) "k=",i
                        if( .not.(                                                         &
                            (abs(                                                          &
                            (xPoints(k, i)-xPoints(k, 1))/dble(xStep_Iter(k))              &
                            -dble(nint((xPoints(k, i)-xPoints(k, 1))/dble(xStep_Iter(k)))) &
                            ))                                                             &
                            < TOLERANCE                                                    &
                            )) then

                            !write(*,*) "(xPoints(k, i)-xPoints(k, 0))/dble(xStep_Iter(k))"
                            !write(*,*) " = ", (xPoints(k, i)-xPoints(k, 0))/dble(xStep_Iter(k))
                            !write(*,*) "dble(nint((xPoints(k, i)-xPoints(k, 0))/dble(xStep_Iter(k))))"
                            !write(*,*) " = ", dble(nint((xPoints(k, i)-xPoints(k, 0))/dble(xStep_Iter(k))))
                            !write(*,*) "Point ", i, "is not in this iteraction"
                            !write(*,*) ""
                            xPoints_mask(:, i, iter)     = .false.
                            randomField_mask(i, :, iter) = .false.
                        end if
                    end do
                end if
            end do
        end do

        !call dispCarvalhol(transpose(xPoints_mask(:,:,1)), "transpose(xPoints_mask(:,:,1))")
        !call dispCarvalhol(randomField_mask(:,:,1), "randomField_mask(:,:,1)")

        write(get_fileId(),*) "-> Creating Standard Gaussian Random field (unstructured)";
        call calculate_random_seed(seed, seedStart)
        allocate(randField(size(xPoints,2), Nmc))
        randField = 0.0
        call create_Std_Gaussian_Field_Unstruct (xPoints, corrL, corrMod, Nmc,  &
                                                  randField, method, seed)

!        call dispCarvalhol(randField, "randField")

        !Calculating statistics of the Gaussian field
        if(rang == 0) write(*,*) "-> Calculating Statistics Before Transformation";
        do i = 1, step_nIter

            xStep_Iter   = step_initial/(dble(step_div)**(dble(i)-1))
            xNStep_Iter  = find_xNStep(xMin, xMax, xStep_Iter);
            step_integer = nint(xStep_Iter/xStep)

            if(rang == 0) write(*,*) "------------------"
            if(rang == 0) write(*,*) "  ITERATION ", i;
            if(rang == 0) write(*,*) "  xStep        = ", xStep_Iter;
            if(rang == 0) write(*,*) "  step_integer = ", step_integer;
            if(rang == 0) write(*,*) "------------------"

            call write_generation_spec(xMinGlob, xMaxGlob, xStep_Iter,         &
                                       corrL, corrMod,                 &
                                       margiFirst, fieldAvg, fieldVar, &
                                       Nmc, method, seed,         &
                                       rang, step_result_path, stringNumb_join("iteration_",i))

!            call dispCarvalhol(                                                      &
!                               transpose(&
!                               reshape(                                              &
!                                       pack(xPoints, mask = xPoints_mask(:,:,i)), &
!                                       shape = [nDim, count(xPoints_mask(1,:,i))] &
!                                      )                                              &
!                               ),&
!                               "Points")

!            call dispCarvalhol(                                                      &
!                               reshape(                                              &
!                                       pack(randField, mask = randomField_mask(:,:,iter)),  &
!                                       shape = [count(randomField_mask(:,1,iter)), Nmc] &
!                                      ),                                              &
!                               "Random Field")

            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            call write_Mono_XMF_h5(reshape(                                           &
                                           pack(xPoints, mask = xPoints_mask(:,:,i)), &
                                           shape = [nDim, count(xPoints_mask(1,:,i))] &
                                           ),                                         &
                                   reshape(                                                    &
                                           pack(randField, mask = randomField_mask(:,:,i)), &
                                           shape = [count(randomField_mask(:,1,i)), Nmc]    &
                                           ),                                                  &
                                   "gauss_", rang, step_result_path, &
                                   MPI_COMM_WORLD, ["it_", "_proc_"], [i, rang], i)

            call calculate_average_and_stdVar_MPI(&
                                        reshape(                                                  &
                                                pack(randField, mask = randomField_mask(:,:,i)),  &
                                                shape = [count(randomField_mask(:,1,i)), Nmc]     &
                                                ),                                                &
                                                avg_Gauss(i), stdDev_Gauss(i),                    &
                                                comm,                                             &
                                                avg_Gauss_evnt(i,:), stdDev_Gauss_evnt(i,:))
        end do

        write(get_fileId(),*) "-> Transforming Random field (unstructured)";
        call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

        !Calculating statistics of the Transformed field
        write(get_fileId(),*) "-> Calculating Statistics After Transformation";
        do i = 1, step_nIter

            xStep_Iter   = step_initial/(dble(step_div)**(dble(i)-1))
            xNStep_Iter  = find_xNStep(xMin, xMax, xStep_Iter);
            step_integer = nint(xStep_Iter/xStep)

            if(rang == 0) write(*,*) "------------------"
            if(rang == 0) write(*,*) "  ITERATION ", i;
            if(rang == 0) write(*,*) "  xStep        = ", xStep_Iter;
            if(rang == 0) write(*,*) "  step_integer = ", step_integer;
            if(rang == 0) write(*,*) "------------------"

!            call dispCarvalhol(                                                      &
!                               reshape(                                              &
!                                       pack(randField, mask = randomField_mask(:,:,iter)),  &
!                                       shape = [count(randomField_mask(:,1,iter)), Nmc] &
!                                      ),                                              &
!                               "Random Field for iteraction")

            call calculate_average_and_stdVar_MPI(&
                                        reshape(                                                     &
                                                pack(randField, mask = randomField_mask(:,:,i)),  &
                                                shape = [count(randomField_mask(:,1,i)), Nmc]     &
                                                ),                                                    &
                                                avg_Trans(i), stdDev_Trans(i),               &
                                                comm,                                        &
                                                avg_Trans_evnt(i,:), stdDev_Trans_evnt(i,:))
        end do

       if(rang == 0) then
           write(*,*) "-----CONVERGENCE TEST CHANGING xStep (no SpeedUp)-----"
           write(*,*) "     By Event "
           call DispCarvalhol(avg_Gauss_evnt   , "avg_Gauss_evnt")
           call DispCarvalhol(stdDev_Gauss_evnt, "stdDev_Gauss_evnt")
           write(*,*) ""
           write(*,*) "     GLOBAL "
           write(*,*) ""
           write(*,'(A4, A10, A2, 4A20)') "Iter", "xNStep", " )", "Avg-Gauss", "StdDev-Gauss", &
                                                              "Avg-Trans", "StdDev-Trans"
           do i = 1, step_nIter
                xStep_Iter   = step_initial/(dble(step_div)**(dble(i)-1))
                xNStep_Iter  = find_xNStep(xMin, xMax, xStep_Iter);
               write(*,'(I4, I10, A2, 4F20.8)') i, product(xNStep_Iter) , " )",              &
                                                  avg_Gauss(i), stdDev_Gauss(i), &
                                                  avg_Trans(i), stdDev_Trans(i)
           end do
           write(*,'(A16, 4F20.8)') "Reference      )", &
                                    0.0d0, 1.0d0,  &
                                    fieldAvg, fieldVar
       end if

        if (allocated(xStep))          deallocate(xStep)
        if (allocated(xStep_Iter))     deallocate(xStep_Iter)
        if (allocated(xNStep_Iter))    deallocate(xNStep_Iter)
        if (allocated(xNStep))         deallocate(xNStep)
        if (allocated(step_integer))   deallocate(step_integer)
        if (allocated(xPoints))        deallocate(xPoints)
        if (allocated(randField))      deallocate(randField)
        if (allocated(seed))           deallocate(seed)
        deallocate(xPoints_mask, randomField_mask)
        deallocate(avg_Gauss, stdDev_Gauss)
        deallocate(avg_Trans, stdDev_Trans)
        deallocate(avg_Gauss_evnt, stdDev_Gauss_evnt)
        deallocate(avg_Trans_evnt, stdDev_Trans_evnt)

    end subroutine func_test_003_conv_xStep
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine func_test_004_conv_Nmc(xMin, xMax,                     &
                                       corrL, corrMod,                  &
                                       margiFirst, fieldAvg, fieldVar,  &
                                       xStep, method,                   &
                                       nmc_initial, nmc_nIter,          &
                                       nmc_mult, nmc_add,             &
                                       xMinGlob, xMaxGlob)

        implicit none

        !INPUTS
        double precision, dimension(:), intent(in) :: xMin, xMax;
        double precision, dimension(:), intent(in) :: corrL;
        character (len=*)             , intent(in) :: corrMod, margiFirst;
        double precision              , intent(in) :: fieldAvg, fieldVar;
        double precision, dimension(:), intent(in) :: xStep;
        integer                       , intent(in) :: method;
        integer                       , intent(in) :: nmc_initial, nmc_nIter
        double precision              , intent(in) :: nmc_mult, nmc_add
        double precision, dimension(:), intent(in), optional :: xMinGlob, xMaxGlob
        integer :: seedStart = 5

        !LOCAL
        integer :: nDim, Nmc, Nmc_iter;
        integer :: i, j, fileId, nColumns;
        integer :: code, rang, error, nb_procs;
        character (len=40) :: doubleFmt;
        double precision, dimension(:, :), allocatable :: randField, xPoints;
        double precision  , dimension(:) , allocatable :: avg_Gauss, stdDev_Gauss, &
                                                          avg_Trans, stdDev_Trans
        double precision :: multInd
        integer :: comm
        integer, dimension(:), allocatable :: seed
        double precision :: testRand
        character (len=200) :: Nmc_result_path

        Nmc_result_path = string_vec_join([results_path,"/",results_folder_name,"/",Nmc_folder_name])

        nDim = size(xMax)
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK (comm ,rang,code)

        allocate(avg_Gauss(nmc_nIter), stdDev_Gauss(nmc_nIter))
        allocate(avg_Trans(nmc_nIter), stdDev_Trans(nmc_nIter))

        if(rang == 0) write(*,*) "-> Creating xPoints";
        call set_XPoints(xMin, xMax, xStep, xPoints)

        Nmc = nmc_initial         &
              *(nmc_mult)**(nmc_nIter-1) &
              + nmc_add*(nmc_nIter-1)

        !Putting everything to the initial condition
        if (allocated(randField)) deallocate(randField)
        call calculate_random_seed(seed, seedStart)
        allocate(randField(size(xPoints,2), Nmc))

        if(rang == 0) write(*,*) "-> Creating Standard Gaussian Random field (unstructured)";
        randField = 0.0
        call create_Std_Gaussian_Field_Unstruct (xPoints, corrL, corrMod, Nmc,  &
                                                  randField, method, seed)

        !Calculating statistics of the Gaussian field
        if(rang == 0) write(*,*) "-> Calculating Statistics Before Transformation";
        do i = 1, nmc_nIter

            Nmc_iter = nmc_initial *(nmc_mult)**(i-1) + nmc_add*(i-1)

            if(rang == 0) write(*,*) "------------------"
            if(rang == 0) write(*,*) "  ITERATION ", i;
            if(rang == 0) write(*,*) "  Nmc = ", Nmc_iter;
            if(rang == 0) write(*,*) "------------------"

            call write_generation_spec(xMinGlob, xMaxGlob, xStep,              &
                                       corrL, corrMod,                 &
                                       margiFirst, fieldAvg, fieldVar, &
                                       Nmc_iter, method, seed,         &
                                       rang, Nmc_result_path, stringNumb_join("iteration_",i))

            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            call write_Mono_XMF_h5(xPoints, randField(:,:Nmc_iter), "gauss_", rang, nmc_result_path, &
                                            MPI_COMM_WORLD, ["it_", "_proc_"], [i, rang], i)

            call calculate_average_and_stdVar_MPI(randField(:,:Nmc_iter),                     &
                                                  avg_Gauss(i), stdDev_Gauss(i),              &
                                                  comm)
        end do

        !Transforming Random field
        if(rang == 0) write(*,*) "-> Transforming Random field (unstructured)";
        call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

        !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) TRANSFORMED")

        if(rang == 0) write(*,*) "-> Calculating Statistics After Transformation";
        do i = 1, nmc_nIter

            Nmc_iter = nmc_initial *(nmc_mult)**(i-1) + nmc_add*(i-1)

            if(rang == 0) write(*,*) "------------------"
            if(rang == 0) write(*,*) "  ITERATION ", i;
            if(rang == 0) write(*,*) "  Nmc = ", Nmc_iter;
            if(rang == 0) write(*,*) "------------------"

            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            call write_Mono_XMF_h5(xPoints, randField(:,:Nmc_iter), "trans_", rang, nmc_result_path, &
                                            MPI_COMM_WORLD, ["it_", "_proc_"], [i, rang], i)

            call calculate_average_and_stdVar_MPI(randField(:,:Nmc_iter),                                   &
                                                  avg_Trans(i), stdDev_Trans(i),               &
                                                  comm)
        end do

       if(rang == 0) then
           write(*,*) "-----CONVERGENCE TEST CHANGING Nmc-----"
           write(*,*) ""
           write(*,*) "     GLOBAL "
           write(*,*) ""
           write(*,'(A4, A10, A2, 4A20)') "Iter", "Nmc", " )", "Avg-Gauss", "StdDev-Gauss", &
                                                              "Avg-Trans", "StdDev-Trans"
           do i = 1, nmc_nIter
               Nmc_iter = nmc_initial *(nmc_mult)**(i-1) + nmc_add*(i-1)
               write(*,'(I4, I10, A2, 4F20.8)') i, Nmc_iter , " )",             &
                                                  avg_Gauss(i), stdDev_Gauss(i), &
                                                  avg_Trans(i), stdDev_Trans(i)
           end do
           write(*,'(A16, 4F20.8)') "Reference      )", &
                                    0.0d0, 1.0d0,      &
                                    fieldAvg, fieldVar
       end if

        if (allocated(xPoints))   deallocate(xPoints)
        if (allocated(randField)) deallocate(randField)
        if (allocated(seed))      deallocate(seed)
        deallocate(avg_Gauss, stdDev_Gauss)
        deallocate(avg_Trans, stdDev_Trans)

    end subroutine func_test_004_conv_Nmc

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine func_test_005_conv_corrL(xMin, xMax,                     &
                                        xStep, corrMod,                 &
                                        margiFirst, fieldAvg, fieldVar, &
                                        Nmc, method,                    &
                                        corrL_initial, corrL_nIter,     &
                                        corrL_mult, corrL_add, corrL_pointsPerCorrL, &
                                        xMinGlob, xMaxGlob)

        implicit none

        !INPUTS
        double precision, dimension(:), intent(inout) :: xMin, xMax;
        double precision, dimension(:), intent(in) :: xStep;
        character (len=*)             , intent(in) :: corrMod, margiFirst;
        double precision              , intent(in) :: fieldAvg, fieldVar;
        integer                       , intent(in) :: Nmc, method;
        integer                       , intent(in) :: corrL_nIter
        double precision, dimension(:), intent(in) :: corrL_mult, corrL_add, corrL_initial
        integer                       , intent(in), optional :: corrL_pointsPerCorrL
        double precision, dimension(:), intent(in), optional :: xMinGlob, xMaxGlob
        integer, parameter :: seedStart = 3

        !LOCAL
        integer :: nDim;
        integer :: i;
        integer :: comm, rang, code;
        double precision, dimension(:, :), allocatable :: randField, xPoints;
        double precision :: t1, t2, t3, t4;
        double precision  , dimension(:) , allocatable :: avg_Gauss, stdDev_Gauss, &
                                                          avg_Trans, stdDev_Trans
        double precision  , dimension(:,:) , allocatable :: avg_Gauss_evnt, stdDev_Gauss_evnt, &
                                                            avg_Trans_evnt, stdDev_Trans_evnt
        double precision, dimension(:), allocatable :: corrL;
        double precision :: multInd
        integer, dimension(:), allocatable :: seed
        double precision :: testRand
        character (len=200) :: corrL_result_path

        corrL_result_path = string_vec_join([results_path,"/",results_folder_name,"/",corrL_folder_name])


        nDim = size(xMax)
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK (comm ,rang,code)

        allocate(corrL(nDim))
        allocate(avg_Gauss(corrL_nIter), stdDev_Gauss(corrL_nIter))
        allocate(avg_Trans(corrL_nIter), stdDev_Trans(corrL_nIter))
        allocate(avg_Gauss_evnt(corrL_nIter, Nmc), stdDev_Gauss_evnt(corrL_nIter, Nmc))
        allocate(avg_Trans_evnt(corrL_nIter, Nmc), stdDev_Trans_evnt(corrL_nIter, Nmc))

        do i = 1, corrL_nIter

            if(rang == 0) write(*,*) ""
            if(rang == 0) write(*,*) "-------------------------------------------------------------------------------";
            if(rang == 0) write(*,*) "ITERATION corrL",       i, "---------------------------------------------------";
            if(rang == 0) write(*,*) "-------------------------------------------------------------------------------";

            !Putting everything to the initial condition
            if (allocated(xPoints))   deallocate(xPoints)
            if (allocated(randField)) deallocate(randField)
            call calculate_random_seed(seed, seedStart)

            corrL = corrL_initial         &
                    *(corrL_mult)**(i-1) &
                    + corrL_add*(i-1)

            if(rang == 0) write(*,*) ""
            if(rang == 0) write(*,*) "  corrL = ", corrL;
            if(rang == 0) write(*,*) ""

            if(rang == 0) write(*,*) "-> Creating xPoints";
            if(present(corrL_pointsPerCorrL)) then
                call set_XPoints_perCorrL(xMin, xMax, corrL_pointsPerCorrL, corrL, xPoints, xMinGlob, xMaxGlob, .true.)
            else
                call cutBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
                call set_XPoints(xMin, xMax, xStep, xPoints)
                call restoreBorders(xMin, xMinGlob , xMax, xMaxGlob, xStep)
            end if
            allocate(randField(size(xPoints,2), Nmc))

            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) RAW")
            !call dispCarvalhol(transpose(xPoints(:, :)), "transpose(xPoints(:, :)) BEG")
            !call dispCarvalhol(transpose(xPoints(:, 1:20)), "transpose(xPoints(:, 1:20)) BEG")
            !call dispCarvalhol(transpose(xPoints(:, size(xPoints,2)-20:)), "transpose(xPoints(:, size(xPoints,2)-20:)) END")

            call write_generation_spec(xMin, xMax, xStep,              &
                                       corrL, corrMod,                 &
                                       margiFirst, fieldAvg, fieldVar, &
                                       Nmc, method, seed,              &
                                       rang, corrL_result_path, stringNumb_join("iteration_",i))

            t1 = MPI_Wtime();

            if(rang == 0) write(*,*) "-> Creating Standard Gaussian Random field (unstructured)";
            randField = 0.0
            call create_Std_Gaussian_Field_Unstruct (xPoints, corrL, corrMod, Nmc,  &
                                                      randField, method, seed)
            t2 = MPI_Wtime();

            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            call write_Mono_XMF_h5(xPoints, randField(:,:), "gauss_", rang, corrL_result_path, &
                                            MPI_COMM_WORLD, ["it_", "_proc_"], [i, rang], i)

            !call write_MatlabTable(randField, "visu/"//trim(stringNumb_join("x_Step_Iteration_",i)))
            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) GAUSSIAN")


            if(rang == 0) write(*,*) "-> Calculating Statistics Before Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                  &
                                                  avg_Gauss(i), stdDev_Gauss(i),              &
                                                  comm,                                       &
                                                  avg_Gauss_evnt(i,:), stdDev_Gauss_evnt(i,:))

            t3 = MPI_Wtime();

            if(rang == 0) write(*,*) "-> Transforming Random field (unstructured)";
            call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

            t4 = MPI_Wtime();

            write(get_fileId(),*) "-> Writing XMF and hdf5 files";
            call write_Mono_XMF_h5(xPoints, randField(:,:), "trans_", rang, corrL_result_path, &
                                            MPI_COMM_WORLD, ["it_", "_proc_"], [i, rang], i)

            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) TRANSFORMED")

            if(rang == 0) write(*,*) "-> Calculating Statistics After Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                   &
                                                  avg_Trans(i), stdDev_Trans(i),               &
                                                  comm,                                        &
                                                  avg_Trans_evnt(i,:), stdDev_Trans_evnt(i,:))

        end do

       if(rang == 0) then
           write(*,*) "-----CONVERGENCE TEST CHANGING corrL-----"
           write(*,*) ""
           write(*,*) " Initial corrL = ", corrL_initial
           write(*,*) "     Iter Mult = ", corrL_mult
           write(*,*) "     Iter Add  = ", corrL_add
           write(*,*) ""
           write(*,*) "     By Event "
           call DispCarvalhol(avg_Gauss_evnt   , "avg_Gauss_evnt")
           call DispCarvalhol(stdDev_Gauss_evnt, "stdDev_Gauss_evnt")
           write(*,*) ""
           write(*,*) "     GLOBAL "
           write(*,*) ""
           write(*,'(A4, A10, A2, 4A20)') "Iter", "nCorrL", " )", "Avg-Gauss", "StdDev-Gauss", &
                                                              "Avg-Trans", "StdDev-Trans"
           do i = 1, corrL_nIter
               multInd = product((xMax - xMin) / &
                                 (corrL_mult**(i-1) + corrL_add*(i-1)))
               write(*,'(I4, F10.1, A2, 4F20.8)') i, multInd , " )",              &
                                                  avg_Gauss(i), stdDev_Gauss(i), &
                                                  avg_Trans(i), stdDev_Trans(i)
           end do
           write(*,'(A16, 4F20.8)') "Reference      )", &
                                    0.0d0, 1.0d0,  &
                                    fieldAvg, fieldVar
       end if

        if (allocated(corrL))     deallocate(corrL)
        if (allocated(xPoints))   deallocate(xPoints)
        if (allocated(randField)) deallocate(randField)
        if (allocated(seed))      deallocate(seed)
        deallocate(avg_Gauss, stdDev_Gauss)
        deallocate(avg_Trans, stdDev_Trans)
        deallocate(avg_Gauss_evnt, stdDev_Gauss_evnt)
        deallocate(avg_Trans_evnt, stdDev_Trans_evnt)

    end subroutine func_test_005_conv_corrL

end module test_func_RF
