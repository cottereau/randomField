module test_func_RF

    use randomFieldND
    use statistics_RF
    use writeResultFile_RF
    use readFile_RF
    use displayCarvalhol
    use mesh_RF
    use mpi

    implicit none



contains

    subroutine func_test_001_conv_xStep(xMin, xMax,                      &
                                        corrL, corrMod,                  &
                                        margiFirst, fieldAvg, fieldVar,  &
                                        Nmc, method,                     &
                                        step_initial, step_nIter,        &
                                        step_mult, step_add)

        implicit none

        !INPUTS
        double precision, dimension(:), intent(in) :: xMin, xMax;
        double precision, dimension(:), intent(in) :: corrL;
        character (len=*)             , intent(in) :: corrMod, margiFirst;
        double precision              , intent(in) :: fieldAvg, fieldVar;
        integer                       , intent(in) :: Nmc, method;
        integer                       , intent(in) :: step_nIter
        double precision, dimension(:), allocatable, intent(in) :: step_mult, step_add, step_initial
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

        nDim = size(xMax)
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK (comm ,rang,code)

        allocate(xStep(nDim))
        allocate(avg_Gauss(step_nIter), stdDev_Gauss(step_nIter))
        allocate(avg_Trans(step_nIter), stdDev_Trans(step_nIter))
        allocate(avg_Gauss_evnt(step_nIter, Nmc), stdDev_Gauss_evnt(step_nIter, Nmc))
        allocate(avg_Trans_evnt(step_nIter, Nmc), stdDev_Trans_evnt(step_nIter, Nmc))

        do i = 1, step_nIter

            if(rang == 0) write(*,*) ""
            if(rang == 0) write(*,*) "-------------------------------------------------------------------------------";
            if(rang == 0) write(*,*) "ITERATION xStep",       i, "---------------------------------------------------";
            if(rang == 0) write(*,*) "-------------------------------------------------------------------------------";

            !Putting everything to the initial condition
            if (allocated(xPoints))   deallocate(xPoints)
            if (allocated(randField)) deallocate(randField)
            call calculate_random_seed(seed, seedStart)

            xStep = step_initial         &
                    *(step_mult)**(i-1) &
                    + step_add*(i-1)

            if(rang == 0) write(*,*) ""
            if(rang == 0) write(*,*) "  xStep = ", xStep;
            if(rang == 0) write(*,*) ""

            if(rang == 0) write(*,*) "-> Creating xPoints";
            call set_XPoints(xMin, xMax, xStep, xPoints)
            allocate(randField(size(xPoints,2), Nmc))

            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) RAW")
            !call dispCarvalhol(transpose(xPoints(:, 1:20)), "transpose(xPoints(:, 1:20)) BEG")
            !call dispCarvalhol(transpose(xPoints(:, size(xPoints,2)-20:)), "transpose(xPoints(:, size(xPoints,2)-20:)) END")

            t1 = MPI_Wtime();

            if(rang == 0) write(*,*) "-> Creating Standard Gaussian Random field (unstructured)";
            randField = 0.0
            call create_Std_Gaussian_Field_Unstruct (xPoints, corrL, corrMod, Nmc,  &
                                                      randField, method, seed)

            call write_MatlabTable(randField, "visu/"//trim(stringNumb_join("x_Step_Iteration_",i)))
            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) GAUSSIAN")

            t2 = MPI_Wtime();

            if(rang == 0) write(*,*) "-> Calculating Statistics Before Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                  &
                                                  avg_Gauss(i), stdDev_Gauss(i),              &
                                                  comm,                                       &
                                                  avg_Gauss_evnt(i,:), stdDev_Gauss_evnt(i,:))

            t3 = MPI_Wtime();

            if(rang == 0) write(*,*) "-> Transforming Random field (unstructured)";
            call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

            t4 = MPI_Wtime();

            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) TRANSFORMED")

            if(rang == 0) write(*,*) "-> Calculating Statistics After Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                   &
                                                  avg_Trans(i), stdDev_Trans(i),               &
                                                  comm,                                        &
                                                  avg_Trans_evnt(i,:), stdDev_Trans_evnt(i,:))

        end do

       if(rang == 0) then
           write(*,*) "-----CONVERGENCE TEST CHANGING xStep-----"
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

    end subroutine func_test_001_conv_xStep

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine func_test_002_conv_Nmc(xMin, xMax,                      &
                                      corrL, corrMod,                  &
                                      margiFirst, fieldAvg, fieldVar,  &
                                      xStep, method,                   &
                                      nmc_initial, nmc_nIter,          &
                                      nmc_mult, nmc_add)

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

        nDim = size(xMax)
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK (comm ,rang,code)

        allocate(avg_Gauss(nmc_nIter), stdDev_Gauss(nmc_nIter))
        allocate(avg_Trans(nmc_nIter), stdDev_Trans(nmc_nIter))

        if(rang == 0) write(*,*) "-> Creating xPoints";
        call set_XPoints(xMin, xMax, xStep, xPoints)

        do i = 1, nmc_nIter

            if(rang == 0) write(*,*) ""
            if(rang == 0) write(*,*) "-------------------------------------------------------------------------------";
            if(rang == 0) write(*,*) "ITERATION Nmc",       i, "-----------------------------------------------------";
            if(rang == 0) write(*,*) "-------------------------------------------------------------------------------";

            !Putting everything to the initial condition
            if (allocated(randField)) deallocate(randField)
            call calculate_random_seed(seed, seedStart)

            Nmc = nmc_initial         &
                    *(nmc_mult)**(i-1) &
                    + nmc_add*(i-1)

            if(rang == 0) write(*,*) ""
            if(rang == 0) write(*,*) "  Nmc = ", Nmc;
            if(rang == 0) write(*,*) ""

            allocate(randField(size(xPoints,2), Nmc))

            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) RAW")
            !call dispCarvalhol(transpose(xPoints(:, 1:20)), "transpose(xPoints(:, 1:20)) BEG")
            !call dispCarvalhol(transpose(xPoints(:, size(xPoints,2)-20:)), "transpose(xPoints(:, size(xPoints,2)-20:)) END")

            t1 = MPI_Wtime();

            if(rang == 0) write(*,*) "-> Creating Standard Gaussian Random field (unstructured)";
            randField = 0.0
            call create_Std_Gaussian_Field_Unstruct (xPoints, corrL, corrMod, Nmc,  &
                                                      randField, method, seed)

            call write_MatlabTable(randField, "visu/"//trim(stringNumb_join("Nmc_Iteration_",i)))
            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) GAUSSIAN")

            t2 = MPI_Wtime();

            if(rang == 0) write(*,*) "-> Calculating Statistics Before Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                  &
                                                  avg_Gauss(i), stdDev_Gauss(i),              &
                                                  comm)

            t3 = MPI_Wtime();

            if(rang == 0) write(*,*) "-> Transforming Random field (unstructured)";
            call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

            t4 = MPI_Wtime();

            !call dispCarvalhol(randField(1:10, :), "randField(1:10, :) TRANSFORMED")

            if(rang == 0) write(*,*) "-> Calculating Statistics After Transformation";
            call calculate_average_and_stdVar_MPI(randField,                                   &
                                                  avg_Trans(i), stdDev_Trans(i),               &
                                                  comm)

        end do

       if(rang == 0) then
           write(*,*) "-----CONVERGENCE TEST CHANGING Nmc-----"
           write(*,*) ""
           write(*,*) "     GLOBAL "
           write(*,*) ""
           write(*,'(A4, A10, A2, 4A20)') "Iter", "multInd", " )", "Avg-Gauss", "StdDev-Gauss", &
                                                              "Avg-Trans", "StdDev-Trans"
           do i = 1, nmc_nIter
               multInd = nmc_mult**(i-1) + nmc_add*(i-1)
               write(*,'(I4, F10.0, A2, 4F20.8)') i, multInd , " )",             &
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

    end subroutine func_test_002_conv_Nmc

end module test_func_RF
