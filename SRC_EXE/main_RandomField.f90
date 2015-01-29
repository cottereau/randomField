program main_RandomField

	use mpi
	use constants_RF
	use readFile_RF
	use randomFieldND
	use test_func_RF
	use statistics_RF
	use mesh_RF
	use writeResultFile_RF
	use displayCarvalhol
	use charFunctions

    implicit none

    !INPUTS
    integer                                        :: Nmc, nDim;
    character (len=30), parameter                  :: mesh_input = "mesh_input"
    character (len=30), parameter                  :: gen_input  = "gen_input"
    character (len=30), parameter                  :: test_input = "test_input"
    character (len=50)                             :: mode;
    character (len=15)                             :: corrMod, margiFirst, meshType, meshMod;
    double precision,   dimension(:),  allocatable :: corrL;
    double precision                               :: fieldAvg, fieldVar;
    logical                                        :: calculateStat = .true.
    integer                                        :: method !1 for Isotropic, 2 for Shinozuka
    integer         , parameter                    :: pointsPerCorrL = 10;
    logical :: step_variate, nmc_variate
    integer :: step_nIter, nmc_nIter, nmc_initial
    double precision, dimension(:), allocatable :: step_mult, step_add, step_initial
    double precision :: nmc_mult, nmc_add


	!OUTPUTS (in this case variables to be filled in in the proces)
	integer         , dimension(:)   , allocatable :: xNStep;
    double precision, dimension(:, :), allocatable :: randField, xPoints;
    double precision, dimension(:)   , allocatable :: xMax, xMin, xStep;

	!LOCAL VARIABLES
    logical            :: file_exist
    integer            :: i, baseStep;
    integer            :: comm, code, rang, error, nb_procs;
    character(len=30)  :: rangChar;
    character(len=200) :: path
    character(len=110), dimension(:)  , allocatable :: HDF5Name
    character(len=30) , dimension(:,:), allocatable :: dataTable;

    !Initializing MPI
    if(rang == 0) write(*,*) "-> Initialize MPI_communications"
    call init_communication()
    comm = MPI_COMM_WORLD
    call MPI_COMM_RANK (comm, rang, code)
    call MPI_COMM_SIZE (comm, nb_procs, code)

    if(rang == 0)then
        write(*,*)
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "************************                             ***********************"
        write(*,*) "************************  RANDOM FIELD LIBRARY TEST  ***********************"
        write(*,*) "************************                             ***********************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*)
    end if

	if(rang == 0) write(*,*) "  Number of procs = ", nb_procs
	if(rang == 0) write(*,*) ""

    !Reading Inputs
    if(rang == 0) write(*,*) "-> Reading inputs"
    call read_input()

    !Initial allocation
    call allocate_init()

    !Generating random fields
    if(rang == 0) write(*,*) "-> Generating random fields"

    if(step_variate) then
        if(rang == 0) write(*,*) "-> Running xStep variation tests"
        call func_test_001_conv_xStep(xMin, xMax,                     &
                                      corrL, corrMod,                 &
                                      margiFirst, fieldAvg, fieldVar, &
                                      Nmc, method,                    &
                                      step_initial, step_nIter,       &
                                      step_mult, step_add)

    end if
    if(nmc_variate) then
        if(rang == 0) write(*,*) "-> Running Nmc variation tests"
        call func_test_002_conv_Nmc(xMin, xMax,                      &
                                    corrL, corrMod,                  &
                                    margiFirst, fieldAvg, fieldVar,  &
                                    xStep, method,                   &
                                    nmc_initial, nmc_nIter,          &
                                    nmc_mult, nmc_add)
    end if

	!Deallocating
	call deallocate_all()


	if(rang == 0) then
	    write(*,*) ""
		write(*,*) "---------------------------------------------------------------------";
	    write(*,*) "-----------------END RANDOM FIELD LIBRARY TEST-----------------------";
		write(*,*) "---------------------------------------------------------------------";
        write(*,*) ""
	end if

    !Finalizing MPI
	call end_communication()

        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
    contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_communication()
            call MPI_INIT(code)
            call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
        end subroutine init_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine allocate_init()
        end subroutine allocate_init

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine deallocate_all()
            if(allocated(corrL))   deallocate(corrL);
            if(allocated(xMax))    deallocate(xMax);
            if(allocated(xMin))    deallocate(xMin);
            if(allocated(xStep))   deallocate(xStep);
            if(allocated(xPoints))   deallocate(xPoints);
            if(allocated(randField)) deallocate(randField);

            if(allocated(step_mult)) deallocate(step_mult);
            if(allocated(step_add)) deallocate(step_add);

        end subroutine deallocate_all

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine end_communication()
            call MPI_FINALIZE(code)
        end subroutine end_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_input()

            !Reading Generation Statistis Input
            write(*,*) "    Reading Generation Statistics Input"
            !path = trim(adjustL(sampleSpecFolder))//"/"//sampleSpecName
            path = gen_input
            path = adjustL(path)
            write(*,*) "        file: ", trim(path)
            inquire(file=path, exist = file_exist)
            if(.not. file_exist) then
                write(*,*) "ERROR - The file ", path, " was not found"
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            call set_DataTable(path, dataTable)
            !if(rang == 0) call DispCarvalhol (dataTable, path);
            call read_DataTable(dataTable, "Nmc"        ,  Nmc)
            call read_DataTable(dataTable, "corrMod"    , corrMod)
            call read_DataTable(dataTable, "margiFirst" , margiFirst)
            call read_DataTable(dataTable, "corrL"      , corrL)
            call read_DataTable(dataTable, "fieldAvg"   , fieldAvg)
            call read_DataTable(dataTable, "fieldVar"   , fieldVar)
            call read_DataTable(dataTable, "method"     , method)
            deallocate(dataTable)

            !Allocation
            nDim = size(corrL)
            allocate (xMax(nDim))
            allocate (xMin(nDim))
            allocate (xStep(nDim))

            !Reading Mesh
            write(*,*) "    Reading Mesh Input"
            path = mesh_input
            path = adjustL(path)
            write(*,*) "        file: ", trim(path)
            inquire(file=path, exist = file_exist)
            if(.not. file_exist) then
                write(*,*) "ERROR - The file ", path, " was not found"
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if

            call set_DataTable(path, dataTable)
            !if(rang == 0) call DispCarvalhol (dataTable, path);
            call read_DataTable(dataTable, "meshType", meshType)
            call read_DataTable(dataTable, "meshMod", meshMod)

            if (meshType == "structured" .or. meshType == "unstructured") then
               select case (meshMod)
                   case("manual")
                       write(rangChar,'(I7)') rang
                       call read_DataTable(dataTable, "Max"//trim(adjustl(rangChar)), xMax)
                       call read_DataTable(dataTable, "Min"//trim(adjustl(rangChar)), xMin)
                   case("automatic")

                        baseStep = nint(dble(nb_procs)**(1.0d0/nDim))
                        if (nb_procs == baseStep**nDim) then
                            call read_DataTable(dataTable, "Max", xMax)
                            call read_DataTable(dataTable, "Min", xMin)
                            call read_DataTable(dataTable, "Step", xStep)
                            call set_Local_Extremes_Mesh(xMin, xMax, rang, nb_procs)
                        else
                            write(*,*) "In meshing automatic mode the number of processors should be an integer like (N)**nDim, where N is an integer and nDim is the number of dimensions"
                            call MPI_ABORT(MPI_COMM_WORLD, error, code)
                        end if
                end select
            else
                write(*,*) "meshType not accepted: ", meshType
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            deallocate(dataTable)

            !Looking for test_input
            write(*,*) "    Checking for Test Input"
            step_variate = .false.
            nmc_variate  = .false.
            path = test_input
            path = adjustL(path)
            write(*,*) "        file: ", trim(path)
            inquire(file=path, exist = file_exist)
            if(file_exist) then
                call set_DataTable(path, dataTable)
                !if(rang == 0) call DispCarvalhol (dataTable, path);
                call read_DataTable(dataTable, "step_variate", step_variate)
                call read_DataTable(dataTable, "nmc_variate", nmc_variate)
                if(step_variate .and. nmc_variate) then
                    write(*,*) "You ask for both Step and Nmc variations, only Step variations will be applied"
                end if
                if(step_variate) then
                    allocate(step_mult(nDim))
                    allocate(step_add(nDim))
                    allocate(step_initial(nDim))
                    call read_DataTable(dataTable, "step_initial", step_initial)
                    call read_DataTable(dataTable, "step_mult", step_mult)
                    call read_DataTable(dataTable, "step_add", step_add)
                    call read_DataTable(dataTable, "step_nIter", step_nIter)
                end if
                if(nmc_variate) then
                    call read_DataTable(dataTable, "nmc_initial", nmc_initial)
                    call read_DataTable(dataTable, "nmc_mult", nmc_mult)
                    call read_DataTable(dataTable, "nmc_add", nmc_add)
                    call read_DataTable(dataTable, "nmc_nIter", nmc_nIter)
                end if
                deallocate(dataTable)
            else
                write(*,*) "The file ", path, " doesn't exist (no tests will be performed)"
            end if

            !Showing read inputs
            if(rang == 0) then
               write(*,*) ""
               write(*,*) "     MESH"
               write(*,*) "     |meshType = ", meshType
               write(*,*) "     |meshMod  = ", meshMod
               write(*,*) "     |xMin     = ", xMin
               write(*,*) "     |xMax     = ", xMax
               write(*,*) "     |xStep    = ", xStep
               write(*,*) ""
               write(*,*) "     EVENTS"
               write(*,*) "     |nDim       = ", nDim
               write(*,*) "     |Nmc        = ", Nmc !number of Monte-Carlo experiments
               write(*,*) "     |corrMod    = ", corrMod
               write(*,*) "     |margiFirst = ", margiFirst
               write(*,*) "     |corrL      = ", corrL
               write(*,*) "     |fieldAvg   = ", fieldAvg
               write(*,*) "     |fieldVar   = ", fieldVar
               write(*,*) "     |method     = ", method
               write(*,*) ""
               if(step_variate) then
                   write(*,*) "     TESTS STEP"
                   write(*,*) "     |step_initial = ", step_initial
                   write(*,*) "     |step_mult    = ", step_mult !number of Monte-Carlo experiments
                   write(*,*) "     |step_add     = ", step_add
                   write(*,*) "     |step_nIter   = ", step_nIter
               end if
               if(nmc_variate) then
                   if(step_variate) then
                        write(*,*) ""
                        write(*,*) "     TESTS NMC (will not be used)"
                   else
                        write(*,*) "     TESTS NMC"
                   end if
                   write(*,*) "     |nmc_initial = ", nmc_initial
                   write(*,*) "     |nmc_mult    = ", nmc_mult !number of Monte-Carlo experiments
                   write(*,*) "     |nmc_add     = ", nmc_add
                   write(*,*) "     |nmc_nIter   = ", nmc_nIter
               end if
            end if

            !Input Validation
            if(Nmc < 1) then
               write(*,*) ""
               write(*,*) "ERROR - Number of events should be a positive integer"
               call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            if(nDim < 1) then
               write(*,*) ""
               write(*,*) "ERROR - nDim should be a positive integer"
               write(*,*) "nDim           = ", nDim
               call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            do i = 1, nDim
               if(corrL(i) <= 0.0d0) then
                   write(*,*) ""
                   write(*,*) "ERROR - corrL should be a positive number greater than 0.0"
                   write(*,*) "corrL ", i, "of proc ", rang, " is ", corrL(i)
                   call MPI_ABORT(MPI_COMM_WORLD, error, code)
               end if
            end do

        end subroutine read_input

end program main_RandomField
