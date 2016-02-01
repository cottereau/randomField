module type_inputRF

    use mpi
    use readFile_RF
    use readUNV_RF

    implicit none

    type :: IPT_RF
        !LOCAL
        integer :: log_ID = -1
        integer :: comm = -1
        integer :: rang = -1
        integer :: nb_procs = -1
        logical :: init=.false.
        logical :: writeDataSet=.true.
        !MESH
        integer :: nDim_mesh
        integer :: meshMod
        double precision, dimension(:), allocatable :: xMaxGlob, xMinGlob;
        integer         , dimension(:), allocatable :: pointsPerCorrL;
        integer         , dimension(:), allocatable :: procPerDim;
        !UNV
        double precision, dimension(:,:), allocatable :: coordList
        integer         , dimension(:,:), allocatable :: connectList
        logical :: monotype
        logical :: unv = .false.
        character(len=1024) :: unv_path

        !GENERATION
        integer :: nDim_gen
        double precision   :: fieldAvg = -1, fieldVar = -1;
        double precision, dimension(:), allocatable :: corrL, overlap
        integer :: corrMod = -1 !1 for Gaussian
        integer :: margiFirst = -1 !1 for Gaussian, 2 for Lognormal
        integer :: method = -1 !1 for Isotropic, 2 for Shinozuka, 3 for Randomization, 4 for FFT
        integer :: Nmc = -1, seedStart
        logical :: independent

    end type IPT_RF

contains
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_IPT_RF(IPT, nDim, log_file_RF_ID, rang)

            implicit none

            !INPUT
            integer, intent(in) :: nDim
            integer, intent(in) :: log_file_RF_ID, rang
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT

            if(.not. allocated(IPT%xMaxGlob)) allocate(IPT%xMaxGlob(nDim))
            if(.not. allocated(IPT%xMinGlob)) allocate(IPT%xMinGlob(nDim))
            if(.not. allocated(IPT%pointsPerCorrL)) allocate(IPT%pointsPerCorrL(nDim))
            if(.not. allocated(IPT%corrL)) allocate(IPT%corrL(nDim))
            if(.not. allocated(IPT%overlap)) allocate(IPT%overlap(nDim))
            if(.not. allocated(IPT%procPerDim)) allocate(IPT%procPerDim(nDim))

            IPT%log_ID = log_file_RF_ID
            IPT%rang   = rang
            IPT%init   = .true.
            IPT%xMaxGlob   = -1.0D0
            IPT%xMinGlob = -1.0D0
            IPT%pointsPerCorrL   = -1.0D0
            IPT%corrL   = -1.0D0
            IPT%overlap = -1.0D0

        end subroutine init_IPT_RF

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine finalize_IPT_RF(IPT)

            implicit none

            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT

            if(allocated(IPT%xMaxGlob)) deallocate(IPT%xMaxGlob)
            if(allocated(IPT%xMinGlob)) deallocate(IPT%xMinGlob)
            if(allocated(IPT%pointsPerCorrL)) deallocate(IPT%pointsPerCorrL)
            if(allocated(IPT%corrL)) deallocate(IPT%corrL)
            if(allocated(IPT%overlap)) deallocate(IPT%overlap)
            if(allocated(IPT%procPerDim)) deallocate(IPT%procPerDim)

        end subroutine finalize_IPT_RF

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_mesh_input(path, IPT)

            implicit none
            !INPUT
            character(len=*), intent(in) :: path
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT
            !LOCAL
            character(len=1024) , dimension(:,:), allocatable :: dataTable;

            call set_DataTable(path, dataTable)

            call wLog("    set_DataTable")

            call read_DataTable(dataTable, "nDim", IPT%nDim_mesh)
            call read_DataTable(dataTable, "meshMod", IPT%meshMod)

            if(.not. IPT%init) call init_IPT_RF(IPT, IPT%nDim_mesh, IPT%log_ID, IPT%rang)

            select case (IPT%meshMod)
                case(1)
                    if(IPT%rang==0) write(*,*) "   Mesh automatic"
                    call wLog("    Mesh automatic")
                    call read_DataTable(dataTable, "Min", IPT%xMinGlob)
                    call read_DataTable(dataTable, "Max", IPT%xMaxGlob)
                    call read_DataTable(dataTable, "pointsPerCorrL", IPT%pointsPerCorrL)
                case(2)
                    !stop("Inside read_mesh_input UNV not updated")
                    if(IPT%rang==0) write(*,*) "   Mesh UNV"
                    call wLog("    Mesh UNV")
                    call wLog("        file: ")
                    !call wLog(trim(IPT%unv_path))
                    call read_DataTable(dataTable, "unv_path", IPT%unv_path)
                    call read_DataTable(dataTable, "pointsPerCorrL", IPT%pointsPerCorrL)
                    IPT%unv = .true.
                    call readUNV(IPT%unv_path, IPT%nDim_mesh, IPT%coordList, IPT%connectList, IPT%monotype, &
                                 IPT%rang, IPT%nb_procs, IPT%comm)
                    !call readUNV_many([IPT%unv_path, IPT%unv_path], IPT%nDim_mesh, IPT%coordList, &
                    !                  IPT%connectList, IPT%monotype, &
                    !                  IPT%rang, IPT%nb_procs, IPT%comm)
                    call wLog("-> defining_UNV_extremes")
                    if(IPT%rang==0) write(*,*) "-> defining_UNV_extremes"
                    call get_Global_Extremes_Mesh(IPT%coordList, IPT%comm, IPT%xMinGlob, IPT%xMaxGlob)
                    if(IPT%rang==0) write(*,*) " IPT%xMinGlob = ", IPT%xMinGlob
                    if(IPT%rang==0) write(*,*) " IPT%xMaxGlob = ", IPT%xMaxGlob
!                    call DispCarvalhol(transpose(IPT%connectList), "transpose(IPT%connectList)", &
!                                       nColumns=8, unit_in=IPT%log_ID)
!                    call DispCarvalhol(transpose(IPT%coordList), "transpose(IPT%coordList)", &
!                                       unit_in=IPT%log_ID)
                case default
                    write(*,*) "meshMod not accepted: ", IPT%meshMod
                    stop(" ")
            end select


            if(allocated(dataTable)) deallocate(dataTable)

        end subroutine read_mesh_input


        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_generation_input(path, IPT)

            implicit none
            !INPUT
            character(len=*), intent(in) :: path
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT
            !LOCAL
            character(len=1024) , dimension(:,:), allocatable :: dataTable;
            integer :: independent

            call set_DataTable(path, dataTable)
            call read_DataTable(dataTable, "nDim", IPT%nDim_gen)

            if(.not. IPT%init) call init_IPT_RF(IPT, IPT%nDim_gen, IPT%log_ID, IPT%rang)

            call read_DataTable(dataTable, "Nmc", IPT%Nmc)
            call read_DataTable(dataTable, "corrMod"    , IPT%corrMod)
            call read_DataTable(dataTable, "margiFirst" , IPT%margiFirst)
            call read_DataTable(dataTable, "fieldAvg"   , IPT%fieldAvg)
            call read_DataTable(dataTable, "fieldVar"   , IPT%fieldVar)
            call read_DataTable(dataTable, "method"     , IPT%method)
            call read_DataTable(dataTable, "seedStart"  , IPT%seedStart)
            call read_DataTable(dataTable, "corrL"      , IPT%corrL)
            call read_DataTable(dataTable, "independent", independent)

            if(independent == 1) then
                IPT%independent = .true.
                call read_DataTable(dataTable, "overlap", IPT%overlap)
            else
                IPT%independent = .false.
            end if

            deallocate(dataTable)

        end subroutine read_generation_input

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine validate_input(IPT)

            !OUTPUT
            type(IPT_RF), intent(in)  :: IPT
            !LOCAL
            integer :: i

            !Input Validation
            if(IPT%Nmc < 1) then
               write(*,*) ""
               write(*,*) "ERROR - Number of events should be a positive integer"
               write(*,*) "IPT%Nmc = ", IPT%Nmc
               stop(" ")
            end if

            if((IPT%nDim_gen < 1) .or. (IPT%nDim_mesh < 1) .or. (IPT%nDim_gen /= IPT%nDim_mesh)) then
               write(*,*) ""
               write(*,*) "ERROR - nDim should be a positive integer and should be the same in both mesh and generation files"
               write(*,*) "IPT%nDim_gen  = ", IPT%nDim_gen
               write(*,*) "IPT%nDim_mesh = ", IPT%nDim_mesh
               stop(" ")
            end if

            do i = 1, IPT%nDim_gen
               if(IPT%corrL(i) <= 0.0d0) then
                   write(*,*) ""
                   write(*,*) "ERROR - corrL should be a positive number greater than 0.0"
                   stop(" ")
               end if
            end do

        end subroutine validate_input

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_IPT_RF(IPT, name, forLog_in, unit_in)
            !INPUT
            type(IPT_RF), intent(in) :: IPT
            character(len=*), intent(in), optional :: name
            integer, intent(in), optional :: unit_in
            logical, intent(in), optional :: forLog_in
            !LOCAL
            character(len = 20) :: dblFmt, matDblFmt, intFmt, matIntFmt
            integer :: unit
            logical :: forLog
            logical :: active

            active = .true.
            unit = 6 !Screen
            if(present(unit_in)) unit = unit_in

            forLog = .false.
            if(present(forLog_in)) then
#ifdef MAKELOG
                if(forLog_in) unit = IPT%log_ID
                forLog = forLog_in
#else
                if(forLog_in) write(*,*) "WARNING!!! Inside show_IPT_RF, forLog_in = .true. but MAKELOG was not defined"
                active = .false.
#endif
            end if

            if(active) then

                write(unit,*) "INPUT------------------------------------------------------------"
                if(present(name)) write(unit,*) "  ", name
                write(unit,*) " "
                write(unit,*) " log_ID = ", IPT%log_ID
                write(unit,*) " rang = ", IPT%rang
                write(unit,*) " init = ", IPT%init
                write(unit,*) " "
                write(unit,*) " MESH -----------------"
                write(unit,*) " nDim_mesh = ", IPT%nDim_mesh
                write(unit,*) " meshMod = ", IPT%meshMod
                write(unit,*) " xMaxGlob = ", IPT%xMaxGlob
                write(unit,*) " xMinGlob = ", IPT%xMinGlob
                write(unit,*) " pointsPerCorrL = ", IPT%pointsPerCorrL
                write(unit,*) " unv = ", IPT%unv
                if(IPT%unv) then
                    write(unit,*) " unv_path = ", IPT%unv_path
                    write(unit,*) " shape(coordList) = "  , shape(IPT%coordList)
                    write(unit,*) " shape(connectList) = ", shape(IPT%connectList)
                    write(unit,*) " monotype = ", IPT%monotype
                end if
                write(unit,*) " "
                write(unit,*) " GENERATION -----------------"
                write(unit,*) " nDim_gen = ", IPT%nDim_gen
                write(unit,*) " fieldAvg = ", IPT%fieldAvg
                write(unit,*) " fieldVar = ", IPT%fieldVar
                write(unit,*) " corrL = ", IPT%corrL
                write(unit,*) " overlap = ", IPT%overlap
                write(unit,*) " corrMod = ", IPT%corrMod
                write(unit,*) " margiFirst = ", IPT%margiFirst
                write(unit,*) " method = ", IPT%method
                write(unit,*) " Nmc = ", IPT%Nmc
                write(unit,*) " independent = ", IPT%independent
                write(unit,*) " seedStart = ", IPT%seedStart
                write(unit,*) " "

            end if

        end subroutine show_IPT_RF

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_Global_Extremes_Mesh(coordList, comm, xMinGlob, xMaxGlob)
        implicit none

        !INPUT
        double precision, dimension(:,:), intent(in) :: coordList
        integer, intent(in) :: comm
        !OUTPUT
        double precision, dimension(:), intent(out) :: xMinGlob, xMaxGlob
        !LOCAL
        integer :: nDim
        double precision, dimension(:), allocatable :: xMinLoc, xMaxLoc
        integer :: code, i

        !call DispCarvalhol(transpose(coordList), "transpose(coordList)")

        nDim = size(coordList,1)
        allocate(xMinLoc(nDim))
        allocate(xMaxLoc(nDim))

        xMinLoc = minval(coordList, 2)
        xMaxLoc = maxval(coordList, 2)

        call wLog("xMinLoc = ")
        call wLog(xMinLoc)
        call wLog("xMaxLoc = ")
        call wLog(xMaxLoc)

        do i = 1, nDim
            call MPI_ALLREDUCE (xMinLoc(i), xMinGlob(i), 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm,code)
            call MPI_ALLREDUCE (xMaxLoc(i), xMaxGlob(i), 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm,code)
        end do

        call wLog(" ")
        call wLog("xMinGlob = ")
        call wLog(xMinGlob)
        call wLog("xMaxGlob = ")
        call wLog(xMaxGlob)

        if(allocated(xMinLoc)) deallocate(xMinLoc)
        if(allocated(xMaxLoc)) deallocate(xMaxLoc)

    end subroutine get_Global_Extremes_Mesh

end module type_inputRF
