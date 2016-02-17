module localization_RF

    use displayCarvalhol
    use math_RF
    use constants_RF
    use mpi
    use write_Log_File
    use type_RF
    use type_MESH
    use common_variables_RF
    use randomFieldND

    implicit none

contains

    !    !-----------------------------------------------------------------------------------------------
    !    !-----------------------------------------------------------------------------------------------
    !    !-----------------------------------------------------------------------------------------------
    !    !-----------------------------------------------------------------------------------------------
    !    subroutine get_neighbours_info (RDF, MSH)
    !        implicit none
    !
    !        !INPUT OUTPUT
    !        type(RF), intent(inout) :: RDF
    !
    !        !INPUT
    !        type(MESH), intent(in) :: MSH
    !
    !        !LOCAL
    !        integer :: neighPos, stage, direction
    !        integer, allocatable, dimension(:) :: tempSeed
    !        integer :: request, code, tag
    !        integer, dimension(MPI_STATUS_SIZE) :: status
    !
    !        do neighPos = 1, size(MSH%neigh)
    !            if(MSH%neigh(neighPos) < 0) cycle !Check if this neighbour exists
    !            call calculate_random_seed(tempSeed, RDF%seedStart+MSH%neigh(neighPos))
    !            RDF%neighSeed(:,neighPos) = tempSeed
    !        end do
    !
    !        do stage = 1, 2 !Sending and then receiving
    !
    !            do neighPos = 1, size(MSH%neigh)
    !
    !                if(MSH%neigh(neighPos) < 0) cycle !Check if this neighbour exists
    !
    !                if(stage == 1) then
    !                    tag = findTag(MSH, neighPos, neighPos, send = .true.)
    !                    call MPI_ISEND (RDF%xRange(:), RDF%nDim, MPI_DOUBLE_PRECISION, &
    !                        MSH%neigh(neighPos), tag, RDF%comm, request, code)
    !                else if(stage == 2) then
    !                    tag = findTag(MSH, neighPos, neighPos, send = .false.)
    !                    call MPI_RECV (RDF%neighRange(:,neighPos), RDF%nDim, MPI_DOUBLE_PRECISION, &
    !                        MSH%neigh(neighPos), tag, RDF%comm, status, code)
    !                end if
    !
    !            end do
    !
    !        end do
    !
    !        deallocate(tempSeed)
    !
    !    end subroutine get_neighbours_info

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine getNeighIndexRange(MSH, minIndexNeigh, maxIndexNeigh)
        !INPUT
        type(MESH), intent(in) :: MSH
        !OUTPUT
        integer, intent(out) :: minIndexNeigh, maxIndexNeigh
        integer :: i
        !logical, dimension(:), intent(out) :: considerNeighbour

        !considerNeighbour = .true.
        !where(MSH%neigh < 0) considerNeighbour = .false.
        !where(minval(MSH%neighShift,1) < 0) considerNeighbour = .false.

        !Global Min and Max positions
        !minIndexNeigh = minval(pack(MSH%indexNeigh(1,:), MSH%considerNeighbour))
        !maxIndexNeigh = maxval(pack(MSH%indexNeigh(2,:), MSH%considerNeighbour))

        minIndexNeigh = -1
        maxIndexNeigh = -1
        do i = 1, size(MSH%indexNeigh,2)
            if(MSH%considerNeighbour(i)) then
                if(minIndexNeigh == -1) then
                    minIndexNeigh = MSH%indexNeigh(1,i)
                    maxIndexNeigh = MSH%indexNeigh(2,i)
                end if
                if(minIndexNeigh > MSH%indexNeigh(1,i)) minIndexNeigh = MSH%indexNeigh(1,i)
                if(maxIndexNeigh < MSH%indexNeigh(2,i)) maxIndexNeigh = MSH%indexNeigh(2,i)
            end if
        end do

    end subroutine getNeighIndexRange

    !    !-----------------------------------------------------------------------------------------------
    !    !-----------------------------------------------------------------------------------------------
    !    !-----------------------------------------------------------------------------------------------
    !    !-----------------------------------------------------------------------------------------------
    !    subroutine applyWeightingFunctions(RDF, MSH, minIndexNeigh, maxIndexNeigh, partitionType)
    !        !INPUT:
    !        type(RF), intent(in) :: RDF
    !        type(MESH), intent(in) :: MSH
    !        integer, intent(in) :: minIndexNeigh, maxIndexNeigh
    !        !logical, dimension(:), intent(in) ::considerNeighbour
    !        integer, intent(in) :: partitionType
    !
    !        !LOCAL
    !        integer :: direction, minPos, maxPos, i
    !        double precision, dimension(minIndexNeigh:maxIndexNeigh) :: unityPartition
    !        double precision, dimension(MSH%nDim) :: originCorner
    !
    !        !Modify extremes of local Random Field-------------------------------------------------------
    !
    !        !Building Shape Functions in all directions
    !        do direction = 1, size(MSH%neigh)
    !
    !            if(.not. MSH%considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction
    !
    !            !Positions in temp Vector
    !            minPos = MSH%indexNeigh(1,direction)
    !            maxPos = MSH%indexNeigh(2,direction)
    !
    !            !Finding origin
    !            originCorner = MSH%xOrNeigh(:, direction)
    !            !originCorner = MSH%xMinNeigh(:, direction)
    !            !where(MSH%neighShift(:, direction) == -1) originCorner = MSH%xMaxNeigh(:, direction)
    !
    !
    !            !Shape Function Generation
    !            call generateUnityPartition(RDF%xPoints(:, minPos:maxPos), originCorner, MSH%overlap, &
    !                                        MSH%neighShift(:, direction), partitionType, &
    !                                        unityPartition(minPos:maxPos))
    !
    !        end do !Direction
    !
    !        RDF%randField(minIndexNeigh:maxIndexNeigh,1) = RDF%randField(minIndexNeigh:maxIndexNeigh,1) &
    !                                                       * sqrt(unityPartition(minIndexNeigh:maxIndexNeigh))
    !
    !        !RDF%randField(minIndexNeigh:maxIndexNeigh,1) = unityPartition(minIndexNeigh:maxIndexNeigh) !TEST
    !
    !    end subroutine applyWeightingFunctions

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine addNeighboursFields(RDF, MSH)
        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF

        !INPUT
        type(MESH), intent(in) :: MSH

        !LOCAL
        integer :: i, direction, neighPos
        integer :: code
        integer, dimension(RDF%nDim) :: minPos, maxPos
        integer :: totalSize, double_size, overHead, overEst, bufferSize
        double precision, dimension(:), allocatable :: buffer
        integer :: testrank = 1, testrank2 = 11
        double precision, dimension(1:MSH%xNTotal), target :: tempRandField
        !integer, dimension(:), allocatable :: request
        integer, dimension(:)  , allocatable :: request
        integer, dimension(:,:), allocatable :: status
        integer :: requestSize, countReq, stage
        integer :: tag
        logical :: sndRcv
        double precision, dimension(:, :), pointer :: TRF_2D
        double precision, dimension(:, :, :), pointer :: TRF_3D


        call wLog ("Inside addNeighboursFields")

        if(RDF%nDim == 2) TRF_2D(1:MSH%xNStep(1),1:MSH%xNStep(2)) => tempRandField
        if(RDF%nDim == 3) TRF_3D(1:MSH%xNStep(1),1:MSH%xNStep(2),1:MSH%xNStep(3)) => tempRandField

        if(RDF%nDim == 1) then
            requestSize = 2* (2*1)
            stop("addNeighboursFields not yet in 1D")
        else if(RDF%nDim == 2) then
            requestSize = 2* (4*1 +  4*3)
        else if(RDF%nDim == 3) then
            requestSize = 2* (8*1 + 12*3 + 6*9)
        else
            stop("No requestSize for this dimension")
        end if

        countReq = 0
        tag = 0

        !Buffer allocation
        call wLog ("Allocating buffer")
        !overEst = MSH%nDim + 1
        overEst = 1
        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,double_size,code)
        overHead = int(1+(MPI_BSEND_OVERHEAD*1.)/double_size)
        bufferSize = overEst*(MSH%nOvlpPoints+overHead)
        call wLog ("bufferSize = ")
        call wLog (bufferSize)
        allocate(buffer(bufferSize))
        call MPI_BUFFER_ATTACH(buffer, double_size*bufferSize,code)
        call wLog ("buffer allocation code= ")
        call wLog (code)

        !Allocation
        allocate(request(requestSize))
        allocate(status(MPI_STATUS_SIZE, size(request)))

        !Allocating Temp Random Field
        tempRandField = 0 !Contributions over the random field
        countReq = 0

        !call show_MESH(MSH)
        !write(*,*) "FROM Proc rank = ", MSH%rang

        call show_MESHneigh(MSH, " ", onlyExisting = .true., forLog = .true.)

        call wLog(" RANG = ")
        call wLog(MSH%rang)


        do stage = 1, 2 !Sending and then receiving


            call wLog(" ")
            call wLog(" ")
            call wLog(" =====================================================")
            call wLog(" ===============================================")
            if(stage == 1) call wLog(" SENDING =========================================")
            if(stage == 2) call wLog(" RECEIVING =========================================")
            call wLog(" =====================================================")
            call wLog(" ===============================================")
            call wLog(" ")

            do direction = 1, size(MSH%neigh)

                if(MSH%neigh(direction) < 0) cycle !Check if this direction exists

                minPos = nint((MSH%xMinNeigh(:, direction) - MSH%xMinBound(:))/MSH%xStep) + 1
                maxPos = nint((MSH%xMaxNeigh(:, direction) - MSH%xMinBound(:))/MSH%xStep) + 1
                totalSize = (product(maxPos - minPos + 1))
                call wLog(" ")
                call wLog(" ")
                call wLog(" DIRECTION =========================================")
                call wLog(direction)
                !                call wLog(" minPos    = ")
                !                call wLog(minPos)
                !                call wLog(" maxPos    = ")
                !                call wLog(maxPos)
                !                call wLog(" totalsize = ")
                !                call wLog(totalSize)
                !                if(MSH%nDim == 2) then
                !                    call wLog("Point in minimal position = ")
                !                    call wLog(RDF%xPoints_2D(:,minPos(1),minPos(2)))
                !                    call wLog("Point in maximal position = ")
                !                    call wLog(RDF%xPoints_2D(:,maxPos(1),maxPos(2)))
                !                else if(MSH%nDim == 3) then
                !                    call wLog("Point in minimal position = ")
                !                    call wLog(RDF%xPoints_3D(:,minPos(1),minPos(2),minPos(3)))
                !                    call wLog("Point in maximal position = ")
                !                    call wLog(RDF%xPoints_3D(:,maxPos(1),maxPos(2),maxPos(3)))
                !                end if
                call wLog(" ")

                do neighPos = 1, size(MSH%neigh)
                    if(MSH%neigh(neighPos) < 0) cycle !Check if this neighbour exists
                    if(neighpos /= direction) cycle

                    !Checking if we should send this part of the field to this neighbour
                    sndRcv = .true.
                    do i = 1, MSH%nDim
                        if(       (MSH%neighShift(i, neighPos) /= 0) &
                            .and. (MSH%neighShift(i, neighPos) /= MSH%neighShift(i, direction))) then
                            sndRcv = .false.
                            exit
                        end if
                    end do

                    !if(MSH%neigh(neighPos) == 0) cycle !FOR TESTS

                    if (sndRcv) then

                        call wLog("countReq")
                        call wLog(countReq)
                        call wLog("neighPos")
                        call wLog(neighPos)

                        if(stage == 1) then

                            !if(MSH%neigh(neighPos) /= testrank) cycle !FOR TESTS
                            !if(MSH%neigh(neighPos) /= testrank2) cycle !FOR TESTS
                            !if(MSH%neigh(neighPos) /= testrank .and. MSH%neigh(neighPos) /= testrank2) cycle
                            !if(MSH%neigh(neighPos) > testrank) cycle !FOR TESTS
                            !if(MSH%neigh(neighPos) < testrank .or. MSH%neigh(neighPos) > testrank2) cycle !FOR TESTS
                            !if(neighPos /= 2) cycle !FOR TESTS
                            !if(neighPos /= 2 .and. neighPos /= 8) cycle !FOR TESTS

                            countReq = countReq + 1
                            tag = findTag(MSH, neighPos, direction, send = .true.)
                            call wLog(" --------------------")
                            call wLog(" Tag send in dir ")
                            call wLog(direction)
                            call wLog(" TAG ")
                            call wLog(tag)
                            call wLog("    TO  rang ")
                            call wLog(MSH%neigh(neighPos))
                            !end if
                            !call MPI_ISEND (RDF%randField(minPos:maxPos,1), totalSize, MPI_DOUBLE_PRECISION, &
                            !                MSH%neigh(neighPos), tag, RDF%comm, request(countReq), code)
                            if(RDF%nDim == 2) then
                                !call dispCarvalhol(RDF%RF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
                                !                    "RDF%RF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2))", &
                                !                    unit_in=RDF%log_ID)
                                !call MPI_ISEND (RDF%RF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
                                !                totalSize, MPI_DOUBLE_PRECISION, &
                                !                MSH%neigh(neighPos), tag, RDF%comm, request(countReq), code)
                                call MPI_IBSEND (RDF%RF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
                                    totalSize, MPI_DOUBLE_PRECISION, &
                                    MSH%neigh(neighPos), tag, RDF%comm, request(countReq), code)

                            else if(RDF%nDim == 3) then
                                call MPI_IBSEND (RDF%RF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)), &
                                    totalSize, MPI_DOUBLE_PRECISION, &
                                    MSH%neigh(neighPos), tag, RDF%comm, request(countReq), code)

                                !call wLog("RDF%RF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3))")
                                !RDF%RF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)) = 4
                                !call wLog("MSH%neigh(neighPos)")
                                !call wLog(MSH%neigh(neighPos))
                                !call wLog("request(countReq)")
                                !call wLog(request(countReq))
                            end if

                        else if(stage == 2) then

                            !if(neighPos /= testrank) cycle !FOR TESTS
                            !if(neighPos /= testrank2) cycle !FOR TESTS
                            !if(neighPos /= testrank .and. neighPos /= testrank2) cycle !FOR TESTS
                            !if(MSH%rang /= testrank) cycle !FOR TESTS
                            !if(MSH%rang > testrank) cycle !FOR TESTS
                            !if(MSH%rang < testrank .or. MSH%rang > testrank2) cycle !FOR TESTS

                            countReq = countReq + 1
                            tag = findTag(MSH, neighPos, direction, send = .false.)
                            call wLog(" --------------------")
                            call wLog("Tag rcv in dir ")
                            call wLog(direction)
                            call wLog(" TAG ")
                            call wLog(tag)
                            call wLog("    FROM  rang ")
                            call wLog(MSH%neigh(neighPos))
                            if(RDF%nDim == 2) then
                                call wLog("    shape(TRF_2D) ")
                                call wLog(shape(TRF_2D))
                                call wLog("    shape(RDF%RF_2D) ")
                                call wLog(shape(RDF%RF_2D))
                            else if(RDF%nDim == 3) then
                                call wLog("    shape(TRF_3D) ")
                                call wLog(shape(TRF_3D))
                                call wLog("    shape(RDF%RF_3D) ")
                                call wLog(shape(RDF%RF_3D))
                            end if

                            !call MPI_RECV (tempRandField(minPos:maxPos,1), totalSize, MPI_DOUBLE_PRECISION, &
                            !               MSH%neigh(neighPos), tag, RDF%comm, status(:,countReq), code)
                            if(RDF%nDim == 2) then
                                call MPI_RECV (TRF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
                                    totalSize, MPI_DOUBLE_PRECISION, &
                                    MSH%neigh(neighPos), tag, RDF%comm, status(:,countReq), code)
                                !RDF%randField(minPos:maxPos,1) = RDF%randField(minPos:maxPos,1) + tempRandField(minPos:maxPos,1)
                                !call dispCarvalhol(TRF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
                                !                   "TRF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2))", &
                                !                   unit_in=RDF%log_ID)
                                RDF%RF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)) = &
                                    TRF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)) &
                                    + RDF%RF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2))

                            else if(RDF%nDim == 3) then
                                call MPI_RECV (TRF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)), &
                                    totalSize, MPI_DOUBLE_PRECISION, &
                                    MSH%neigh(neighPos), tag, RDF%comm, status(:,1), code)
                                !RDF%randField(minPos:maxPos,1) = RDF%randField(minPos:maxPos,1) + tempRandField(minPos:maxPos,1)
                                RDF%RF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)) = &
                                    TRF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)) &
                                    + RDF%RF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3))
                            end if
                        end if
                    end if
                end do
            end do
        end do

        !write(*,*) "WAITING"
        !call MPI_WAITALL (countReq, request(1:countReq), status(:,1:countReq), code)
        !if(MSH%rang == testrank) write(*,*) "RDF%randField(:,1) = ", RDF%randField(:,1)

        call MPI_BUFFER_DETACH (buffer,double_size*bufferSize,code)
        if(allocated(buffer)) deallocate(buffer)

        if(allocated(request)) deallocate(request)
        if(allocated(status))  deallocate(status)

        if(associated(TRF_2D)) nullify(TRF_2D)
        if(associated(TRF_3D)) nullify(TRF_3D)


        end subroutine addNeighboursFields

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine addNeighboursFieldsV2(RDF, MSH)
        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF

        !INPUT
        type(MESH), intent(in) :: MSH

        !LOCAL
        integer :: i, direction, myRank
        integer :: op_neighRank, neighPos, neighRank
        integer, dimension(RDF%nDim) :: dirShift, neighShift, tempShift
        integer, dimension(RDF%nDim) :: op_dirShift, op_neighShift
        double precision, dimension(RDF%nDim) ::xMinDir, xMaxDir
        integer :: op_direction, op_neighPos
        integer :: code
        integer, dimension(RDF%nDim) :: minPos, maxPos
        logical, dimension(RDF%nDim) :: shiftMask
        integer :: totalSize, double_size, overHead, overEst, bufferSize
        double precision, dimension(:), allocatable :: buffer
        integer :: testrank = 1, testrank2 = 11
        double precision, dimension(1:MSH%xNTotal), target :: tempRandField, tempRandField_Now
        !integer, dimension(:), allocatable :: request
        integer :: request
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer :: tag
        logical :: sndRcv, snd, rcv
        double precision, dimension(:, :), pointer :: TRF_2D, TRF_2D_Now
        double precision, dimension(:, :, :), pointer :: TRF_3D, TRF_3D_Now


        call wLog ("Inside addNeighboursFieldsV2")


        !Buffer allocation
        call wLog ("Allocating buffer")
        !overEst = MSH%nDim + 1
        overEst = 1
        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,double_size,code)
        overHead = int(1+(MPI_BSEND_OVERHEAD*1.)/double_size)
        bufferSize = overEst*(MSH%nOvlpMax+overHead)
        call wLog ("bufferSize = ")
        call wLog (bufferSize)
        allocate(buffer(bufferSize))
        call MPI_BUFFER_ATTACH(buffer, double_size*bufferSize,code)
        call wLog ("buffer allocation code= ")
        call wLog (code)

        !Allocation
        if(RDF%nDim == 2) TRF_2D(1:MSH%xNStep(1),1:MSH%xNStep(2)) => tempRandField
        if(RDF%nDim == 2) TRF_2D_Now(1:MSH%xNStep(1),1:MSH%xNStep(2)) => tempRandField_Now
        if(RDF%nDim == 3) TRF_3D(1:MSH%xNStep(1),1:MSH%xNStep(2),1:MSH%xNStep(3)) => tempRandField
        if(RDF%nDim == 3) TRF_3D_Now(1:MSH%xNStep(1),1:MSH%xNStep(2),1:MSH%xNStep(3)) => tempRandField_Now

        if(RDF%nDim == 2) then
            call wLog("    shape(TRF_2D) ")
            call wLog(shape(TRF_2D))
            call wLog("    shape(RDF%RF_2D) ")
            call wLog(shape(RDF%RF_2D))
        else if(RDF%nDim == 3) then
            call wLog("    shape(TRF_3D) ")
            call wLog(shape(TRF_3D))
            call wLog("    shape(RDF%RF_3D) ")
            call wLog(shape(RDF%RF_3D))
        end if

        !Initialization
        tempRandField = 0 !Contributions over the random field
        tag = 0

        !call show_MESH(MSH)
        !write(*,*) "FROM Proc rank = ", MSH%rang

        !call show_MESHneigh(MSH, " ", onlyExisting = .true., forLog = .true.)
        call show_MESHneigh(MSH, " ", onlyExisting = .false., forLog = .true.)

        call wLog(" RANG = ")
        call wLog(MSH%rang)
        myRank = MSH%rang

        !call DispCarvalhol(MSH%mappingFromShift, "MSH%mappingFromShift", unit_in=MSH%log_ID)
        !call DispCarvalhol(MSH%neighShift, "MSH%neighShift", unit_in=MSH%log_ID)

        call wLog(" ")

        do direction = 1, size(MSH%neigh)

            call MPI_BARRIER(MSH%comm, code) !This is necessary so we can reduce the size of the buffer

            call wLog(" ")
            call wLog(" ")
            call wLog(" DIRECTION =========================================")
            call wLog(direction)
            call wLog("    Shift")
            call wLog(MSH%neighShift(:, direction))

            !SENDING---------------------------------------------------------------

            snd = .true.
            dirShift   = MSH%neighShift(:, direction)
            neighShift = dirShift
            neighRank  = MSH%neigh(direction)

            if(MSH%neigh(direction) < 0) snd = .false. !Check if this direction exists

            if(snd) then
                call wLog(" ")
                call wLog(" SENDING ==============")
                !Defining slice that should be sent
                !call wLog("DIMENSIONING GOING SLICE!!!")

                xMinDir = MSH%xMinNeigh(:, direction)
                xMaxDir = MSH%xMaxNeigh(:, direction)
                minPos = nint((xMinDir - MSH%xMinBound(:))/MSH%xStep) + 1
                maxPos = nint((xMaxDir - MSH%xMinBound(:))/MSH%xStep) + 1
                totalSize = (product(maxPos - minPos + 1))
                tag       = neighRank

                call wLog("totalSize")
                call wLog(totalSize)
                call wLog("minPos")
                call wLog(minPos)
                call wLog("maxPos")
                call wLog(maxPos)
                call wLog("dirShift")
                call wLog(dirShift)
                call wLog(" ")
                call wLog(" TAG ")
                call wLog(tag)
                call wLog("    TO  rang ")
                call wLog(neighRank)

                if(MSH%nDim == 2) then
                    call wLog("Point in minimal position = ")
                    call wLog(RDF%xPoints_2D(:,minPos(1),minPos(2)))
                    call wLog("Point in maximal position = ")
                    call wLog(RDF%xPoints_2D(:,maxPos(1),maxPos(2)))
    !                call MPI_ISEND (RDF%RF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
    !                    totalSize, MPI_DOUBLE_PRECISION, &
    !                    neighRank, tag, RDF%comm, request, code)

                    call MPI_IBSEND (RDF%RF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
                            totalSize, MPI_DOUBLE_PRECISION, &
                            neighRank, tag, RDF%comm, request, code)

                end if
                if(MSH%nDim == 3) then
                    call wLog("Point in minimal position = ")
                    call wLog(RDF%xPoints_3D(:,minPos(1),minPos(2),minPos(3)))
                    call wLog("Point in maximal position = ")
                    call wLog(RDF%xPoints_3D(:,maxPos(1),maxPos(2),maxPos(3)))
    !                call MPI_ISEND (RDF%RF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)), &
    !                    totalSize, MPI_DOUBLE_PRECISION, &
    !                    neighRank, tag, RDF%comm, request, code)

                    call MPI_IBSEND (RDF%RF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)), &
                        totalSize, MPI_DOUBLE_PRECISION, &
                        neighRank, tag, RDF%comm, request, code)
                end if
            else

                call wLog(" not sending ===========")

            end if

            !RECEIVING---------------------------------------------------------------
            rcv = .true.
            dirShift   = -dirShift
            neighShift = -neighShift
            neighRank  = MSH%op_neigh(direction)

            if(MSH%op_neigh(direction) < 0) rcv = .false. !Check if this direction exists

            !Finding opposite direction index
            do op_direction = 1, size(MSH%neigh)
                if(all(dirShift == MSH%neighShift(:, op_direction))) exit
            end do

            call wLog("op_direction = ")
            call wLog(op_direction)

            if(rcv) then
                call wLog(" ")
                call wLog(" RECEIVING ==============")
                !Defining slice that should be sent

                xMinDir = MSH%xMinNeigh(:, op_direction)
                xMaxDir = MSH%xMaxNeigh(:, op_direction)

                minPos = nint((xMinDir - MSH%xMinBound(:))/MSH%xStep) + 1
                maxPos = nint((xMaxDir - MSH%xMinBound(:))/MSH%xStep) + 1
                totalSize = (product(maxPos - minPos + 1))
                neighRank = MSH%op_neigh(direction)
                tag       = myRank

                call wLog("totalSize")
                call wLog(totalSize)
                call wLog("minPos")
                call wLog(minPos)
                call wLog("maxPos")
                call wLog(maxPos)
                call wLog("dirShift")
                call wLog(dirShift)
                call wLog(" ")
                call wLog(" TAG ")
                call wLog(tag)
                call wLog("    FROM  rang ")
                call wLog(neighRank)

                tempRandField_Now = 0.0D0

                if(MSH%nDim == 2) then
                    call wLog("Point in minimal position = ")
                    call wLog(RDF%xPoints_2D(:,minPos(1),minPos(2)))
                    call wLog("Point in maximal position = ")
                    call wLog(RDF%xPoints_2D(:,maxPos(1),maxPos(2)))
                    call MPI_RECV (TRF_2D_Now(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
                        totalSize, MPI_DOUBLE_PRECISION, &
                        neighRank, tag, RDF%comm, status, code)
                        TRF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)) = &
                        TRF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2))   &
                        + TRF_2D_Now(minPos(1):maxPos(1),minPos(2):maxPos(2))
                end if
                if(MSH%nDim == 3) then
                    call wLog("Point in minimal position = ")
                    call wLog(RDF%xPoints_3D(:,minPos(1),minPos(2),minPos(3)))
                    call wLog("Point in maximal position = ")
                    call wLog(RDF%xPoints_3D(:,maxPos(1),maxPos(2),maxPos(3)))
                    call MPI_RECV (TRF_3D_Now(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)), &
                        totalSize, MPI_DOUBLE_PRECISION, &
                        neighRank, tag, RDF%comm, status, code)
                        TRF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)) = &
                        TRF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3))   &
                        + TRF_3D_Now(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3))
                end if
            else
                call wLog(" not receiving ==============")
            end if
        end do

        RDF%randField(:,1) = RDF%randField(:,1) + tempRandField

        call MPI_BUFFER_DETACH (buffer,double_size*bufferSize,code)
        if(allocated(buffer)) deallocate(buffer)

        if(associated(TRF_2D_Now)) nullify(TRF_2D_Now)
        if(associated(TRF_3D_Now)) nullify(TRF_3D_Now)
        if(associated(TRF_2D)) nullify(TRF_2D)
        if(associated(TRF_3D)) nullify(TRF_3D)


    end subroutine addNeighboursFieldsV2


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine applyWeightingFunctions_OnMatrix(RDF, MSH, partitionType)
    !INPUT
    type(RF), intent(in) :: RDF
    type(MESH), intent(in) :: MSH
    integer, intent(in) :: partitionType

    !LOCAL
    integer :: direction, i
    double precision, dimension(1:MSH%xNTotal) :: unityPartition
    double precision, dimension(MSH%nDim) :: originCorner
    integer, dimension(MSH%nDim) :: minPos, maxPos

    !call wLog("shape(RDF%RF_2D) = ")
    !call wLog(shape(RDF%RF_2D))
    !call wLog("shape(RDF%xPoints_2D) = ")
    !call wLog(shape(RDF%xPoints_2D))
    !call wLog("RDF%xPoints_2D(:,1,1) = ")
    !call wLog(RDF%xPoints_2D(:,1,1))
    !call wLog("RDF%xPoints_2D(:,30,29) = ")
    !call wLog(RDF%xPoints_2D(:,30,29))

    call wLog("Inside 'applyWeightingFunctions_OnMatrix'")


    !Modify extremes of local Random Field-------------------------------------------------------
    unityPartition(:) = 1

    !Building Shape Functions in all directions
    do direction = 1, size(MSH%neigh)

        !if(.not. MSH%considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction
        if(MSH%neigh(direction) < 0) cycle
        if(sum(abs(MSH%neighShift(:, direction))) /= 1) cycle !We need only the "pure directions"

        !if(direction/=4) cycle !FOR TESTS

        !Finding range
        originCorner = MSH%xOrNeigh(:, direction)
        minPos = nint((MSH%xMinNeigh(:, direction) - MSH%xMinBound(:))/MSH%xStep) + 1
        maxPos = nint((MSH%xMaxNeigh(:, direction) - MSH%xMinBound(:))/MSH%xStep) + 1

        call wLog("DIRECTION")
        call wLog(direction)
        call wLog("originCorner")
        call wLog(originCorner)
        call wLog("MSH%neighShift(:, direction)")
        call wLog(MSH%neighShift(:, direction))
        call wLog("minPos = ")
        call wLog(minPos)
        call wLog("maxPos")
        call wLog(maxPos)

        !Shape Function Generation
        call generateUnityPartition_OnMatrix(RDF, MSH%xNStep, originCorner, MSH%overlap, &
            MSH%neighShift(:, direction), partitionType, &
            unityPartition, minPos, maxPos)

    end do !Direction

    do i = 1, size(RDF%randField(:,1),1)
        RDF%randField(i,1) = RDF%randField(i,1) * sqrt(unityPartition(i))
    end do
    !RDF%randField(:,1) = RDF%randField(:,1) * sqrt(unityPartition(:)) !This syntax demands too much memory for Occygen

    !RDF%randField(:,1) = unityPartition(:) !TEST
    !RDF%randField(:,1) = sqrt(unityPartition(:)) !TEST

end subroutine applyWeightingFunctions_OnMatrix

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
subroutine generateUnityPartition_OnMatrix(RDF, xNStep, originCorner, overlap, neighShift, &
    partitionType, unityPartition, minPos, maxPos)

    implicit none

    !INPUT
    type(RF), intent(in) ::  RDF
    integer, dimension(:), intent(in) ::  xNStep
    double precision, dimension(:)  , intent(in) :: originCorner, overlap
    integer, dimension(:)  , intent(in) :: neighShift
    integer, intent(in) :: partitionType
    integer, dimension(:), intent(in) :: minPos, maxPos

    !OUTPUT
    double precision, dimension(:), intent(out), target :: unityPartition

    !LOCAL
    integer :: i
    double precision, dimension(:, :), pointer :: UP_2D
    double precision, dimension(:, :, :), pointer :: UP_3D

    if(RDF%nDim == 2) UP_2D(1:xNStep(1),1:xNStep(2)) => unityPartition
    if(RDF%nDim == 3) UP_3D(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => unityPartition

    !call wLog("shape(unityPartition) = ")
    !call wLog(shape(unityPartition))
    !call wLog("shape(UP_2D) = ")
    !call wLog(shape(UP_2D))
    !call wLog("shape(RDF%xPoints_2D) = ")
    !call wLog(shape(RDF%xPoints_2D))
    !call wLog("shape(UP_3D) = ")
    !call wLog(shape(UP_3D))
    !call wLog("shape(RDF%xPoints_3D) = ")
    !call wLog(shape(RDF%xPoints_3D))

    !call wLog("Point in minimal position = ")
    !call wLog(RDF%xPoints_3D(:,22,1,1))
    !call wLog(RDF%xPoints_3D(:,minPos(1),minPos(2),minPos(3)))



    do i = 1, RDF%nDim
        !do i = 1, 1 !For Tests
        if(neighShift(i) == 0) cycle
        !            call wLog("Dimension = ")
        !call wLog(i)

        if(partitionType == 1) then
            if(RDF%nDim == 2) then
                !                    call wLog("Point in minimal position = ")
                !                    call wLog(RDF%xPoints_2D(i,minPos(1),minPos(2)))
                !                    call wLog("Point in maximal position = ")
                !                    call wLog(RDF%xPoints_2D(i,maxPos(1),maxPos(2)))
                !                    call wLog("originCorner(i) = ")
                !                    call wLog(originCorner(i))
                !UP_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)) = 1 !For Tests
                UP_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)) = &
                    ((1.0D0 + cos(PI*(RDF%xPoints_2D(&
                    i,minPos(1):maxPos(1),minPos(2):maxPos(2)) &
                    - originCorner(i))/overlap(i)))&
                    / 2.0D0) &
                    * UP_2D(minPos(1):maxPos(1),minPos(2):maxPos(2))
            else if (RDF%nDim == 3) then
                !                    call wLog("Point in minimal position = ")
                !                    call wLog(RDF%xPoints_3D(i,minPos(1),minPos(2),minPos(3)))
                !                    call wLog("Point in maximal position = ")
                !                    call wLog(RDF%xPoints_3D(i,maxPos(1),maxPos(2),maxPos(3)))
                !                    call wLog("originCorner(i) = ")
                !                    call wLog(originCorner(i))
                UP_3D(minPos(1):maxPos(1),minPos(2):maxPos(2), minPos(3):maxPos(3)) = &
                    ((1.0D0 + cos(PI*(RDF%xPoints_3D(&
                    i,minPos(1):maxPos(1),minPos(2):maxPos(2), minPos(3):maxPos(3)) &
                    - originCorner(i))/overlap(i)))&
                    / 2.0D0) &
                    * UP_3D(minPos(1):maxPos(1),minPos(2):maxPos(2), minPos(3):maxPos(3))
            else
                stop("Unity Partition not implemented in this Dimension (generateUnityPartition_OnMatrix)")
            end if
        else
        stop('ERROR!! Inside "generateUnityPartition_OnMatrix" - partition Type not defined')
    end if
end do

if(associated(UP_2D)) nullify(UP_2D)
if(associated(UP_3D)) nullify(UP_3D)

end subroutine generateUnityPartition_OnMatrix

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
subroutine generateUnityPartition(xPoints, originCorner, overlap, neighShift, partitionType, unityPartition)

    implicit none

    !INPUT
    double precision, dimension(:,:), intent(in) ::  xPoints
    double precision, dimension(:)  , intent(in) :: originCorner, overlap
    integer, dimension(:)  , intent(in) :: neighShift
    integer, intent(in) :: partitionType

    !OUTPUT
    double precision, dimension(:), intent(out) :: unityPartition

    !LOCAL
    integer :: i, nDim

    nDim = size(originCorner)


    unityPartition = 1.0D0
    do i = 1, nDim
        if(neighShift(i) == 0) cycle

        if(partitionType == 1) then
            unityPartition = ((1.0D0 + cos(PI*(xPoints(i,:) - originCorner(i))/overlap(i)))&
                / 2.0D0) &
                * unityPartition
        else
            stop('ERROR!! Inside "generateUnityPartition" - partition Type not defined')
        end if
    end do

end subroutine generateUnityPartition

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

    orig  = (dble(MSH%origin-1)*MSH%xStep + MSH%xMinGlob)
    xStep = MSH%xStep

    call wLog("xStep")
    call wLog(xStep)
    call wLog("orig")
    call wLog(orig)
    call wLog("MSH%xNStep")
    call wLog(MSH%xNStep)

    offset(1) = 1
    do i = 2, MSH%nDim
        offset(i) = product(MSH%xNStep(1:i-1))
    end do

    call wLog("offset")
    call wLog(offset)
    call wLog("shape(RDF%randField)")
    call wLog(shape(RDF%randField))

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
    integer :: i


    if(MSH%neigh(neighPos) < 0 .or. MSH%neigh(direction) < 0) then
        !write(*,*) "Inside findTag , Invalid Neighbour"
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

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine copy_RF_xPoints(MSH, orRDF, destRDF, destXPoints, direction)
!        implicit none
!
!        !INPUT AND OUTPUT
!        type(MESH), intent(in) :: MSH
!        type(RF), intent(in)   :: orRDF
!        type(RF) ::destRDF
!        integer, intent(in) :: direction
!        double precision, dimension(:, :), allocatable, intent(out), target :: destXPoints;
!
!        !LOCAL
!        integer :: i, j, totalSize
!
!        totalSize = MSH%indexNeigh(2,direction) - MSH%indexNeigh(1,direction) + 1
!
!        if(allocated(destXPoints)) then
!            if(size(destXPoints) /= totalSize) then
!                nullify(destRDF%xPoints)
!                deallocate(destXPoints)
!            end if
!        end if
!
!        if(.not.allocated(destXPoints)) allocate(destXPoints(MSH%nDim, 1:totalSize))
!
!        destRDF%xPoints => destXPoints
!
!        destRDF%xPoints(:,:) = orRDF%xPoints(:, MSH%indexNeigh(1,direction):MSH%indexNeigh(2,direction))
!        destRDF%xNTotal = totalSize
!
!    end subroutine copy_RF_xPoints

end module localization_RF
