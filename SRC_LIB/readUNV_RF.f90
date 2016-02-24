module readUNV_RF
    use displayCarvalhol
    use write_Log_File

    implicit none

    integer :: lineNb

contains
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine readUNV(path, nDim, coordList, connectList, monotype, rang, nb_procs)
        !INPUT
        character (len=*), intent(in) :: path
        integer, intent(in)           :: nDim, rang, nb_procs

        !OUTPUT
        double precision, dimension(:,:), allocatable, intent(out), optional :: coordList
        integer         , dimension(:,:), allocatable, intent(out), optional :: connectList
        logical, optional, intent(out) :: monoType

        !LOCAL
        integer :: fileID = 600
        character (len=200) :: line
        integer :: stat
        logical :: inside
        integer :: nNodes
        integer :: nNodesLoc, nodeStart, nodeEnd
        integer :: nElemLoc, elemStart, elemEnd
        integer :: nElem, maxConnect
        !integer, allocatable, dimension(:) :: elemSizes
        integer, dimension(:)  , allocatable :: sizeList

        !integer, dimension(:,:), allocatable :: startEnd

        stat     = 0;
        lineNb = 0
        inside = .false.

        open (unit = fileID , file = path, action = 'read')

            !Dimensioning lecture
            call wLog("--------------------------------------------------")
            call wLog("DIMENSIONING LECTURE------------------------------")
            call wLog("--------------------------------------------------")
            call wLog("")

            !allocate(startEnd(2, nb_procs))

            read(fileID, fmt=*, IOSTAT = stat) line
            lineNb = lineNb + 1

            do while (stat == 0)

                if(trim(adjustL(line)) == "-1") then
                    !!!write(get_fileId(),*) "Line ", lineNb, "is -1"
                else if((trim(adjustL(line)) == "2411" .or. trim(adjustL(line)) == "781"))then
                    !!!write(get_fileId(),*) " "
                    !!!write(get_fileId(),*) "Dimensioning Coordinates"
                    !!!write(get_fileId(),*) "Start Line = ", lineNb
                    call prepareCoordinates(fileID, nNodes)
                else if((trim(adjustL(line)) == "2412" .or. trim(adjustL(line)) == "780")) then
                    !!!write(get_fileId(),*) " "
                    !!!write(get_fileId(),*) "Dimensioning Connectivity"
                    !!!write(get_fileId(),*) "Start Line = ", lineNb
                    call prepareConnectivity(fileID, nElem, maxConnect, monoType)
                !else if(trim(adjustL(line)) == "2477") then
                    !write(*,*) " "
                    !write(*,*) "Dimensioning Physical Volumes (NOT YET)"
                    !write(*,*) "Start Line = ", lineNb
                end if

                read(fileID, fmt=*, IOSTAT = stat) line
                lineNb = lineNb + 1

            end do

            !Sharing nodes between processors
            if((nNodes < nb_procs .and. present(coordList)) .or. &
                nElem < nb_procs .and. present(connectList)) then
                !!write(get_fileId(),*) "ERROR!!! Too little Nodes/Elements for this number of processors"
                !!write(get_fileId(),*) " nNodes   = ", nNodes
                !!write(get_fileId(),*) " nb_procs = ", nb_procs
                !!write(get_fileId(),*) " nElem    = ", nElem
                stop
            end if

            nodeStart = (rang*(nNodes/nb_procs)) + 1
            nodeEnd   = ((rang+1)*(nNodes/nb_procs))
            if((rang+1) == nb_procs) nodeEnd = nNodes
            nNodesLoc = nodeEnd - nodeStart + 1

            elemStart = (rang*(nElem/nb_procs)) + 1
            elemEnd   = ((rang+1)*(nElem/nb_procs))
            if((rang+1) == nb_procs) elemEnd = nElem
            nElemLoc = elemEnd - elemStart + 1

            !!write(get_fileId(),*) "nodeStart = ", nodeStart
            !!write(get_fileId(),*) "nodeEnd   = ", nodeEnd
            !!write(get_fileId(),*) "elemStart = ", elemStart
            !!write(get_fileId(),*) "elemEnd   = ", elemEnd

            !Allocation
            if(present(coordList)) then
                allocate(coordList(nDim, nNodesLoc))
                coordList(:,:)   = -1
                !!write(get_fileId(),*) "shape(coordList) = ", shape(coordList)
            end if
            if(present(connectList)) then
                allocate(connectList(1:maxConnect, nElemLoc))
                allocate(sizeList(nElemLoc))
                connectList(:,:) = -1
                !!write(get_fileId(),*) "shape(connectList) = ", shape(connectList)
            end if

!
!            end if

            !Real Lecture
            !!write(get_fileId(),*) "--------------------------------------------------"
            !!write(get_fileId(),*) "EFFECTIVE LECTURE------------------------------"
            !!write(get_fileId(),*) "--------------------------------------------------"
            !!write(get_fileId(),*) ""
            rewind(fileID)
            lineNb = 0

            read(fileID, fmt=*, IOSTAT = stat) line
            lineNb = lineNb + 1

            do while (stat == 0)

                if(trim(adjustL(line)) == "-1") then
                    !!!write(get_fileId(),*) "Line ", lineNb, "is -1"
                else if((trim(adjustL(line)) == "2411" .or. trim(adjustL(line)) == "781") &
                        .and. present(coordList) )then
                    !!!write(get_fileId(),*) " "
                    !!!write(get_fileId(),*) "Line = ", lineNb
                    !!write(get_fileId(),*) "Reading Coordinates"
                    call readCoordinates(fileID, coordList, nodeStart, nodeEnd)
                else if((trim(adjustL(line)) == "2412" .or. trim(adjustL(line)) == "780") &
                        .and. present(connectList))then
                    !!!write(get_fileId(),*) " "
                    !!!write(get_fileId(),*) "Line = ", lineNb
                    !!write(get_fileId(),*) "Reading Connectivity"
                    call readConnectivity(fileID, connectList, sizeList, elemStart, elemEnd, monoType)

                !else if(trim(adjustL(line)) == "2477") then
                    !write(*,*) " "
                    !write(*,*) "Line = ", lineNb
                    !write(*,*) "Reading Physical Volumes (NOT YET)"
                end if


                read(fileID, fmt=*, IOSTAT = stat) line
                lineNb = lineNb + 1

            end do

        close(fileID)

        if(allocated(sizeList)) deallocate(sizeList)

        !!write(get_fileId(),*) "--------------------------------------------------"
        !!write(get_fileId(),*) "END OF UNV LECTURE------------------------------"
        !!write(get_fileId(),*) "--------------------------------------------------"
        !!write(get_fileId(),*) ""

        !if(present(coordList)) call dispCarvalhol(transpose(coordList(:, :10)), "coord List", "(F15.5)")
        !if(present(connectList)) call dispCarvalhol(transpose(connectList(:, :10)), "connect List", "(I8)", 8)
        !if(present(sizeList)) call dispCarvalhol(sizeList(:10), "size List", "(I8)")

    end subroutine readUNV

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine readUNV_many(pathList, nDim, coordList, connectList, monotype, rang, nb_procs)
        !INPUT
        character (len=*), dimension(:), intent(in) :: pathList
        integer, intent(in)           :: nDim, rang, nb_procs

        !OUTPUT
        double precision, dimension(:,:), allocatable, intent(out), optional :: coordList
        integer         , dimension(:,:), allocatable, intent(out), optional :: connectList
        logical, optional, intent(out) :: monoType

        !LOCAL
        integer :: fileID = 600
        character (len=200) :: line
        integer :: stat
        integer :: f
        logical :: inside
        integer :: nNodes
        integer :: nNodesLoc, nodeStart, nodeEnd
        integer :: nElemLoc, elemStart, elemEnd
        integer :: nElem, maxConnect
        !integer, allocatable, dimension(:) :: elemSizes
        integer, dimension(:)  , allocatable :: sizeList
        integer :: nFiles
        integer, dimension(size(pathList)) :: nElemFile, nNodesFile
        !integer, dimension(:,:), allocatable :: startEnd
        logical :: readNodes, readElem
        integer :: startPointNodes, endPointNodes, startPointElem, endPointElem
        integer :: nodeMaxFile, nodeMinFile, elemMaxFile, elemMinFile
        integer :: nodeOffset, elemOffset

        nFiles = size(pathList)

        stat       = 0
        lineNb     = 0
        inside     = .false.
        monoType   = .true.
        maxConnect = 0

        !Dimensioning lecture
        call wLog("--------------------------------------------------")
        call wLog("DIMENSIONING LECTURE------------------------------")
        call wLog("--------------------------------------------------")
        call wLog("")

        !allocate(startEnd(2, nb_procs))

        do f = 1, nFiles

            open (unit = fileID , file = pathList(f), action = 'read')

            read(fileID, fmt=*, IOSTAT = stat) line
            lineNb = lineNb + 1

            do while (stat == 0)

                if(trim(adjustL(line)) == "-1") then
                    !!!write(get_fileId(),*) "Line ", lineNb, "is -1"
                else if((trim(adjustL(line)) == "2411" .or. trim(adjustL(line)) == "781"))then
                    !!!write(get_fileId(),*) " "
                    !!!write(get_fileId(),*) "Dimensioning Coordinates"
                    !!!write(get_fileId(),*) "Start Line = ", lineNb
                    call prepareCoordinates(fileID, nNodes)
                else if((trim(adjustL(line)) == "2412" .or. trim(adjustL(line)) == "780")) then
                    !!!write(get_fileId(),*) " "
                    !!!write(get_fileId(),*) "Dimensioning Connectivity"
                    !!!write(get_fileId(),*) "Start Line = ", lineNb
                    call prepareConnectivity(fileID, nElem, maxConnect, monoType)
                !else if(trim(adjustL(line)) == "2477") then
                    !write(*,*) " "
                    !write(*,*) "Dimensioning Physical Volumes (NOT YET)"
                    !write(*,*) "Start Line = ", lineNb
                end if

                read(fileID, fmt=*, IOSTAT = stat) line
                lineNb = lineNb + 1

            end do

            nElemFile(f)  = nElem
            nNodesFile(f) = nNodes

            close(fileID)

        end do

        !Sharing nodes and Elements between processors
        nElem  = sum(nElemFile)
        nNodes = sum(nNodesFile)

        call wLog(" nElemFile   = ")
        call wLog(nElemFile)
        call wLog(" nElem    = ")
        call wLog(nElem)
        call wLog(" nNodesFile    = ")
        call wLog(nNodesFile)
        call wLog(" nNodes   = ")
        call wLog(nNodes)
        call wLog(" nb_procs = ")
        call wLog(nb_procs)
        call wLog(" monoType = ")
        call wLog(monoType)
        call wLog(" maxConnect = ")
        call wLog(maxConnect)

        if((nNodes < nb_procs .and. present(coordList)) .or. &
            nElem < nb_procs .and. present(connectList)) then
            write(*,*) "ERROR!!! Too little Nodes/Elements for this number of processors"

            stop
        end if

        nodeStart = (rang*(nNodes/nb_procs)) + 1
        nodeEnd   = ((rang+1)*(nNodes/nb_procs))
        if((rang+1) == nb_procs) nodeEnd = nNodes
        nNodesLoc = nodeEnd - nodeStart + 1

        elemStart = (rang*(nElem/nb_procs)) + 1
        elemEnd   = ((rang+1)*(nElem/nb_procs))
        if((rang+1) == nb_procs) elemEnd = nElem
        nElemLoc = elemEnd - elemStart + 1

        call wLog("nodeStart = ")
        call wLog(nodeStart)
        call wLog("nodeEnd   = ")
        call wLog(nodeEnd)
        call wLog("elemStart = ")
        call wLog(elemStart)
        call wLog("elemEnd   = ")
        call wLog(elemEnd)

        !Allocation
        if(present(coordList)) then
            allocate(coordList(nDim, nNodesLoc))
            coordList(:,:)   = -1
            !!write(get_fileId(),*) "shape(coordList) = ", shape(coordList)
        end if
        if(present(connectList)) then
            allocate(connectList(1:maxConnect, nElemLoc))
            allocate(sizeList(nElemLoc))
            connectList(:,:) = -1
            !!write(get_fileId(),*) "shape(connectList) = ", shape(connectList)
        end if

        !Real Lecture
        call wLog("--------------------------------------------------")
        call wLog("EFFECTIVE LECTURE------------------------------")
        call wLog("--------------------------------------------------")
        !rewind(fileID)

        nodeOffset = 0
        elemOffset = 0

        do f = 1, nFiles

            call wLog("FILE =  ---------------------")
            call wLog(f)

            !Nodes
            if(f == 1) then
                startPointNodes = 1
            else
                startPointNodes = sum(nNodesFile(1:f-1)) + 1
            end if
            endPointNodes   = sum(nNodesFile(1:f))

           call wLog("startPointNodes = ")
           call wLog(startPointNodes)
           call wLog("endPointNodes = ")
           call wLog(endPointNodes)

            if(                                                           &
               (nodeStart <= endPointNodes .and. nodeStart >= startPointNodes) .or. &
               (nodeEnd <= endPointNodes .and. nodeEnd >= startPointNodes)          &
               ) then
               readNodes = .true.

               nodeMinFile = maxval([nodeStart, startPointNodes]) - startPointNodes + 1
               nodeMaxFile = minval([nodeEnd, endPointNodes]) - startPointNodes + 1

               call wLog("nodeMinFile = ")
               call wLog(nodeMinFile)
               call wLog("nodeMaxFile = ")
               call wLog(nodeMaxFile)

               call wLog("readNodes")
               call wLog(readNodes)
            else
                call wLog("NO NODES OF INTEREST IN THIS FILE")
                readNodes = .false.
            end if

            !Elements
            if(f == 1) then
                startPointElem = 1
            else
                startPointElem = sum(nElemFile(1:f-1)) + 1
            end if
            endPointElem   = sum(nElemFile(1:f))

            call wLog("startPointElem = ")
            call wLog(startPointElem)
            call wLog("endPointElem = ")
            call wLog(endPointElem)

            if(                                                           &
               (elemStart <= endPointElem .and. elemStart >= startPointElem) .or. &
               (elemEnd <= endPointElem .and. elemEnd >= startPointElem)          &
               ) then
               readElem = .true.

               elemMinFile = maxval([elemStart, startPointElem]) - startPointElem + 1
               elemMaxFile = minval([elemEnd, endPointElem]) - startPointElem + 1

               call wLog("elemMinFile = ")
               call wLog(elemMinFile)
               call wLog("elemMaxFile = ")
               call wLog(elemMaxFile)

               call wLog("readElem")
               call wLog(readElem)

            else
                call wLog("NO ELEMENTS OF INTEREST IN THIS FILE")
                readElem = .false.
            end if

            if((.not. readNodes) .and. (.not. readElem)) then
                call wLog("Nothing to read in this file by this processor")
                cycle
            end if

            open (unit = fileID , file = pathList(f), action = 'read')

            lineNb = 0
            stat   = 0

            !First Lecture
            read(fileID, fmt=*, IOSTAT = stat) line
            lineNb = lineNb + 1

            do while (stat == 0)

                call wLog("lineNb = ")
                call wLog(lineNb)

                if(trim(adjustL(line)) == "-1") then
                    !!!write(get_fileId(),*) "Line ", lineNb, "is -1"
                else if((trim(adjustL(line)) == "2411" .or. trim(adjustL(line)) == "781") &
                        .and. present(coordList) )then
                    !!!write(get_fileId(),*) " "
                    !!!write(get_fileId(),*) "Line = ", lineNb
                    call wLog("elemMaxFile = ")
                    call wLog("Reading Coordinates")
                    if(readNodes) call readCoordinates(fileID, coordList, nodeMinFile, nodeMaxFile, &
                                                       nodeOffset)
                    call wLog("First Coords = ")
                    call wLog(coordList(:, 1:2))
                    !call dispCarvalhol(transpose(coordList(:, 1:10)), "coord List", "(F15.5)")
                else if((trim(adjustL(line)) == "2412" .or. trim(adjustL(line)) == "780") &
                        .and. present(connectList))then
                    !!!write(get_fileId(),*) " "
                    !!!write(get_fileId(),*) "Line = ", lineNb
                    !!write(get_fileId(),*) "Reading Connectivity"
                    if(readElem) call readConnectivity(fileID, connectList, sizeList, elemMinFile, elemMaxFile, monoType, &
                                                       elemOffset, startPointNodes - 1)

                !else if(trim(adjustL(line)) == "2477") then
                    !write(*,*) " "
                    !write(*,*) "Line = ", lineNb
                    !write(*,*) "Reading Physical Volumes (NOT YET)"
                end if

                !Reading Next Line
                read(fileID, fmt=*, IOSTAT = stat) line
                lineNb = lineNb + 1

            end do

            close(fileID)
        end do

        if(allocated(sizeList)) deallocate(sizeList)

        call wLog("--------------------------------------------------")
        call wLog("END OF UNV LECTURE--------------------------------")
        call wLog("--------------------------------------------------")
        call wLog(" ")

        !if(present(coordList)) call dispCarvalhol(transpose(coordList(:, :10)), "coord List", "(F15.5)")
        !if(present(connectList)) call dispCarvalhol(transpose(connectList(:, :10)), "connect List", "(I8)", 8)
        !if(present(sizeList)) call dispCarvalhol(sizeList(:10), "size List", "(I8)")

    end subroutine readUNV_many

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine prepareCoordinates(fileID, nNodes)
        !INPUT
        integer, intent(in) :: fileID
        !OUTPUT
        integer, intent(out) :: nNodes
        !LOCAL
        character (len=200) :: line
        logical :: inside

        inside = .true.
        nNodes = 0

        !Jump initials -1
        read(fileID, fmt=*) line
        do while(trim(adjustL(line)) == "-1")
            read(fileID, fmt=*) line
            lineNb = lineNb + 1
        end do
        backspace(fileID)
        lineNb = lineNb - 1

        !Reading coordinates
        do while(inside)

            !Read coordinates header
            read(fileID, fmt=*) line
            lineNb = lineNb + 1

            if(trim(adjustL(line)) == "-1") then
                !!!write(get_fileId(),*) "Line ", lineNb, "is -1 (inside coord)"
                inside = .false.
            else
                !Counting number of nodes
                nNodes = nNodes + 1
                read(fileID, fmt=*) line
                lineNb = lineNb + 1
            end if

        end do

        !!write(get_fileId(),*) "nNodes = ", nNodes
        !!!write(get_fileId(),*) "exit Line = ", lineNb

    end subroutine prepareCoordinates

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine prepareConnectivity(fileID, nElem, maxConnect, monoType)
        !INPUT
        integer, intent(in) :: fileID
        !OUTPUT
        integer, intent(out) :: nElem
        integer, intent(inout) :: maxConnect
        logical, intent(inout) :: monoType
        !LOCAL
        integer, dimension(50) :: nElemByType
        character (len=20) :: line
        integer, dimension(6) :: connectInfo
        integer :: i
        logical :: inside

        inside = .true.
        nElem = 0
        nElemByType(:) = 0

        !Jump initials -1
        read(fileID, fmt=*) line
        do while(trim(adjustL(line)) == "-1")
            read(fileID, fmt=*) line
            lineNb = lineNb + 1
        end do
        backspace(fileID)
        lineNb = lineNb - 1

        do while(inside)

            read(fileID, fmt=*) line
            lineNb = lineNb + 1

            !Read the connectivity header information
            if(trim(adjustL(line)) == "-1") then
                !!!write(get_fileId(),*) "Line ", lineNb, "is -1 (inside connect)"
                inside = .false.
                cycle
            else
                backspace(fileID)
                lineNb = lineNb - 1
                !read(fileID, fmt=*) strVec_6
                read(fileID, *) connectInfo
                lineNb = lineNb + 1
                !write(*,*) "connectInfo = ", connectInfo
                nElem = nElem + 1
                !if(nElemByType(connectInfo(6))==0) !!write(get_fileId(),*) "First triangle is element ", nElem
                nElemByType(connectInfo(6)) = nElemByType(connectInfo(6)) + 1
            end if

            !Process max connectivity
            if(maxConnect == 0) then
                maxConnect = connectInfo(6)
                call wLog("Initial maxConnect = ")
                call wLog(maxConnect)
            end if

            if(maxConnect /= connectInfo(6)) then
                if(maxConnect < connectInfo(6)) maxConnect = connectInfo(6)
                monoType = .false.
                call wLog("Changing maxConnect = ")
                call wLog(maxConnect)
            end if

            !Jump the connectivity data lines
            do i = 1, ceiling(dble(connectInfo(6))/8)
                read(fileID, fmt=*) line
                lineNb = lineNb + 1
            end do

        end do

        !!write(get_fileId(),*) "nElem        = ", nElem
        !!write(get_fileId(),*) "maxConnect   = ", maxConnect
        !!write(get_fileId(),*) "monoType     = ", monoType
        !!write(get_fileId(),*) "Elements type"
        do i = 1, size(nElemByType)
            if (nElemByType(i) /= 0) then
                !!write(get_fileId(),*) i, " Nodes :", nElemByType(i), " Elements"
            end if
        end do
        !!!write(get_fileId(),*) "exit Line = ", lineNb

    end subroutine prepareConnectivity

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine preparePhysVol()
    end subroutine

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine readCoordinates(fileID, coordList, nS_in, nE_in, offset_in)
        !INPUT
        integer, intent(in) :: fileID
        integer, intent(in), optional :: nS_in, nE_in
        !OUTPUT
        double precision, dimension(:,:), intent(out) :: coordList
        integer, intent(inout), optional :: offset_in
        !LOCAL
        character (len=200) :: line
        integer :: nNode
        logical :: inside
        !character (len=50), dimension(3) :: strVec_3
        integer :: nS, nE, offset

        !Optional Arguments
        if(.not.present(nS_in)) then
            nS = 1
        else
            nS = nS_in
        end if

        if(.not.present(nE_in)) then
            nE = size(coordList, 2)
        else
            nE = nE_in
        end if

        if(.not.present(offset_in)) then
            offset = 0
        else
            offset = offset_in
        end if

        !Init
        inside = .true.
        nNode = 0

        !Jump initials -1
        read(fileID, fmt=*) line
        lineNb = lineNb + 1
        do while(trim(adjustL(line)) == "-1")
            read(fileID, fmt=*) line
            lineNb = lineNb + 1
        end do
        backspace(fileID)
        lineNb = lineNb - 1

        call wLog("lineNb = ")
        call wLog(lineNb)
        call wLog("shape(coordList) = ")
        call wLog(shape(coordList))
        call wLog("offset = ")
        call wLog(offset)

        !Reading coordinates
        do while(inside)

            !Read coordinates header
            read(fileID, fmt=*) line
            lineNb = lineNb + 1

            if(trim(adjustL(line)) == "-1") then
                !!!write(get_fileId(),*) "Line ", lineNb, "is -1 (inside coord)"
                inside = .false.
            else
                !Reading coordinates values
                nNode = nNode + 1
                if(nNode >= nS .and. nNode <= nE) then
                    !read(fileID, fmt=*) strVec_3
!                    call wLog("Reading Coords ")
!                    call wLog("Node:")
!                    call wLog(nNode)
!                    call wLog("Pos in Vector")
!                    call wLog(nNode-nS+1+offset)
                    read(fileID, fmt=*) coordList(:,nNode-nS+1+offset)
!                    call wLog(coordList(:,nNode-nS+1+offset))
                    lineNb = lineNb + 1
                else
                    read(fileID, fmt=*) line
                    lineNb = lineNb + 1
                end if
            end if
        end do

        if(present(offset_in)) offset_in = offset_in + nNode


    end subroutine readCoordinates

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine readConnectivity(fileID, connectList, sizeList, nS_in, nE_in, monoType, &
                                offset_in, connectStart_in)
        !INPUT
        integer, intent(in) :: fileID
        integer, intent(in), optional :: nS_in, nE_in
        logical, intent(in) :: monoType
        integer, intent(in), optional :: connectStart_in
        !OUTPUT
        integer, dimension(:,:), intent(out) :: connectList
        integer, dimension(:), optional, intent(out) :: sizeList
        integer, intent(inout), optional :: offset_in
        !LOCAL
        character (len=20) :: line
        integer, dimension(6) :: connectInfo
        integer :: i, n, start, end
        logical :: inside
        integer :: nS, nE, offset
        integer :: connectStart

        !Optional Arguments
        if(.not.present(nS_in)) then
            nS = 1
        else
            nS = nS_in
        end if

        if(.not.present(nE_in)) then
            nE = size(connectList, 2)
        else
            nE = nE_in
        end if

        if(.not.present(offset_in)) then
            offset = 0
        else
            offset = offset_in
        end if

        if(.not.present(connectStart_in)) then
            connectStart = 0
        else
            connectStart = connectStart_in
        end if

        connectStart = connectStart_in

        !Init
        inside = .true.
        n = 0

        !Jump initials -1
        read(fileID, fmt=*) line
        lineNb = lineNb + 1
        do while(trim(adjustL(line)) == "-1")
            read(fileID, fmt=*) line
            lineNb = lineNb + 1
        end do
        backspace(fileID)
        lineNb = lineNb - 1

        !write(*,*) "READ CONNECT"
        if((monoType) .and. present(sizeList)) then
            sizeList(:) = size(connectList, 1)
            !!write(get_fileId(),*) "Is monotype"
        end if

        do while(inside)

            read(fileID, fmt=*) line
            lineNb = lineNb + 1

            !!!write(get_fileId(),*) "INSIDE"

            !!!write(get_fileId(),*) "lineNb = ", lineNb
            !!!write(get_fileId(),*) "line   = ", line


            !Read the connectivity header information
            if((trim(adjustL(line)) == "-1")) then
                !!!write(get_fileId(),*) "Line ", lineNb, "is the exit line"
                inside = .false.
                cycle
            else
                n = n + 1
                !!!write(get_fileId(),*) "nElem = ", n
                if(n < nS) then
                    read(fileID, fmt=*) line
                    lineNb = lineNb + 1
                    cycle
                end if
                if(n > nE) then
                    !!!write(get_fileId(),*) "Line ", lineNb, "is an ignored element (Exit line)"
                    inside = .false.
                    cycle
                end if
            end if

            backspace(fileID)
            lineNb = lineNb - 1

            !Read size
            !!!write(get_fileId(),*) "Elem Read"
            read(fileID, fmt=*) connectInfo
            lineNb = lineNb + 1

            !!!write(get_fileId(),*) "Starting on line ", lineNb, "----------------"
            !!!write(get_fileId(),*) "connectInfo = ", connectInfo

            if((.not. monoType) .and. present(sizeList)) then
                sizeList(n-nS+1) = connectInfo(6)
            end if

            !Read the connectivity data lines
            do i = 1, ceiling(dble(connectInfo(6))/8.0D0)
                start = (i-1)*8 + 1
                end   = start + 7
                if(end > connectInfo(6)) end = connectInfo(6)
                !!!write(get_fileId(),*) "Line ", lineNb, "Read "
                read(fileID, fmt=*) connectList(start:end,n-nS+1+offset)
                connectList(start:end,n-nS+1+offset) = connectList(start:end,n-nS+1+offset) + connectStart
                lineNb = lineNb + 1
            end do


        end do

        if(present(offset_in)) offset_in = offset_in + n

        !!!write(get_fileId(),*) "n            = ", n
        !!!write(get_fileId(),*) "exit Line = ", lineNb

        !call dispCarvalhol(sizeList, "sizeList",unit_in = get_fileId())
        !call dispCarvalhol(connectList, "connectList",unit_in = get_fileId())

    end subroutine readConnectivity

end module readUNV_RF
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
