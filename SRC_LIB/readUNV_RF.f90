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
    subroutine readUNV(path, nDim, coordList, connectList, monotype, rang, nb_procs, comm)
        !INPUT
        character (len=*), intent(in) :: path
        integer, intent(in)           :: nDim, rang, nb_procs, comm

        !OUTPUT
        double precision, dimension(:,:), allocatable, intent(out), optional :: coordList
        integer         , dimension(:,:), allocatable, intent(out), optional :: connectList
        logical, optional, intent(out) :: monoType

        !LOCAL
        integer :: fileID = 600
        character (len=200) :: line
        integer :: stat, code
        integer :: i
        logical :: inside
        integer :: nNodes
        integer :: nNodesLoc, nodeStart, nodeEnd
        integer :: nElemLoc, elemStart, elemEnd
        integer :: nElem, maxConnect
        integer, allocatable, dimension(:) :: elemSizes
        integer, dimension(:)  , allocatable :: sizeList

        integer, dimension(:,:), allocatable :: startEnd

        stat     = 0;
        lineNb = 0
        inside = .false.

        open (unit = fileID , file = path, action = 'read')

            !Dimensioning lecture
            call wLog("--------------------------------------------------")
            call wLog("DIMENSIONING LECTURE------------------------------")
            call wLog("--------------------------------------------------")
            call wLog("")

            allocate(startEnd(2, nb_procs))

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
        integer, intent(out) :: maxConnect
        logical, intent(out) :: monoType
        !LOCAL
        integer, dimension(50) :: nElemByType
        character (len=20) :: line
        integer, dimension(6) :: connectInfo
        integer :: i
        logical :: inside

        inside = .true.
        nElem = 0
        monoType = .true.
        maxConnect = 0
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
            if(maxConnect == 0) maxConnect = connectInfo(6)
            if(maxConnect /= connectInfo(6)) then
                if(maxConnect < connectInfo(6)) maxConnect = connectInfo(6)
                monoType = .false.
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
    subroutine readCoordinates(fileID, coordList, nS_in, nE_in)
        !INPUT
        integer, intent(in) :: fileID
        integer, intent(in), optional :: nS_in, nE_in
        !OUTPUT
        double precision, dimension(:,:), intent(out) :: coordList
        !LOCAL
        character (len=200) :: line
        integer :: nNode
        logical :: inside
        !character (len=50), dimension(3) :: strVec_3
        integer :: nS, nE

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
                    read(fileID, fmt=*) coordList(:,nNode-nS+1)
                    lineNb = lineNb + 1
                else
                    read(fileID, fmt=*) line
                    lineNb = lineNb + 1
                end if
            end if
        end do


    end subroutine readCoordinates

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine readConnectivity(fileID, connectList, sizeList, nS_in, nE_in, monoType)
        !INPUT
        integer, intent(in) :: fileID
        integer, intent(in), optional :: nS_in, nE_in
        logical, intent(in) :: monoType
        !OUTPUT
        integer, dimension(:,:), intent(out) :: connectList
        integer, dimension(:), optional, intent(out) :: sizeList
        !LOCAL
        character (len=20) :: line
        integer, dimension(6) :: connectInfo
        integer :: i, n, start, end
        logical :: inside
        integer :: nS, nE

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
                read(fileID, fmt=*) connectList(start:end,n-nS+1)
                lineNb = lineNb + 1
            end do


        end do

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
