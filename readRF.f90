module readRF
    use displayCarvalhol
contains
    subroutine readFile()
        implicit none

        !File specific
        character (len=30) :: path = 'input06.txt';
        character (len=30) :: comment = '!';
        character :: tagID = '$'
        integer :: fileID;

        !Stocking
        double precision, dimension(:,:), allocatable :: dataTable;
        integer, dimension(:), allocatable :: tagPattern;
        integer, dimension(:), allocatable :: lastLine;

        !Others
        character (len=30) :: foundedTags(10);
        character (len=30) :: empty='';
        double precision :: dataLine;
        character (len=30) :: tag;
        integer :: i, j, k, stat, nLines;
        integer :: tagPos, tagCount, tagTotal;
        integer :: dataRow, dataColumn, dataCount, dataTotal;
        logical :: posTaken, dataPassed, labelPassed;

        stat = 0;
        nLines = 0;
        tagPos = 0;
        tagCount = 0;
        tagTotal = 0;
        posTaken = .FALSE.;
        dataPassed = .FALSE.;
        dataRow = 0;
        dataColumn = 0;
        dataCount = 0;
        dataTotal = 0;
        labelPassed = .FALSE.;
        foundedTags = 'notUsed';

        fileID=2;

        open (unit = fileID , file = path, action = 'read')
            do while (stat.eq.0) !Just to mesure

                nLines = nLines+1;

!                !write(*,*) "nLines = ", nLines;

                read(fileID,'(A)',IOSTAT = stat) tag

                if (tag(1:1).eq.tagID) then !If a tag i founded
                    k=0;
!                    !write(*,*) "Is a tag";
                    posTaken = .FALSE.

                    do j = 1, size(foundedTags) !Founding tag position
!                        !write(*,*) "j = ", j;
                        if(tag.eq.foundedTags(j)) then
                            tagPos = j;
                            posTaken = .TRUE.
!                            !write(*,*) "tag founded in = ", j;
                            exit
                        endif
                    enddo
                    !write(*,*) "posTaken = ", posTaken;

                    if(posTaken .eqv. .FALSE.) then
                        tagTotal = tagTotal +1
                        foundedTags(tagTotal) = tag
                        posTaken = .TRUE.
!                        !write(*,*) "tagTotal = ", tagTotal;
!                        !write(*,*) "foundedTags = ", foundedTags;
                    endif

                else if (tag(1:1).eq.comment .or. tag(1:1).eq.empty) then
                    k=0;
                else
                    dataTotal = dataTotal+1
                    k=1;
!                    !write(*,*) "dataTotal = ", dataTotal;
                endif

            enddo

            if(mod((dataTotal-k),tagTotal).eq.0) then
                dataTotal = (dataTotal-k)/tagTotal; !obs: -k because of the last line
                !write(*,*) "dataTotalLines = ", dataTotal;
                rewind(fileID)
            else
                write(*,*) "In readRF, Tags and data dimensions don't match ";
                !Put abort routine here
            endif

            write(*,*) "dataTotal = ", dataTotal;
            write(*,*) "tagTotal = ", tagTotal;

            !Initializing
            allocate(dataTable(dataTotal, tagTotal));
            !allocate(foundedTags(40));
            allocate(tagPattern(tagTotal));
            allocate(lastLine(size(dataTable,2)));
            stat = 0;
            tagTotal = 0;
            dataTotal = 0;
            tagPattern = 0;
            lastLine = 0;
            dataTable = 0;
            tagPos = 0;
            posTaken = .FALSE.;

            call Disp2D(dataTable)

            do i = 1, (nLines-1)
                read(fileID,'(A)',IOSTAT = stat) tag

                if (tag(1:1).eq.tagID) then !If a tag i founded
                    posTaken = .FALSE.; !tag colomn index is yet to be decided
                    labelPassed = .TRUE.; !says we have put a label (so now we can put data)


                    if(dataPassed.eqv..TRUE.) then !If it's a new pattern we reset the variables
                        tagPattern = 0; !Vector indicating all tags columns position
                        tagCount = 0; !Integer indicating this tag column position
                        dataPassed = .FALSE.; !Says that it'll be a new data packet
                        dataCount = 0;
                    endif

                    tagCount = tagCount +1;

                    do j = 1, size(foundedTags) !Founding tag position
                        if(tag.eq.foundedTags(j)) then
                            tagPos = j;
                            posTaken = .TRUE.
                            exit
                        endif
                    enddo

                    if(posTaken.eqv..FALSE.) then
                        tagPos = tagTotal + 1; !Take the next position
                        foundedTags(tagPos) = tag;
                        tagTotal = tagTotal + 1; !Last tag index (global)
                        posTaken = .TRUE.;
                    endif

                    if(posTaken.eqv..TRUE.) then
                        tagPattern(tagCount) = tagPos !Put the position found in the pattern (in case we have more than one tag)
                    else
                        print*, "In ReadRF - ERROR - Position was not taken"
                    endif

                else if (tag(1:1).eq.comment .or. tag(1:1).eq.empty) then !If the first character is a space or it's a comment we ignore it

                else if(labelPassed.eqv..TRUE.)then !Its consider data and we know its tag
                    backspace(fileID)
                    dataLine = 0.0;
                    read(fileID,*,IOSTAT = stat) dataLine
                    dataPassed = .TRUE.;
                    dataCount = dataCount + 1;

                    !write(*,*) "modCyclic(dataCount,tagCount)= ", modCyclic(dataCount,tagCount);
                    !write(*,*) "tagPattern(modCyclic(dataCount,tagCount)) = ", tagPattern(modCyclic(dataCount,tagCount));

                    dataColumn = tagPattern(modCyclic(dataCount,tagCount));
                    lastLine(dataColumn) =  lastLine(dataColumn) +1;
                    dataTable(lastLine(dataColumn), dataColumn) = dataLine

                    !write(*,*) "Data line = ", dataLine;
                    !write(*,*) "dataRow = ", lastLine(dataColumn);
                    !write(*,*) "dataColumn = ", dataColumn;
                    call Disp2D(dataTable)
                endif
            enddo
        close(fileID)

    end subroutine

    function modCyclic(dividend,divisor) result(cyclicRest)
        implicit none
        integer :: dividend, divisor, cyclicRest

        if (mod(dividend, divisor).eq.0) then
            cyclicRest = divisor;
        else
            cyclicRest = mod(dividend, divisor);
        endif
    end function modCyclic

end module readRF
