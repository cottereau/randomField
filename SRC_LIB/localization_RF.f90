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

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_neighbours_info (RDF, MSH)
        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF

        !INPUT
        type(MESH), intent(in) :: MSH

        !LOCAL
        integer :: neighPos, stage, direction
        integer, allocatable, dimension(:) :: tempSeed
        integer :: request, code, tag
        integer, dimension(MPI_STATUS_SIZE) :: status

        do neighPos = 1, size(MSH%neigh)
            if(MSH%neigh(neighPos) < 0) cycle !Check if this neighbour exists
            call calculate_random_seed(tempSeed, RDF%seedStart+MSH%neigh(neighPos))
            RDF%neighSeed(:,neighPos) = tempSeed
        end do

        do stage = 1, 2 !Sending and then receiving

            do neighPos = 1, size(MSH%neigh)

                if(MSH%neigh(neighPos) < 0) cycle !Check if this neighbour exists

                if(stage == 1) then
                    tag = findTag(MSH, neighPos, neighPos, send = .true.)
                    call MPI_ISEND (RDF%xMaxBound(:)-RDF%xMinBound(:), RDF%nDim, MPI_DOUBLE_PRECISION, &
                        MSH%neigh(neighPos), tag, RDF%comm, request, code)
                else if(stage == 2) then
                    tag = findTag(MSH, neighPos, neighPos, send = .false.)
                    call MPI_RECV (RDF%neighRange(:,neighPos), RDF%nDim, MPI_DOUBLE_PRECISION, &
                        MSH%neigh(neighPos), tag, RDF%comm, status, code)
                end if

            end do

        end do

        deallocate(tempSeed)

    end subroutine get_neighbours_info

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine getNeighIndexRange(MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)
        !INPUT
        type(MESH), intent(in) :: MSH
        !OUTPUT
        integer, intent(out) :: minIndexNeigh, maxIndexNeigh
        logical, dimension(:), intent(out) :: considerNeighbour

        considerNeighbour = .true.
        where(MSH%neigh < 0) considerNeighbour = .false.
        where(minval(MSH%neighShift,1) < 0) considerNeighbour = .false.

        !Global Min and Max positions
        minIndexNeigh = minval(pack(MSH%indexNeigh(1,:), considerNeighbour))
        maxIndexNeigh = maxval(pack(MSH%indexNeigh(2,:), considerNeighbour))

    end subroutine getNeighIndexRange

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine applyWeightingFunctions(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)
        !INPUT:
        type(RF), intent(in) :: RDF
        type(MESH), intent(in) :: MSH
        integer, intent(in) :: minIndexNeigh, maxIndexNeigh
        logical, dimension(:), intent(in) ::considerNeighbour

        !LOCAL
        integer :: direction, minPos, maxPos, i
        double precision, dimension(minIndexNeigh:maxIndexNeigh) :: powerWF
        double precision, dimension(MSH%nDim) :: originCorner

        !Modify extremes of local Random Field-------------------------------------------------------

        !Loop on directions
        do direction = 1, size(MSH%neigh)

            if(.not. considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction

            !Positions in temp Vector
            minPos = MSH%indexNeigh(1,direction)
            maxPos = MSH%indexNeigh(2,direction)

            !Finding origin
            originCorner = MSH%xMinNeigh(:, direction)
            where(MSH%neighShift(:, direction) == -1) originCorner = MSH%xMaxNeigh(:, direction)
            !!write(get_fileId(), *) "  originCorner = ", originCorner

            !Shape Function multiplication (Obs: redefinition of power, should come after Normalization values calculation)
            powerWF(minPos:maxPos) = 0.0D0
            do i = 1, MSH%nDim
                if(MSH%neighShift(i, direction) == 0) cycle

                powerWF(minPos:maxPos) = ((RDF%xPoints(i, minPos:maxPos) - originCorner(i))/MSH%overlap(i))**2 &
                                         + powerWF(minPos:maxPos)
            end do

            powerWF(minPos:maxPos) = sqrt(powerWF(minPos:maxPos))
            where (powerWF(minPos:maxPos) > 1.0D0) powerWF(minPos:maxPos) = 1.0D0

        end do !Direction

        RDF%randField(minIndexNeigh:maxIndexNeigh,1) = RDF%randField(minIndexNeigh:maxIndexNeigh,1) &
                                                 * sqrt((cos(PI*(powerWF))+1.0D0)/2.0D0)

    end subroutine applyWeightingFunctions

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine normalizeOverlap(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)

        implicit none

        !INPUT
        type(RF), intent(inout) :: RDF
        type(MESH), intent(in) :: MSH
        integer, intent(in) :: minIndexNeigh, maxIndexNeigh
        logical, dimension(:), intent(in) ::considerNeighbour

        !OUTPUT
        double precision, dimension(minIndexNeigh:maxIndexNeigh) :: normFactor

        !LOCAL
        integer :: i, j, n
        double precision, dimension(MSH%nDim, 2**MSH%nDim) :: originList
        double precision, dimension(minIndexNeigh:maxIndexNeigh) :: power
        integer, dimension(MSH%nDim) :: tempXNStep
        double precision :: distance
        integer :: minPos, maxPos
        integer :: direction, zeroVal

        !Find normalization values-------------------------------------------------------
        normFactor(:) = 0
        power(:) = 0

        !Loop on directions
        !write(get_fileId(), *) "   Finding normalization Factors "
        do direction = 1, size(MSH%neigh)

            if(.not. considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction

            !Positions in temp Vector
            minPos = MSH%indexNeigh(1,direction)
            maxPos = MSH%indexNeigh(2,direction)



            !Finding all vertices
            tempXNStep = 2
            do i = 1, 2**MSH%nDim
                call get_Permutation(i, MSH%xMaxNeigh(:,direction), tempXNStep, originList(:, i), &
                                     MSH%xMinNeigh(:,direction), snapExtremes = .true.);
            end do

            !Normalization values
            normFactor(minPos:maxPos) = 0 !TO DELETE
            zeroVal = MSH%nDim - (sum((MSH%neighShift(:, direction))**2))

            !We sweep all vertices in this direction
            do j = 1, size(originList, 2)
                power(minPos:maxPos) = 0
                do i = 1, MSH%nDim
                    if(MSH%neighShift(i, direction) == 0) cycle

                    power(minPos:maxPos) = ((RDF%xPoints(i, minPos:maxPos) - originList(i, j))/MSH%overlap(i))**2 &
                                               + power(minPos:maxPos);
                end do

                power(minPos:maxPos) = sqrt(power(minPos:maxPos))
                where (power(minPos:maxPos) > 1.0D0) power(minPos:maxPos) = 1.0D0
                normFactor(minPos:maxPos) = sqrt((cos((PI)*(power(minPos:maxPos)))+1.0D0)/2.0D0) + normFactor(minPos:maxPos)


            end do !Neighbours

            normFactor(minPos:maxPos) = normFactor(minPos:maxPos)/(2.0**zeroVal)

        end do !Direction

        !RDF%randField!Normalizing Field---------------------------------------------------------------
        RDF%randField(minIndexNeigh:maxIndexNeigh,1) = RDF%randField(minIndexNeigh:maxIndexNeigh,1)/normFactor(:)

    end subroutine normalizeOverlap

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

        orig  = (dble(RDF%origin-1)*MSH%xStep + MSH%xMinGlob)
        xStep = MSH%xStep

        offset(1) = 1
        do i = 2, MSH%nDim
          offset(i) = product(MSH%xNStep(1:i-1))
        end do

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
            write(*,*) "Inside findTag , Invalid Neighbour"
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

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine takeNeighboursContribution(RDF, MSH, minIndexNeigh, maxIndexNeigh, considerNeighbour)
        !INPUT
        type(MESH), intent(in) :: MSH
        integer, intent(in) :: minIndexNeigh, maxIndexNeigh
        logical, dimension(:), intent(in) ::considerNeighbour
        !OUTPUT
        type(RF), intent(inout) :: RDF
        !LOCAL
        double precision, dimension(MSH%nDim) :: neighOrCorner, originCorner
        logical :: sndRcv
        integer :: i, direction, neighPos, minPos, maxPos
        double precision, dimension(minIndexNeigh:maxIndexNeigh) :: power
        type(RF) :: tmpRDF


        !Taking the contributions from neighbours------------------------------------------------------
        !write(get_fileId(), *) "   Taking the contributions from neighbours "
        call init_RF(tmpRDF, RDF%nDim, RDF%Nmc, RDF%comm, RDF%rang, RDF%nb_procs)
        call copy_RF_properties(RDF, tmpRDF)

        do direction = 1, size(MSH%neigh)

            if(.not. considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction

            !write(get_fileId(), *) " "
            !write(get_fileId(), *) "  DIRECTION      = ", direction
            !write(get_fileId(), *) "  Neighbour Rank = ", MSH%neigh(direction)

            originCorner = MSH%xMinNeigh(:, direction)
            where(MSH%neighShift(:, direction) == -1) originCorner = MSH%xMaxNeigh(:, direction)


            !!write(get_fileId(),*) "Direction = ", direction

            !Preparing the xPoints of a given direction
            call copy_RF_xPoints(MSH, RDF, tmpRDF, tmpRDF%xPoints_Local, direction)
            call allocate_randField(tmpRDF, tmpRDF%randField_Local)
            minPos = MSH%indexNeigh(1,direction)
            maxPos = MSH%indexNeigh(2,direction)

            do neighPos = 1, size(MSH%neigh)

                if(.not. considerNeighbour(neighPos)) cycle !Don't consider Neighbours in this direction

                sndRcv = .true.

                !Checking if this neighbour should be taking into account in this part of the field
                do i = 1, MSH%nDim
                    if(       (MSH%neighShift(i, neighPos) /= 0) &
                        .and. (MSH%neighShift(i, neighPos) /= MSH%neighShift(i, direction))) then
                        sndRcv = .false.
                        exit
                    end if
                end do

                !From this point we know that we want the contribution of this neighbour in this direction
                if (sndRcv) then

                    tmpRDF%xMinBound = 0.0D0
                    tmpRDF%xMaxBound = RDF%neighRange(:,neighPos)
                    tmpRDF%seed      = RDF%neighSeed(:,neighPos)

                    !Generating Standard Gaussian Field
                    select case (tmpRDF%method)
                        case(SHINOZUKA)
                            call gen_Std_Gauss_Shinozuka(tmpRDF)
                        case(ISOTROPIC)
                            call gen_Std_Gauss_Isotropic(tmpRDF)
                        case(RANDOMIZATION)
                            call gen_Std_Gauss_Randomization(tmpRDF)
                    end select

                    !Finding origin for Shape Function
                    neighOrCorner = originCorner + MSH%overlap*MSH%neighShift(:, neighPos)

                    !Shape Function multiplication and sum of the contribution
                    power(minPos:maxPos) = 0
                    do i = 1, MSH%nDim
                        if(MSH%neighShift(i, direction) == 0) cycle

                        power(minPos:maxPos) = ((tmpRDF%xPoints(i, :) - neighOrCorner(i))/MSH%overlap(i))**2 &
                                                   + power(minPos:maxPos)
                    end do

                    power(minPos:maxPos) = sqrt(power(minPos:maxPos))
                    where (power > 1.0D0) power = 1.0

                    RDF%randField(minPos:maxPos,1) = RDF%randField(minPos:maxPos,1) + &
                                                     (tmpRDF%randField(:,1) * sqrt((cos((PI)*(power(minPos:maxPos)))+1.0D0)/2.0D0))
                end if
            end do !Neighbours
        end do !Directions
        call finalize_RF(tmpRDF)

    end subroutine takeNeighboursContribution

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine copy_RF_xPoints(MSH, orRDF, destRDF, destXPoints, direction)
        implicit none

        !INPUT AND OUTPUT
        type(MESH), intent(in) :: MSH
        type(RF), intent(in)   :: orRDF
        type(RF) ::destRDF
        integer, intent(in) :: direction
        double precision, dimension(:, :), allocatable, intent(out), target :: destXPoints;

        !LOCAL
        integer :: i, j, totalSize

        totalSize = MSH%indexNeigh(2,direction) - MSH%indexNeigh(1,direction) + 1

        if(allocated(destXPoints)) then
            if(size(destXPoints) /= totalSize) then
                nullify(destRDF%xPoints)
                deallocate(destXPoints)
            end if
        end if

        if(.not.allocated(destXPoints)) allocate(destXPoints(MSH%nDim, 1:totalSize))

        destRDF%xPoints => destXPoints

        destRDF%xPoints(:,:) = orRDF%xPoints(:, MSH%indexNeigh(1,direction):MSH%indexNeigh(2,direction))
        destRDF%xNTotal = totalSize

    end subroutine copy_RF_xPoints

end module localization_RF
