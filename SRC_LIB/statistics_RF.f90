module statistics_RF

    use displayCarvalhol
    use math_RF
    use mpi
    use write_Log_File
    use type_STAT

    implicit none
    include 'fftw3.f'
    !WARNING just before this line we have include 'fftw3.f'

contains
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine calculate_average_and_stdVar_MPI(STA)

        implicit none

        !INPUT
        type(STAT) :: STA

        !LOCAL
        double precision, dimension(:), allocatable :: sumRF, sumRFsquare
        double precision, dimension(:), allocatable :: totalSumRF, totalSumRFsquare;
        double precision, dimension(:), allocatable :: sumRF_point, sumRFsquare_point
        double precision, dimension(:), allocatable :: totalSumRF_point, totalSumRFsquare_point
        integer :: Nmc, xNTotal
        !integer, dimension(:), allocatable :: xNTotal_Vec, deplacement
        integer :: code, nb_procs, rang, comm
        integer :: i

        write(*,*) "Calculating Average and stdVar"

        Nmc     = size(STA%randField, 2)
        xNTotal = size(STA%randField, 1)
        comm    = STA%comm

        !Allocating
        allocate(sumRF(Nmc))
        allocate(sumRFsquare(Nmc))
        allocate(totalSumRF(Nmc))
        allocate(totalSumRFsquare(Nmc))

        allocate(sumRF_point(xNTotal))
        allocate(sumRFsquare_point(xNTotal))
        allocate(totalSumRF_point(xNTotal))
        allocate(totalSumRFsquare_point(xNTotal))

        !allocate(xNTotal_Vec(nb_procs))
        !allocate(deplacement(nb_procs))

        !Calculating
        write(*,*) "sumRF(1) = ", sumRF(1)
        write(*,*) "sumRFsquare(1) = ", sumRFsquare(1)
        !write(*,*) "sumRF_point(1) = ", sumRF_point(1)
        !write(*,*) "sumRFsquare_point(1) = ", sumRFsquare_point(1)

        !Total Number of Points
        write(*,*) "xNTotal = ", xNTotal

        !By Event
        sumRF(:)       = sum( STA%randField    , dim = 1)
        sumRFsquare(:) = sum((STA%randField)**2, dim = 1)
        call MPI_ALLREDUCE (sumRF,totalSumRF,Nmc,MPI_DOUBLE_PRECISION, &
            MPI_SUM,comm,code)
        call MPI_ALLREDUCE (sumRFsquare,totalSumRFsquare,Nmc,MPI_DOUBLE_PRECISION, &
            MPI_SUM,comm,code)
        STA%evntAvg      = totalSumRF/xNTotal;
        STA%evntStdDev   = sqrt(totalSumRFsquare/dble(xNTotal) &
            - (STA%evntAvg)**2)

        !Global
        STA%globalAvg    = sum(totalSumRF)/dble(xNTotal*Nmc);
        STA%globalStdDev = sqrt(sum(totalSumRFsquare)/dble(xNTotal*Nmc) &
            - (STA%globalAvg)**2)

        !        !By Point
        !        sumRF_point(:)       = sum( STA%randField    , dim = 2)
        !        sumRFsquare_point(:) = sum((STA%randField)**2, dim = 2)
        !        call MPI_ALLREDUCE (sumRF_point,totalSumRF_point,xNTotal,MPI_DOUBLE_PRECISION, &
        !                            MPI_SUM,comm,code)
        !        call MPI_ALLREDUCE (sumRFsquare_point,totalSumRFsquare_point,xNTotal,MPI_DOUBLE_PRECISION, &
        !                            MPI_SUM,comm,code)


        !        write(*,*) "HERE Point"
        !        write(*,*) "shape(sumRF_point) = ", shape(sumRF_point)
        !        write(*,*) "shape(totalSumRF_point) = ", shape(totalSumRF_point)
        !        write(*,*) "xNTotal = ", xNTotal
        !        !Number of points
        !        !write(*,*) "HERE Number of Points"
        !        call MPI_ALLREDUCE (xNTotal, all_xNTotal,1,MPI_INTEGER, &
        !                            MPI_SUM,comm,code)
        !
        !
        !
        !        !write(*,*) "xNTotal_Vec = ", xNTotal_Vec
        !        call MPI_ALLREDUCE (sumRFsquare_point,totalSumRFsquare_point,xNTotal,MPI_DOUBLE_PRECISION, &
        !                            MPI_SUM,comm,code)
        !

        !        !write(get_fileId(),*) "totalSumRF(1) = ", totalSumRF(1)
        !        !write(get_fileId(),*) "totalSumRFsquare(1) = ", totalSumRFsquare(1)
        !        !write(get_fileId(),*) "totalSumRF_point(1) = ", totalSumRF_point(1)
        !        !write(get_fileId(),*) "totalSumRFsquare_point(1) = ", totalSumRFsquare_point(1)
        !        !write(get_fileId(),*) "all_xNTotal = ", all_xNTotal
        !
        !        !by Event
        !        if(present(evntAvg))    evntAvg      = totalSumRF/dble(all_xNTotal);
        !        if(present(evntStdDev)) evntStdDev   = sqrt(totalSumRFsquare/dble(all_xNTotal) &
        !                                               - (evntAvg)**2)
        !
        !        !by Point
        !        if(present(pointAvg) .or. present(pointStdDev)) then
        !            deplacement(1) = 0
        !            do i = 2, nb_procs
        !                deplacement(i) = sum(xNTotal_Vec(1:i-1))
        !            end do
        !
        !            if(rang == 0) write(*,*) "xNTotal_Vec =", xNTotal_Vec
        !            if(rang == 0) write(*,*) "deplacement =", deplacement
        !        end if
        !
        !        if(present(pointAvg)) then
        !            !allocate(pointAvg(all_xNTotal))
        !            call MPI_ALLGATHERV (sumRF_point, xNTotal, MPI_DOUBLE_PRECISION, &
        !                                 pointAvg, xNTotal_Vec, deplacement, MPI_DOUBLE_PRECISION, &
        !                                 comm, code)
        !            pointAvg(:) =  pointAvg(:)/Nmc
        !        end if
        !
        !        if(present(pointStdDev).and.present(pointAvg)) then
        !            !allocate(pointStdDev(all_xNTotal))
        !            call MPI_ALLGATHERV (sumRFsquare_point, xNTotal, MPI_DOUBLE_PRECISION, &
        !                                 pointStdDev, xNTotal_Vec, deplacement, MPI_DOUBLE_PRECISION, &
        !                                 comm, code)
        !
        !            pointStdDev   = sqrt(pointStdDev/dble(Nmc) &
        !                            - (pointAvg)**2)
        !        end if
        !
        !


        !Deallocating
        if(allocated(sumRF))       deallocate(sumRF)
        if(allocated(sumRFsquare)) deallocate(sumRFsquare)
        if(allocated(totalSumRF))       deallocate (totalSumRF)
        if(allocated(totalSumRFsquare)) deallocate (totalSumRFsquare)

        if(allocated(sumRF_point)) deallocate(sumRF_point)
        if(allocated(sumRFsquare_point)) deallocate(sumRFsquare_point)
        if(allocated(totalSumRF_point)) deallocate(totalSumRF_point)
        if(allocated(totalSumRFsquare_point)) deallocate(totalSumRFsquare_point)

        !if(allocated(xNTotal_Vec)) deallocate(xNTotal_Vec)
    !        if(allocated(deplacement)) deallocate(deplacement)
    !
    end subroutine calculate_average_and_stdVar_MPI

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine rebuild_Sk(STA)
        implicit none
        !INPUT
        type(STAT) :: STA
        !LOCAL
        integer :: i, comm, code
        integer, dimension(STA%nDim) :: nPointsMin
        integer :: delta
        double precision, dimension(:), allocatable :: Temp_SkTot_Dir

        write(*,*) "rebuild_Sk"
        do i = 1, STA%nDim
            call rebuild_Sk_FFT(STA%randField(:,1), STA%xNStep_Loc, i, STA%nDim, &
                                STA%Sk_Dir(STA%Sk_Ind(i,1):STA%Sk_Ind(i,2)),     &
                                STA%globalAvg)
        end do

        write(*,*) "Reducing Sk to one proc"
        write(*,*) "STA%xNStep_Loc = ", STA%xNStep_Loc
        do i = 1, STA%nDim
            call MPI_ALLREDUCE (STA%xNStep_Loc(i), nPointsMin(i), 1, MPI_INTEGER, &
                                MPI_MIN, STA%comm, code)
        end do

        write(*,*) "nPointsMin BEFORE = ", nPointsMin

        nPointsMin = nPointsMin/2; !We'll only save half of the correlation function (symmetric)
        write(*,*) "nPointsMin = ", nPointsMin

        allocate(STA%SkTot_Dir(sum(nPointsMin)))

        do i = 1, STA%nDim
            if(i == 1) then
                STA%SkTot_Ind(i,1) = 1
            else
                STA%SkTot_Ind(i,1) = sum(nPointsMin(1:i-1)) + 1
            end if
            STA%SkTot_Ind(i,2) = sum(nPointsMin(1:i))

            delta = STA%SkTot_Ind(i,2) - STA%SkTot_Ind(i,1)

            STA%SkTot_Dir(STA%SkTot_Ind(i,1):STA%SkTot_Ind(i,2)) = STA%Sk_Dir(STA%Sk_Ind(i,1):STA%Sk_Ind(i,1)+delta)

            !call dispCarvalhol(STA%SkTot_Dir(STA%SkTot_Ind(i,1):STA%SkTot_Ind(i,2)), "STA%SkTot_Dir(STA%SkTot_Ind(i,1):STA%SkTot_Ind(i,2))")
        end do

        !write(*,*) "shape(STA%SkTot_Dir) = ", shape(STA%SkTot_Dir)
        !if(STA%rang == 0) allocate(Temp_SkTot_Dir(size(STA%SkTot_Dir)))

        allocate(Temp_SkTot_Dir(size(STA%SkTot_Dir)))
        call MPI_ALLREDUCE (STA%SkTot_Dir, Temp_SkTot_Dir, size(STA%SkTot_Dir), &
                            MPI_DOUBLE_PRECISION, MPI_SUM, STA%comm, code)
        STA%SkTot_Dir = Temp_SkTot_Dir/dble(STA%nb_procs)

        !call MPI_REDUCE (STA%SkTot_Dir, Temp_SkTot_Dir, size(STA%SkTot_Dir), &
        !                 MPI_DOUBLE_PRECISION, MPI_SUM, 0, STA%comm, code)
        !if(STA%rang == 0) then
        !    STA%SkTot_Dir = Temp_SkTot_Dir/dble(STA%nb_procs)
        !    write(*,*) "STA%SkTot_Dir = ", STA%SkTot_Dir
        !end if

        if(allocated(Temp_SkTot_Dir)) deallocate(Temp_SkTot_Dir)

!        !Printing
!        if(STA%rang == 0) then
!            do i = 1, STA%nDim
!                call dispCarvalhol(STA%SkTot_Dir(STA%SkTot_Ind(i,1):STA%SkTot_Ind(i,2)), &
!                                   "STA%SkTot_Dir(STA%SkTot_Ind(i,1):STA%SkTot_Ind(i,2))", unit_in=6)
!            end do
!        end if

    end subroutine rebuild_Sk

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine rebuild_Sk_FFT(randFieldVec, xNStep, dir, nDim, Sk_part, avg)
        implicit none
        !INPUT
        double precision, dimension(:), intent(in), target :: randFieldVec
        integer(kind=8), dimension(:), intent(in) :: xNStep
        integer, intent(in) :: dir, nDim
        double precision, intent(in) :: avg
        !OUTPUT
        double precision, dimension(:), intent(out) :: Sk_part
        !LOCAL
        double complex  , dimension(size(randFieldVec)), target :: SkVec
        !double precision, dimension(:,:)  , pointer :: RF_2D
        !double precision, dimension(:,:,:), pointer :: RF_3D
        double complex  , dimension(:,:)  , pointer :: Sk_2D
        double complex  , dimension(:,:,:), pointer :: Sk_3D
        integer*8 plan
        integer :: howMany, sizeDir, i, j, k
        !integer, dimension(4,6,8) :: testeMat

        !write(*,*) "Calculating Sk"

        SkVec = randFieldVec - avg;

        if(nDim == 2) then
            !RF_2D(1:xNStep(1),1:xNStep(2)) => randFieldVec
            Sk_2D(1:xNStep(1),1:xNStep(2)) => SkVec
        else if(nDim == 3) then
            !RF_3D(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => randFieldVec
            Sk_3D(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => SkVec
        end if
        !
        !Making FFT

        howMany = product(xNStep)/xNStep(dir)
        !sizeDir = size(Sk_2D,dir)
        sizeDir = xNStep(dir)

        !write(*,*) "Making FFT"
        if(nDim == 1) then
            call dfftw_plan_dft_1d(plan, sizeDir, SkVec(:), SkVec(:), &
                FFTW_FORWARD, FFTW_ESTIMATE)
            call dfftw_execute_dft(plan, SkVec(:), SkVec(:))
            call dfftw_destroy_plan(plan)

        else if(nDim == 2) then
            if(dir == 1) then
                do j = 1, xNStep(2)
                    call dfftw_plan_dft_1d(plan, sizeDir, Sk_2D(:,j), Sk_2D(:,j), &
                        FFTW_FORWARD, FFTW_ESTIMATE)
                    call dfftw_execute_dft(plan, Sk_2D(:, j), Sk_2D(:, j))
                    call dfftw_destroy_plan(plan)
                end do

            else if(dir == 2) then
                do i = 1, xNStep(1)
                    call dfftw_plan_dft_1d(plan, sizeDir, Sk_2D(i,:), Sk_2D(i,:), &
                        FFTW_FORWARD, FFTW_ESTIMATE)
                    call dfftw_execute_dft(plan, Sk_2D(i,:), Sk_2D(i,:))
                    call dfftw_destroy_plan(plan)
                end do
            end if


        else if(nDim == 3) then
            if(dir == 1) then
                do j = 1, xNStep(2)
                    do k = 1, xNStep(3)
                        call dfftw_plan_dft_1d(plan, sizeDir, Sk_3D(:,j,k), Sk_3D(:,j,k), &
                            FFTW_FORWARD, FFTW_ESTIMATE)
                        call dfftw_execute_dft(plan, Sk_3D(:,j,k), Sk_3D(:,j,k))
                        call dfftw_destroy_plan(plan)
                    end do
                end do

            else if(dir == 2) then
                do i = 1, xNStep(1)
                    do k = 1, xNStep(3)
                        call dfftw_plan_dft_1d(plan, sizeDir, Sk_3D(i,:,k), Sk_3D(i,:,k), &
                            FFTW_FORWARD, FFTW_ESTIMATE)
                        call dfftw_execute_dft(plan, Sk_3D(i,:,k), Sk_3D(i,:,k))
                        call dfftw_destroy_plan(plan)
                    end do
                end do

            else if(dir == 3) then
                do i = 1, xNStep(1)
                    do j = 1, xNStep(2)
                        call dfftw_plan_dft_1d(plan, sizeDir, Sk_3D(i,j,:), Sk_3D(i,j,:), &
                            FFTW_FORWARD, FFTW_ESTIMATE)
                        call dfftw_execute_dft(plan, Sk_3D(i,j,:), Sk_3D(i,j,:))
                        call dfftw_destroy_plan(plan)
                    end do
                end do
            end if

        end if

        !write(*,*) "Multiplying by the conjugate"
        if(nDim == 1) then
            SkVec = SkVec * conjg(SkVec)
        else if (nDim == 2) then
            Sk_2D = Sk_2D * conjg(Sk_2D)
        else if (nDim == 3) then
            Sk_3D = Sk_3D * conjg(Sk_3D)
        else
            stop("ERROR!!! Inside rebuild_Sk_FFT nDim not accepted")
        end if

        !write(*,*) "Making Average"
        !OBS: Redefinition of SkVec for simplicity

        if(nDim == 1) then
            SkVec(1:xNStep(dir)) = SkVec
        else if (nDim == 2) then
            if(dir == 1) then
                SkVec(1:xNStep(dir)) = sum(Sk_2D, 2)/size(Sk_2D,2)
            else if(dir == 2) then
                SkVec(1:xNStep(dir)) = sum(Sk_2D, 1)/size(Sk_2D,1)
            end if
        else if (nDim == 3) then
            if(dir == 1) then
                SkVec(1:xNStep(dir)) = sum((sum(Sk_3D, 2)/size(Sk_3D,2)), 2)/size(Sk_3D,3)
            else if(dir == 2) then
                SkVec(1:xNStep(dir)) = sum((sum(Sk_3D, 1)/size(Sk_3D,1)), 2)/size(Sk_3D,3)
            else if(dir == 3) then
                SkVec(1:xNStep(dir)) = sum((sum(Sk_3D, 1)/size(Sk_3D,1)), 1)/size(Sk_3D,2)
            end if
        else
            stop("ERROR!!! Inside rebuild_Sk_FFT nDim not accepted")
        end if

        !Reverting FFT
        !write(*,*) "Reverting FFT"

        call dfftw_plan_dft_1d(plan, size(SkVec(1:xNStep(dir))), SkVec(1:xNStep(dir)), SkVec(1:xNStep(dir)), &
            FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, SkVec(1:xNStep(dir)), SkVec(1:xNStep(dir)))
        call dfftw_destroy_plan(plan)

        !Taking real part
        Sk_part = SkVec(1:xNStep(dir))
        Sk_part = Sk_part/Sk_part(1)

        !write(*,*) "Sk_part = ", Sk_part

        if(associated(Sk_2D)) nullify(Sk_2D)
        if(associated(Sk_3D)) nullify(Sk_3D)

    end subroutine rebuild_Sk_FFT

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine rebuild_corrL(STA, corrL_out)
        implicit none
        !INPUT
        type(STAT) :: STA
        !OUTPUT
        double precision, dimension(:), intent(out) :: corrL_out
        !LOCAL
        integer :: i

        do i = 1, STA%nDim
            corrL_out(i) = 2.0D0 * &
                           sum(STA%SkTot_Dir(STA%SkTot_Ind(i,1):STA%SkTot_Ind(i,2))) * &
                           STA%xStep(i)
        end do

    end subroutine rebuild_corrL

end module statistics_RF
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
