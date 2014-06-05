program mainCA1D

    use displayCarvalhol
    use vectorOpCarvalhol
    use randField1D
    use readRF
    !regis LIBS = -L/pathtolibs -lblas

!    include "blas.f90"
!    write (*,*) "mainCA Test"

    implicit none

    !Pre-INPUTS (with initialization)
    integer, parameter :: Nmc=1; !number of Monte-Carlo experiments
    logical :: randInit=.TRUE.; !if it's set to false each event will have the same random numbers
    integer :: printOpt= 0; !0-nothing, 1-a little, 2, everything
    integer :: fileOpt= 1; !0-standard, 1-For Matlab

    !INPUTS (without initialization)
    double precision, dimension(Nmc) :: xMax;
    double precision, dimension(Nmc) :: kMax;
    double precision, dimension(Nmc) :: corrL; !correlation length
    character (len=20) :: corrMod
    integer, dimension(Nmc) :: xNStep;
    integer, dimension(Nmc) :: kNStep;

    !OTHERS
    integer :: i, j, file; !Counters
    double precision :: zero = 0.0;
    double precision, dimension(:), allocatable :: locAverage;
    double precision, dimension(:,:), allocatable :: xVec, kVec,randField;
    double precision :: pi = 3.1415926535898;

    write(*,*) "------------START mainCA1D-----------------------------------------";
    write(*,*) "------------Inside mainCA1D > Variables Init-----------------------";

    !Variables Init
    xMax   = (/(2*pi    , i=1, Nmc)/);
    xNStep = (/(100     , i=1, Nmc)/);
    kMax   = (/(2*pi    , i=1, Nmc)/);
    kNStep = (/(100     , i=1, Nmc)/);
    corrL  = (/(1       , i=1, Nmc)/); !correlation length
    corrMod = "gaussian";

    if(randInit.eqv..TRUE.) then
        call random_seed() !Used for the random field generation
    endif

    allocate(randField(int(maxval(xNStep)),Nmc))
    randField = 0;
    allocate(xVec(int(maxval(xNStep)),Nmc))
    xVec = 0;
    allocate(kVec(int(maxval(kNStep)),Nmc))
    kVec = 0;

    write(*,*) "------------Inside mainCA1D > Starting Random Field Creation-----------------------";

    do i=1, Nmc
        if(randInit.eqv..FALSE.) then
            call random_seed() !Reinitialize rand in each step
        endif

        xVec(:,i) = createLinearProg(zero, xMax(i), xNStep(i)) !Creating xVec
        kVec(:,i) = createLinearProg(zero, kMax(i), kNStep(i)) !Creating kVec
        randField(:,i) = createRandField1D(xVec(:,i), kVec(:,i), corrL(i), corrMod); !Random Field

        if(printOpt.gt.0) then
            write(*,*) "Nmc = ", i;
            if(printOpt.eq.2) then
                write(*,*) "xVec";
                call Disp1D(xVec(:,i))
                write(*,*) "kVec";
                call Disp1D(kVec(:,i))
            endif
        endif
    enddo

    if(printOpt.gt.2) then
        write(*,*) "RandomFields";
        call Disp2D(randField)
    endif


    !Writing Result File
    write(*,*) "------------Inside mainCA1D > Writing result file-----------------------";
    file=1;
    open (unit = file , file = "result1Dv2.csv", action = 'write')
        if (fileOpt.eq.0) then
            do i = 1, Nmc
                write (file,'(A15)', advance='NO') "xMax"
                write (file,'(F15.7)', advance='NO') xMax(i)
            enddo
            write (file,*) "" !Jumps Line
            do i = 1, Nmc
                write (file,'(A15)', advance='NO') "xNStep"
                write (file,'(I15)', advance='NO') xNStep(i)
            enddo
            write (file,*) "" !Jumps Line
            do i = 1, Nmc
                write (file,'(A15)', advance='NO') "kMax"
                write (file,'(F15.7)', advance='NO') kMax(i)
            enddo
            write (file,*) "" !Jumps Line
            do i = 1, Nmc
                write (file,'(A15)', advance='NO') "kNStep"
                write (file,'(I15)', advance='NO') kNStep(i)
            enddo
            write (file,*) "" !Jumps Line
            do i = 1, Nmc
                write (file,'(A15)', advance='NO') "corrL"
                write (file,'(F15.7)', advance='NO') corrL(i)
            enddo
            write (file,*) "" !Jumps Line
            do i = 1, Nmc
                write (file,'(A15)', advance='NO') "corrMod"
                write (file,'(A15)', advance='NO') corrMod
            enddo
            write (file,*) "" !Jumps Line
            write (file,*) "" !Jumps Line

            !Titles
            do i = 1, 2*Nmc
                if (mod(i,2).eq.1) then
                    write (file,'(A15)', advance='NO') "xVec"
               else
                    write (file,'(A15)', advance='NO') "randField"
               endif
            enddo
            write (file,*) "" !Jumps Line
        endif
        !Values
        do j = 1, int(maxval(xNStep))
            do i = 1, Nmc
                if (j.gt.xNStep(i)) then
                  write (file,'(F15.7)', advance='NO') xVec(xNStep(i),i)
                  write (file,'(F15.7)', advance='NO') randField(xNStep(i),i)
                else
                  write (file,'(F15.7)', advance='NO') xVec(j,i)
                  write (file,'(F15.7)', advance='NO') randField(j,i)
                endif
            enddo
            write (file,*) "" !Jumps Line
        enddo
    close(file)

    !Statistics
    write(*,*) "------------Inside mainCA1D > Statistics-----------------------";
    allocate(locAverage(Nmc))
    do i = 1, Nmc
        locAverage(i) = averageValue(randField(:,i))
    enddo

    if(printOpt.gt.0) then
        write(*,*) "Local Mean";
        call Disp1D(locAverage)
        write(*,*) "Global Mean";
        write(*,*) averageValue(locAverage);
    endif

    call readfile()

    write(*,*) "------------Inside mainCA1D > END-----------------------";
contains

end program mainCA1D
