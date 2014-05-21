module randField1D

    use displayCarvalhol
    use vectorOpCarvalhol
    !regis LIBS = -L/pathtolibs -lblas

!    include "blas.f90"
!    write (*,*) "mainCA Test"

contains

    function createRandField1D(xVec, kVec, corrL, corrMod) result (randField)
        implicit none
        double precision, dimension(:), intent(in) :: xVec, kVec;
        double precision, intent(in) :: corrL;
        character (len=20), intent(in) :: corrMod;

        integer :: i, j;        !Counters
        double precision :: kStepDelta;
        double precision, dimension(:,:), allocatable :: phiN, angMat;
        double precision, dimension(:), allocatable :: Sk, randField;
        double precision :: pi = 3.1415926535898
        integer :: printOpt= 1; !0-nothing, 1-a little, 2, everything

        if (printOpt.gt.0) then
            write(*,*) "------------Inside randField1D > createRandField1D-----------------------";
        endif

        kStepDelta = kVec(2)-kVec(1)

        !Random Phase (phiN)
        allocate (phiN(size(kVec), size(xVec)));
        phiN = 0;
        call random_number(phiN(:,1))
        phiN(:,1) = 2*pi*phiN(:,1)
        do i = 2, size(phiN,2)
            phiN(:,i) = phiN(:,1)
        enddo
    
        !Angle Matrix (angMat)
        allocate (angMat(size(kVec), size(xVec)));
        !call gemm(kVec,xVec, angMat)
        angMat = 0
        do i = 1, size(xVec) !Matrix Multiplication (a changer par "gemm")
            do j = 1, size(kVec)
                    angMat(j, i) = angMat(j, i) + kVec(j)*xVec(i)
            enddo
        enddo

    
        !Spectrum creation
        Sk = createSpectrum(kVec, corrMod, corrL)
    
        !Random field generation
        allocate(randField(size(xVec)))
        randField = 0;
        !call gemm(sqrt(2*kDeltaFac*Sk, cos(angMat), randField)
        do i = 1, size(angMat, 2) !Matrix Multiplication (a changer par "gemm")
            do j = 1, size(angMat, 1)
                randField(i) = randField(i) + sqrt(2*kStepDelta*Sk(j)) * cos(angMat(j,i)+phiN(j,i))
            enddo
        enddo
        
        if (printOpt.gt.1) then
            write(*,*) "phiN";
            call Disp2D(phiN);
            write(*,*) "angMat";
            call Disp2D(angMat);
            write(*,*) "Sk";
            call Disp1D(Sk);
            write(*,*) "randField";
            call Disp1D(randField);
        endif

        if (printOpt.gt.0) then
            write(*,*) "------------End randField1D > createRandField1D-----------------------";
        endif

    end

    function createSpectrum(kVector, corrMod, corrL) result (Sk)
!        ! Return Spectrum from a chosen correlation model
        implicit none
        integer :: j;
        double precision :: eta, corrL;
        double precision, dimension(:) :: kVector;
        double precision, dimension(:), allocatable :: Sk;
        character (len=20) :: corrMod
        double precision :: pi = 3.1415926535898

        allocate (Sk(size(kVector)));

        !write(*,*) "corrMod = ", corrMod;

        select case(corrMod)

            case("gaussian")
                do j = 1, size(kVector)
                    !write(*,*) "Ski = ", Sk(j);
                    eta = corrL/sqrt(pi)
                    Sk(j) = ((1./(2.*pi)) * exp((-(eta**2)/4)*(kVector(j)**2)));
                enddo

        end select
    end

end module randField1D
