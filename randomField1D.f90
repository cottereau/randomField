module randomField1D

    use displayCarvalhol
    use spectra_RF

contains

!----------------------------------------------------------------------------------------------------------------------
    subroutine createRandomField1D (Nmc, randInit, xMax, kMax, xNStep, kNStep, corrL, corrMod, randField)
        implicit none

        !INPUT
        integer, intent(in) :: Nmc;
        logical, intent(in) :: randInit;
        double precision, dimension(:), allocatable, intent(in) :: xMax, kMax, corrL;
        integer, dimension(:), allocatable, intent(in) :: xNStep, kNStep;
        character (len=20), intent(in) :: corrMod;

        !OUTPUT
        double precision, dimension(:, :), allocatable, intent(out) :: randField;

        !LOCAL VARIABLES
        integer :: i, j;        !Counters
        double precision :: zero = 0.0;
        double precision, dimension(:,:), allocatable :: xVec, kVec;
        double precision :: pi = 3.1415926535898;

        write(*,*) "------------START randomField1D-----------------------------------------";

        if(randInit) then
            call random_seed() !Used for the random field generation
        endif

        allocate(randField(int(maxval(xNStep)),Nmc))
        randField = 0;
        allocate(xVec(int(maxval(xNStep)),Nmc))
        xVec = 0;
        allocate(kVec(int(maxval(kNStep)),Nmc))
        kVec = 0;

        do i=1, Nmc
            if(randInit.eqv..FALSE.) then
                call random_seed() !Reinitialize rand in each step
            endif

            xVec(:,i) = createLinearProg(zero, xMax(i), xNStep(i)) !Creating xVec
            kVec(:,i) = createLinearProg(zero, kMax(i), kNStep(i)) !Creating kVec
            randField(:,i) = make_RandField1D(xVec(:,i), kVec(:,i), corrL(i), corrMod);            !Random Field

            call Disp2D(xVec, "xVec")
            call Disp2D(kVec, "kVec")

        enddo

		deallocate(kVec);
		deallocate(xVec);
        write(*,*) "------------END randomField1D-----------------------------------------";
    end



!----------------------------------------------------------------------------------------------------------------------
    function createLinearProg(qmin, qmax, nSteps) result (lpVec)

        implicit none

        !INPUT
        double precision, intent(in) :: qmin, qmax;
        integer, intent(in) :: nSteps;

        !OUTPUT
        double precision, dimension(:), allocatable :: lpVec;

        !LOCAL VARIABLES
        integer :: i;

        allocate (lpVec(nSteps))

        lpVec = (/((i-1)*(qmax-qmin)/((DBLE(nSteps))-1.0), i=1, nSteps)/)

    end



!----------------------------------------------------------------------------------------------------------------------
    function make_RandField1D(xVec, kVec, corrL, corrMod) result (randField)

        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: xVec, kVec;
        double precision, intent(in) :: corrL;
        character (len=20), intent(in) :: corrMod;

        !OUTPUT
		double precision, dimension(:), allocatable :: randField;

		!LOCAL VARIABLES
        integer :: i, j;        !Counters
        double precision :: kStepDelta;
        double precision, dimension(:,:), allocatable :: phiN, angMat;
        double precision, dimension(:), allocatable :: spectrum;
        double precision :: pi = 3.1415926535898

		write(*,*) "------------IN make_RandField1D-----------------------";

        kStepDelta = kVec(2)-kVec(1)

		!Allocation
		allocate (phiN(size(kVec), size(xVec)));
		allocate (angMat(size(kVec), size(xVec)));
		allocate (randField(size(xVec)));

        !Random Phase (phiN)
        phiN = 0;
        call random_number(phiN(:,1))
        phiN(:,1) = 2*pi*phiN(:,1)
        do i = 2, size(phiN,2)
            phiN(:,i) = phiN(:,1)
        enddo

        !Angle Matrix (angMat)
        angMat = 0
        do j = 1, size(kVec)
            angMat(j, :) = kVec(j)*xVec(:)
        enddo

        !Spectrum creation
        spectrum = make_Spectrum1D(kVec, corrMod, corrL)

        !Random field generation!
        randField = 0;
		randField = matmul(sqrt(2*kStepDelta*spectrum), cos(angMat+phiN));

		call Disp2D(phiN, "phiN")
		call Disp2D(angMat, "angMat")
        call Disp1D(spectrum, "spectrum")

		deallocate (phiN);
		deallocate (angMat);
		deallocate (spectrum);

		write(*,*) "------------OUT make_RandField1D-----------------------";
    end

end module randomField1D
