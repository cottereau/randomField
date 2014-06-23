module randomField1D

    use displayCarvalhol
    use spectra_RF

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine createRandomField1D (Nmc, randInit, corrMod, corrL, xMax, xPeriod, &
    								kMax, xNStep, kNStep, randField);

        implicit none

        !INPUT
        integer,                             intent(in) :: Nmc;
        logical,                             intent(in) :: randInit;
        character (len=*),  dimension(:),    intent(in) :: corrMod;
        double precision,   dimension(:, :), intent(in) :: corrL, xMax, xPeriod;

        !OUTPUT
        double precision, dimension(:, :), allocatable, intent(out) :: kMax;
        integer,          dimension(:, :), allocatable, intent(out) :: xNStep, kNStep;
        double precision, dimension(:, :), allocatable, intent(out) :: randField;

        !LOCAL VARIABLES
        integer :: i, j, nDim;
        double precision :: zero = 0.0;
        double precision, dimension(:,:), allocatable :: xVec, kVec;
        double precision :: pi = 3.1415926535898;

        write(*,*) "------------START randomField1D-----------------------------------------";

		nDim = 1;

		!Allocating
		allocate(kMax  (Nmc, nDim));
		allocate(xNStep(Nmc, nDim));
		allocate(kNStep(Nmc, nDim));

		!Calculating kMax, kNStep and xNStep
		call set_kMax1D(corrMod, corrL, kMax)
		xNStep = ceiling(xMax/(pi/kMax))+1;
		kNStep = ceiling(kMax/(2.0d0*pi/xPeriod))+1;

		!ONLY FOR TESTS (Overwriting)
		kMax  (:, 1) = (/(2*pi,  i=1, Nmc)/);
	    xNStep(:, 1) = (/(10*i,  i=1, Nmc)/);
	    kNStep(:, 1) = (/(10*i,    i=1, Nmc)/);

		call Disp2D(corrL, "corrL");
		call Disp2D(xMax, "xMax");
		call Disp2D(xPeriod, "xPeriod");
		call Disp2D(kMax, "kMax");
		write(*,*) "xNStep", xNStep;
		write(*,*) "kNStep", kNStep;

		!Calculating
        allocate(randField(int(maxval(xNStep(:,1))),Nmc))
        randField = 0;
        allocate(xVec(int(maxval(xNStep(:,1))),Nmc))
        xVec = 0;
        allocate(kVec(int(maxval(kNStep(:,1))),Nmc))
        kVec = 0;

        if(randInit) then
            call random_seed() !Used for the random field generation
        endif

        do i=1, Nmc
            if(.not.randInit) then
                call random_seed() !Reinitialize rand in each step
            endif

            xVec(:,i) = createLinearProg(zero, xMax(i,1), xNStep(i,1)) !Creating xVec
            kVec(:,i) = createLinearProg(zero, kMax(i,1), kNStep(i,1)) !Creating kVec

			call Disp2D(xVec, "xVec");
			call Disp2D(kVec, "kVec");
			call Disp2D(randField, "randField");
            randField(1:xNStep(i,1),i) = make_RandField1D       &
            							(xVec(1:xNStep(i,1),i), &
            							 kVec(1:kNStep(i,1),i), &
            							 corrL(i,1),            &
            							 corrMod(i));            !Random Field
            call Disp2D(randField, "randField");
        enddo

		deallocate(kVec);
		deallocate(xVec);
        write(*,*) "------------END randomField1D-----------------------------------------";
    end



!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    function createLinearProg(qmin, qmax, nSteps) result (lpVec)

        implicit none

        !INPUT
        double precision, intent(in) :: qmin, qmax;
        integer,          intent(in) :: nSteps;

        !OUTPUT
        double precision, dimension(:), allocatable :: lpVec;

        !LOCAL VARIABLES
        integer :: i;

        allocate (lpVec(nSteps))

        lpVec = (/((i-1)*(qmax-qmin)/((DBLE(nSteps))-1.0), i=1, nSteps)/)

    end



!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    function make_RandField1D(xVec, kVec, corrL, corrMod) result (randField)

        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: xVec, kVec;
        double precision,               intent(in) :: corrL;
        character (len=*),              intent(in) :: corrMod;

        !OUTPUT
		double precision, dimension(:), allocatable :: randField;

		!LOCAL VARIABLES
        integer                                       :: i, j;
        double precision                              :: kStepDelta;
        double precision, dimension(:,:), allocatable :: phiN, angMat;
        double precision, dimension(:), allocatable   :: spectrum;
        double precision                              :: pi = 3.1415926535898;

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
		randField = sqrt(2.0d0)*matmul(sqrt(2*kStepDelta*spectrum), cos(angMat+phiN));

		call Disp2D(phiN, "phiN")
		call Disp2D(angMat, "angMat")
        call Disp1D(spectrum, "spectrum")

		deallocate (phiN);
		deallocate (angMat);
		deallocate (spectrum);

		write(*,*) "------------OUT make_RandField1D-----------------------";
    end

end module randomField1D
