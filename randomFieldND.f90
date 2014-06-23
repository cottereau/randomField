module randomFieldND

    use displayCarvalhol
    use spectra_RF

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine createRandomFieldND(Nmc, randInit, corrMod, corrL, xMax, xPeriod, &
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
        integer :: i, j, k, nDim, kNS, xNS;
        double precision, dimension(:,:), allocatable :: xMin, kMin;
        integer,          dimension(:),   allocatable :: xNStepTotal, kNStepTotal;
        double precision                              :: Sk, phiN;
        double precision, dimension(:),   allocatable :: xVec, kVec; !Allocated in function
        double precision, dimension(:),   allocatable :: deltaK, randLoc;
        double precision                              :: pi = 3.1415926535898, zero = 0d0;

        write(*,*) "------------START randomFieldND-----------------------------------------";

		nDim = size(xMax,1);

		!Allocating
		allocate(kMax       (nDim, Nmc));
		allocate(xNStep     (nDim, Nmc));
		allocate(kNStep     (nDim, Nmc));
		allocate(xMin       (nDim, Nmc));
		allocate(kMin       (nDim, Nmc));
		allocate(xNStepTotal(      Nmc));
		allocate(kNStepTotal(      Nmc));

		write(*,*) "nDim = ", nDim;
		write(*,*) "corrMod", corrMod;
		call Disp2D(corrL   , "corrL");

		!Initializing
		call set_kMaxND(corrMod, corrL, kMax)
		xNStep = ceiling(xMax/(pi/kMax))+1;
		kNStep = ceiling(kMax/(2.0d0*pi/xPeriod))+1;
		xMin = 0d0;
		kMin = 0d0;

		call Disp2D(xMin, "xMin");
		call Disp2D(xMax, "xMax");
		call Disp2D(kMin, "kMin");
		call Disp2D(kMax, "kMax");

		!ONLY FOR TESTS (Overwriting)
		kMax  (1, :) = (/(2*pi,  i=1, Nmc)/);
	    xNStep(1, :) = (/(10*i,  i=1, Nmc)/);
	    kNStep(1, :) = (/(10*i,  i=1, Nmc)/);

!		kMax  (2, :) = (/(2*pi, i=1, Nmc)/);
!	    xNStep(2, :) = (/(4,    i=1, Nmc)/);
!		kNStep(2, :) = (/(4,    i=1, Nmc)/);
!
!		kMax  (3, :) = (/(2*pi, i=1, Nmc)/);
!	    xNStep(3, :) = (/(3,    i=1, Nmc)/);
!		kNStep(3, :) = (/(2,    i=1, Nmc)/);

		call Disp2Dint(xNStep, "xNStep");
		call Disp2Dint(kNStep, "kNStep");

		!Working on Permutations
		call set_stepsTotal(xNStep, xNStepTotal)
		call set_stepsTotal(kNStep, kNStepTotal);

		call Disp1Dint(xNStepTotal, "xNStepTotal");
		call Disp1Dint(kNStepTotal, "kNStepTotal");

		write(*,*) "Permutation X (Nmc = 1)"
		do i = 1, xNStepTotal(1)
			write(*,*) i, ">", get_Permutation(i, xMin(:,1), xMax(:,1), xNStep(:,1));
		end do

		write(*,*) "Permutation K (Nmc = 1)"
		do i = 1, kNStepTotal(1)
			write(*,*) i, ">", get_Permutation(i, kMin(:,1), kMax(:,1), kNStep(:,1));
		end do

		!Spectrum
		write(*,*) "Spectrum (Nmc = 1)"
		do i = 1, kNStepTotal(1)
			write(*,*) get_SpectrumND(get_Permutation(i, kMin(:,1), kMax(:,1), kNStep(:,1)), &
			                          corrMod(1),                                               &
			                          corrL(:,1));
		end do

		!Random Field
		allocate(randField(int(maxval(xNStepTotal)),Nmc));
		allocate(randLoc  (int(maxval(xNStepTotal))    ));
		allocate(deltaK   (nDim                        ));
		randField = 0;
		randLoc   = 0;
		deltaK    = 0;

		do k = 1, Nmc
			write(*,*) "k = ",k
            if(.not.randInit) call random_seed() !Reinitialize rand for each Monte-Carlo event
		    deltaK = (kMax(:,k)-kMin(:,k))/(kNStep(:,k)-1); !Defines deltaK
		    kNS    = kNStepTotal(k);
		    xNS    = xNStepTotal(k);
		    call disp1D(deltaK, "deltaK")

			do j = 1, kNS
				call random_number(phiN)
				write(*,*) "phiN", phiN
			    kVec           = get_Permutation(j, kMin(:,k), kMax(:,k), kNStep(:,k));
				Sk             = get_SpectrumND(kVec, corrMod(k), corrL(:,k));
				randLoc        = 0;
				randLoc(1:xNS) = 2*pi*phiN !Random part of the angle matrix for this k permutation

!				write(*,*) "j = ",j
!				write(*,*) "Sk = ",Sk
!				call disp1D(kVec,"kVec")
!				call disp1D(randLoc(1:5),"randLoc(1:5) (only random part)")

				do i = 1, xNS
					xVec       = get_Permutation(i, xMin(:,k), xMax(:,k), xNStep(:,k));
					randLoc(i) = randLoc(i) + dot_product(kVec, xVec); !Not-random part of the angle matrix

					!write(*,*) "i = ",i
				end do
!				call disp1D(randLoc(1:5),"randLoc(1:5) (after non-random part)")
				randLoc  (1:xNS)   = (2*Sk*product(deltaK))**(0.5d0)*cos(randLoc(1:xNS))
!				call disp1D(randLoc(1:5),"randLoc(1:5) (contribution to randField)")
				randField(1:xNS,k) = randField(1:xNS,k) + randLoc(1:xNS)
!				call disp1D(randField(:,k),"randField")
			end do
		end do

		randField(:,:) = 2d0**(0.5d0+(nDim-1.0d0))*randField(:,:)

		deallocate(xNStepTotal);
		deallocate(kNStepTotal);
		deallocate(deltaK);

        write(*,*) "------------END randomFieldND-----------------------------------------";
    end
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine set_stepsTotal(NStep, NStepTotal);

        implicit none

        !INPUT
        integer, dimension(:, :), intent(in) :: NStep;
        !OUTPUT
        integer, dimension(:), intent(out) :: NStepTotal;
        !LOCAL VARIABLES
        integer :: i;

		do i = 1, size(NStep,2)
			NStepTotal(i) = product(NStep(:,i));
		end do

    end

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    function get_Permutation(pos, qmin, qmax, nStep) result(pVec)

        implicit none

        !INPUT
        integer                                    :: pos;
        double precision, dimension(:), intent(in) :: qmin, qmax;
        integer,          dimension(:), intent(in) :: nStep;
        !OUTPUT
        double precision, dimension(:), allocatable :: pVec;
        !LOCAL VARIABLES
        integer :: i, j;
        integer :: seedStep, nDim;

		nDim = size(nStep);
		allocate(pVec(nDim))

		do j = 1, nDim
		    seedStep = product(nStep(j+1:));
		    if (j == nDim) seedStep = 1;
			i = cyclicMod(int((pos-0.1)/seedStep)+1, nStep(j))
			pVec(j) = (i-1)*(qmax(j)-qmin(j))/(nStep(j)-1);
		end do
    end

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    function cyclicMod(pos, base) result(resPos)

        implicit none
        !INPUT
        integer, intent(in) :: pos, base; !desired position
        !OUTPUT
        integer :: resPos;!result Position

		resPos = mod(pos,base);
		if(resPos == 0) resPos = base;

    end
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!    function make_RandField2D(xVec, kVec, corrL, corrMod) result (randField)
!
!        implicit none
!
!        !INPUT
!        double precision, dimension(:), intent(in) :: xVec, kVec;
!        double precision, intent(in) :: corrL;
!        character (len=*), intent(in) :: corrMod;
!
!        !OUTPUT
!		double precision, dimension(:), allocatable :: randField;
!
!		!LOCAL VARIABLES
!        integer :: i, j;        !Counters
!        double precision :: kStepDelta;
!        double precision, dimension(:,:), allocatable :: phiN, angMat;
!        double precision, dimension(:), allocatable :: spectrum;
!        double precision :: pi = 3.1415926535898
!
!		write(*,*) "------------IN make_RandField1D-----------------------";
!
!        kStepDelta = kVec(2)-kVec(1)
!
!		!Allocation
!		allocate (phiN(size(kVec), size(xVec)));
!		allocate (angMat(size(kVec), size(xVec)));
!		allocate (randField(size(xVec)));
!
!        !Random Phase (phiN)
!        phiN = 0;
!        call random_number(phiN(:,1))
!        phiN(:,1) = 2*pi*phiN(:,1)
!        do i = 2, size(phiN,2)
!            phiN(:,i) = phiN(:,1)
!        enddo
!
!        !Angle Matrix (angMat)
!        angMat = 0
!        do j = 1, size(kVec)
!            angMat(j, :) = kVec(j)*xVec(:)
!        enddo
!
!        !Spectrum creation
!        spectrum = make_Spectrum2D(kVec, corrMod, corrL)
!
!        !Random field generation!
!        randField = 0;
!		randField = matmul(sqrt(2*kStepDelta*spectrum), cos(angMat+phiN));
!
!		call Disp2D(phiN, "phiN")
!		call Disp2D(angMat, "angMat")
!        call Disp1D(spectrum, "spectrum")
!
!		deallocate (phiN);
!		deallocate (angMat);
!		deallocate (spectrum);
!
!		write(*,*) "------------OUT make_RandField1D-----------------------";
!    end

end module randomFieldND
