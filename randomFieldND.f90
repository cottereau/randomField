module randomFieldND

    use displayCarvalhol
    use spectra_RF
    use math_RF

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

		write(*,*) "";
        write(*,*) "------------START randomFieldND-----------------------------------------";
		write(*,*) "";

		nDim = size(xMax,1);

		!Allocating
		allocate(kMax       (nDim, Nmc));
		allocate(xNStep     (nDim, Nmc));
		allocate(kNStep     (nDim, Nmc));
		allocate(xMin       (nDim, Nmc));
		allocate(kMin       (nDim, Nmc));
		allocate(xNStepTotal(      Nmc));
		allocate(kNStepTotal(      Nmc));

!		write(*,*) "nDim = ", nDim;
!		write(*,*) "corrMod", corrMod;
!		call Disp2D(corrL   , "corrL");

		write(*,*) ">>>>>>>>> Variables initialization: kMax, xNStep and kNStep";
		!Initializing
		call set_kMaxND(corrMod, corrL, kMax)
		xNStep = ceiling(xMax/(pi/kMax))+1;
		kNStep = ceiling(kMax/(2.0d0*pi/xPeriod))+1;
		xMin = 0d0;
		kMin = 0d0;

!		call Disp2D(xMin, "xMin");
!		call Disp2D(xMax, "xMax");
!		call Disp2D(kMin, "kMin");
!		call Disp2D(kMax, "kMax");
!		write(*,*) "xNStep = ", xNStep*1
!		write(*,*) "kNStep = ", kNStep*1

		!ONLY FOR TESTS (Overwriting)----------------------
		kMax  (1, :) = (/(2*pi, i=1, Nmc)/);
	    xNStep(1, :) = (/(50,  i=1, Nmc)/);
	    kNStep(1, :) = (/(50,  i=1, Nmc)/);

		if(nDim > 1) then
			kMax  (2, :) = (/(2*pi, i=1, Nmc)/);
		    xNStep(2, :) = (/(50,    i=1, Nmc)/);
			kNStep(2, :) = (/(50,    i=1, Nmc)/);
		end if

		if(nDim > 2) then
			kMax  (3, :) = (/(2*pi, i=1, Nmc)/);
		    xNStep(3, :) = (/(50,    i=1, Nmc)/);
			kNStep(3, :) = (/(50,    i=1, Nmc)/);
		end if
		!---------------------------------------------------

		!Working on Total Number of Permutations
		do i = 1, size(xNStep,2)
			xNStepTotal(i) = product(xNStep(:,i));
		end do

		do i = 1, size(kNStep,2)
			kNStepTotal(i) = product(kNStep(:,i));
		end do

!		call Disp1Dint(xNStepTotal, "xNStepTotal");
!		call Disp1Dint(kNStepTotal, "kNStepTotal");



		!Random Field
		allocate(randField(int(maxval(xNStepTotal)),Nmc));
		allocate(randLoc  (int(maxval(xNStepTotal))    ));
		allocate(deltaK   (nDim                        ));
		randField = 0;
		randLoc   = 0;
		deltaK    = 0;

		write(*,*) "";
		write(*,*) ">>>>>>>>> Random Field Creation";
		write(*,*) "";
		write(*,*) "Number of dimensions = ", nDim;
		write(*,*) "Number of events     = ", Nmc;
		write(*,*) "Number x points      = ", maxval(xNStepTotal);
		write(*,*) "Number k points      = ", maxval(kNStepTotal);
		write(*,*) "";


		do k = 1, Nmc
            if(.not.randInit) call random_seed() !Reinitialize rand for each Monte-Carlo event
		    deltaK = (kMax(:,k)-kMin(:,k))/(kNStep(:,k)-1); !Defines deltaK
		    kNS    = kNStepTotal(k);
		    xNS    = xNStepTotal(k);
		    !call disp1D(deltaK, "deltaK")

			do j = 1, kNS
				call random_number(phiN)
				!write(*,*) "phiN", phiN
			    kVec           = get_Permutation(j, kMin(:,k), kMax(:,k), kNStep(:,k));
				Sk             = get_SpectrumND(kVec, corrMod(k), corrL(:,k));
				randLoc        = 0;
				randLoc(1:xNS) = 2*pi*phiN !Random part of the angle matrix for this k permutation

				do i = 1, xNS
					xVec       = get_Permutation(i, xMin(:,k), xMax(:,k), xNStep(:,k));
					randLoc(i) = randLoc(i) + dot_product(kVec, xVec); !Not-random part of the angle matrix
				end do

				randLoc  (1:xNS)   = (2*Sk*product(deltaK))**(0.5d0)*cos(randLoc(1:xNS))

				randField(1:xNS,k) = randField(1:xNS,k) + randLoc(1:xNS)
			end do
			write(*,*) "Event ",k, "of", Nmc, "completed"
		end do

		randField(:,:) = 2d0**(0.5d0+(nDim-1.0d0))*randField(:,:)

		!Only printing------------------------------------------
!		call Disp2Dint(xNStep, "xNStep");
!		call Disp2Dint(kNStep, "kNStep");
!
!		write(*,*) "Permutation X (Nmc = 1)"
!		do i = 1, xNStepTotal(1)
!			write(*,*) i, ">", get_Permutation(i, xMin(:,1), xMax(:,1), xNStep(:,1));
!		end do
!
!		write(*,*) "Permutation K (Nmc = 1)"
!		do i = 1, kNStepTotal(1)
!			write(*,*) i, ">", get_Permutation(i, kMin(:,1), kMax(:,1), kNStep(:,1));
!		end do
!
!		!Spectrum
!		write(*,*) "Spectrum (Nmc = 1)"
!		do i = 1, kNStepTotal(1)
!			write(*,*) get_SpectrumND(get_Permutation(i, kMin(:,1), kMax(:,1), kNStep(:,1)), &
!			                          corrMod(1),                                               &
!			                          corrL(:,1));
!		end do

		!---------------------------------------------------

		deallocate(xMin);
		deallocate(kMin);
		deallocate(xNStepTotal);
		deallocate(kNStepTotal);
		deallocate(deltaK);

		write(*,*) "";
        write(*,*) "------------END randomFieldND-----------------------------------------";
		write(*,*) "";
    end

end module randomFieldND
