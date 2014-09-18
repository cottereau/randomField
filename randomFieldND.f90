module randomFieldND

    use displayCarvalhol
    use spectra_RF
    use math_RF
    use obsolete_RF
	use mpi

    interface createRandomField
		module procedure createRandomFieldUnstruct,   &
		                 createRandomFieldStructured
	end interface createRandomField

	double precision :: periodMult = 1.1 !"range" multiplier
	double precision :: kAdjust = 10 !"kNStep minimum" multiplier

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine createRandomFieldUnstruct(xPoints, corrMod, corrL, fieldAvg, fieldVar, Nmc, randField);
        implicit none

        !INPUT
        double precision, dimension(:, :), intent(in) :: xPoints;
        double precision, dimension(:)   , intent(in) :: corrL;
        character (len=*)                , intent(in) :: corrMod;
        integer                          , intent(in) :: Nmc;
        double precision                 , intent(in) :: fieldAvg, fieldVar;

        !OUTPUT
        double precision, dimension(:, :), allocatable, intent(out) :: randField;

        !LOCAL VARIABLES
        integer         , dimension(:)  , allocatable :: kNStep;
        double precision, dimension(:)  , allocatable :: kMax;
        double precision, dimension(:,:), allocatable :: kSign, phiN;
        double precision, dimension(:)  , allocatable :: kVec, kVecUnsigned; !Allocated in function
        double precision, dimension(:)  , allocatable :: deltaK, angleVec, kDelta;
        double precision, dimension(:)  , allocatable :: xMaxGlob, xMinGlob;
        integer          :: i, j, k, m, nDim;
        integer          :: xNTotal, kNTotal;
        integer          :: nb_procs, rang, code, error;
        double precision :: Sk, randPhase, deltaKprod, normalAvg, normalVar;
        !double precision :: periodMult = 1.1, kAdjust = 5;
        double precision :: pi = 3.1415926535898, zero = 0d0;

		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
		call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)

		!if(rang == 0) write(*,*) "";
        !if(rang == 0) write(*,*) "------------START randomFieldND (proc = ", rang, ")---------";
		!if(rang == 0) write(*,*) "";

		nDim    = size(xPoints, 2);
		xNTotal = size(xPoints, 1);

		!if(rang == 0) write(*,*) ">>>>>>>>> Variables initialization: kMax, kDelta, kNStep, xMinGlob, xMaxGlob";
		!Allocating
		allocate(kMax   (nDim));
		allocate(kNStep (nDim));
		allocate(kDelta (nDim));


		call set_Extremes(xPoints, xMinGlob, xMaxGlob) !Communicating extremes
		call set_kMaxND(corrMod, corrL, kMax) !Defining kMax according to corrMod and corrL
		kDelta     = 2*pi/(periodMult*(xMaxGlob - xMinGlob)) !Delta max in between two wave numbers to avoid periodicity
		kNStep     = kAdjust*(ceiling(kMax/kDelta) + 1);
		kNTotal    = product(kNStep);



		!Random Field
		allocate(angleVec    (xNTotal));
		allocate(deltaK      (nDim));
		allocate(kVec        (nDim));
		allocate(kVecUnsigned(nDim));
		allocate(kSign       (2**(nDim-1), nDim));
		allocate(phiN        (size(kSign,1), kNTotal));
		allocate(randField   ((xNTotal),Nmc));

		if(rang == 0) then
			write(*,*) "";
			write(*,*) ">>>>>>>>> Random Field Creation (only showing proc 0)";
			write(*,*) "";
			write(*,*) "Number of dimensions = ", nDim;
			write(*,*) "Number of events     = ", Nmc;
			write(*,*) "Number x points ";
			write(*,*) "      by Dimension   = --- (unstructured)";
			write(*,*) "      Total          = ", xNTotal;
			write(*,*) "Number k points ";
			write(*,*) "      by Dimension   = ", kNStep;
			write(*,*) "      Total          = ", kNTotal;
			write(*,*) "";
		end if

		randField  = 0;
		angleVec   = 0;
		deltaK     = 0;
	    deltaK     = (kMax)/(kNStep-1); !Defines deltaK
		call set_kSign(kSign)

		select case (corrMod)
			case("gaussian")
				normalVar = fieldVar
				normalAvg = fieldAvg
			case("lognormal")
				if(fieldAvg <= 0) then
					write(*,*) ""
					write(*,*) "ERROR - when using lognormal fieldAvg should be a positive number greater than 0.001"
					call MPI_ABORT(MPI_COMM_WORLD, error, code)
				end if
				normalVar = log(1 + fieldVar/(fieldAvg**2))
				normalAvg = log(fieldAvg) - normalVar/2

		end select

		if(rang == 0) write(*,*) "normalVar = ", normalVar
		if(rang == 0) write(*,*) "normalAvg = ", normalAvg

		do k = 1, Nmc
			call random_number(phiN(:,1:kNTotal))
			do j = 1, kNTotal
				call get_Permutation(j, kMax, kNStep, kVecUnsigned);
			    do m = 1, size(kSign,1)
			    	kVec      = kVecUnsigned * kSign(m, :)
					Sk        = get_SpectrumND(kVec, corrMod, corrL);
					randPhase = 2*pi*phiN(m, j)
					do i = 1, xNTotal
						angleVec(i) = dot_product(kVec, xPoints(i,:)) + randPhase !angle + random phase
					end do
					randField(:,k) = sqrt(Sk) * cos(angleVec(:)) + randField(:,k)
				end do
			end do
			!if(rang == 0) write(*,*) "Event ",k, "of", Nmc, "completed (counting only in proc 0)";
		end do

		randField(:,:) = sqrt(2/pi)*sqrt(product(corrL)*product(deltaK))&
						 * randField(:,:) * sqrt(normalVar)&
						 + normalAvg; !Obs: sqrt(2*pi) division
!		randField(:,:) = (2.0d0)**(2.0d0-nDim)*sqrt(product(corrL)*product(deltaK)/pi)&
!						 * randField(:,:) * sqrt(normalVar)&
!						 + normalAvg;

		if (corrMod == "lognormal") then
			randField(:,:) = exp(randField(:,:))
		end if


		if(allocated(deltaK))       deallocate(deltaK);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
		if(allocated(kVec))         deallocate(kVec);
		if(allocated(kMax))         deallocate(kMax);
		if(allocated(kNStep))       deallocate(kNStep);
		if(allocated(kDelta))       deallocate(kDelta);
		if(allocated(angleVec))     deallocate(angleVec);
		if(allocated(kSign))        deallocate(kSign);
		if(allocated(phiN))         deallocate(phiN);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
		if(allocated(xMinGlob))     deallocate(xMinGlob);
        if(allocated(xMaxGlob))     deallocate(xMaxGlob);


		!if(rang == 0) write(*,*) "";
        !if(rang == 0) write(*,*) "------------END randomFieldND (proc = ", rang, "---------";
		!if(rang == 0) write(*,*) "";

    end subroutine createRandomFieldUnstruct

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine createRandomFieldStructured(xMin, xMax, corrMod, corrL, Nmc, &
										   xNStep, randField);

        implicit none

        !INPUT
        integer                       , intent(in) :: Nmc;
        character (len=*)             , intent(in) :: corrMod;
        double precision, dimension(:), intent(in) :: corrL, xMax, xMin;

        !OUTPUT
        integer         , dimension(:)   , allocatable, intent(out) :: xNStep;
        double precision, dimension(:, :), allocatable, intent(out) :: randField;

        !LOCAL VARIABLES
        integer         , dimension(:)  , allocatable :: kNStep;
        double precision, dimension(:)  , allocatable :: kMax;
        double precision, dimension(:,:), allocatable :: kSign, phiN;
        double precision, dimension(:),   allocatable :: xVec, kVec, kVecUnsigned; !Allocated in function
        double precision, dimension(:),   allocatable :: deltaK, angleVec, kDelta;
        double precision, dimension(:)  , allocatable :: xMaxGlob, xMinGlob;
        integer          :: i, j, k, m, nDim;
        integer          :: xNTotal, kNTotal;
        integer          :: testDim; !Only for tests
        integer          :: nb_procs, rang, code !, xStart, xEnd, sizeUnif, sizeLoc;
        double precision :: xPerCorrL = 5; !number of points per Correlation Length
        double precision :: Sk;
        double precision :: pi = 3.1415926535898, zero = 0d0;

		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
		call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)

		!if(rang == 0) write(*,*) "";
        !if(rang == 0) write(*,*) "------------START randomFieldND (proc = ", rang, ")---------";
		!if(rang == 0) write(*,*) "";

		nDim = size(xMax,1);

		!if(rang == 0) write(*,*) ">>>>>>>>> Variables initialization: kMax, xNStep and kNStep";
		!Allocating
		allocate(kMax   (nDim));
		allocate(xNStep (nDim));
		allocate(kNStep (nDim));
		allocate(kDelta (nDim));

		!Communicating extremes
		call set_Extremes(xMin, xMax, xMinGlob, xMaxGlob)

		!Initializing (obs: parameters FOR EACH PROC)
		call set_kMaxND(corrMod, corrL, kMax)
		xNStep  = xPerCorrL*ceiling((xMax-xMin)/corrL);
		kDelta  = 2*pi/(periodMult*(xMaxGlob - xMinGlob)) !Delta max in between two wave numbers to avoid periodicity
		kNStep  = kAdjust*(ceiling(kMax/kDelta) + 1);
		xNTotal = product(xNStep);
		kNTotal = product(kNStep);

!		if (rang == 0) then
!			call DispCarvalhol (xMax, "xMax")
!			call DispCarvalhol (xMin, "xMin")
!			call DispCarvalhol (corrL, "corrL")
!			call DispCarvalhol ((xMax-xMin), "(xMax-xMin)")
!			call DispCarvalhol (ceiling((xMax-xMin)/corrL), "ceiling((xMax-xMin)/corrL)")
!			call DispCarvalhol (kMax, "kMax")
!			call DispCarvalhol (xPeriod, "xPeriod")
!			call DispCarvalhol (kNStep, "kNStep")
!		end if



		!Random Field
		allocate(randField((xNTotal),Nmc));
		allocate(angleVec   (xNTotal));
		allocate(deltaK    (nDim));
		allocate(kSign     (2**(nDim-1), nDim));
		allocate(phiN      (size(kSign,1), kNTotal));
		allocate(kVecUnsigned(nDim));
		allocate(kVec(nDim));
		allocate(xVec(nDim));

		if(rang == 0) then
			write(*,*) "";
			write(*,*) ">>>>>>>>> Random Field Creation (only showing proc 0)";
			write(*,*) "";
			write(*,*) "Number of dimensions = ", nDim;
			write(*,*) "Number of events     = ", Nmc;
			write(*,*) "Number x points ";
			write(*,*) "      by Dimension   = ", xNStep;
			write(*,*) "      Total          = ", xNTotal;
			write(*,*) "Number k points ";
			write(*,*) "      by Dimension   = ", kNStep;
			write(*,*) "      Total          = ", kNTotal;
			write(*,*) "";
		end if

		randField = 0;
		angleVec   = 0;
		deltaK    = 0;
	    deltaK    = (kMax)/(kNStep-1); !Defines deltaK
		call set_kSign(kSign)

		do k = 1, Nmc
			call random_number(phiN(:,1:kNTotal))
			do j = 1, kNTotal
				call get_Permutation(j, kMax, kNStep, kVecUnsigned);
			    do m = 1, size(kSign,1)
			    	kVec    = kVecUnsigned * kSign(m, :)
					Sk      = get_SpectrumND(kVec, corrMod, corrL);
					angleVec(:) = 0;
					angleVec(:) = 2 * pi * phiN(m, j) !TO VERIFY - Random part of the angle matrix for this k permutation
					do i = 1, xNTotal
						call get_Permutation(i, xMax, xNStep, xVec, xMin);
						angleVec(i) = angleVec(i) + dot_product(kVec, xVec); !Not-random part of the angle matrix
					end do
					angleVec(:)     = sqrt(Sk*2*product(deltaK)) &
										     *cos(angleVec(:))
					randField(:,k) = randField(:,k) + angleVec(:)
				end do
			end do
			randField(:,k) = product(corrL)*randField(:,k)
			!if(rang == 0) write(*,*) "Event ",k, "of", Nmc, "completed (counting only in proc 0)";
		end do

		randField(:,:) = sqrt(2d0)*randField(:,:)

		!Only printing------------------------------------------

		if(rang == 0) then
			write (*,*) "randField = ", rang
			do i = 1, 10! size(randField, 1)
				write(*,*) randField(i,:)
			end do
		end if

!		call DispCarvalhol(xNStep, "xNStep");
!		call DispCarvalhol(kNStep, "kNStep");
!
!
!		if (rang == 0) then
!			write(*,*) "Permutation X (Nmc = 1), rang ", rang
!			do i = 1, xNTotal
!				call get_Permutation(i, xMax, xNStep, xVec, xMin);
!				write(*,'(I,A,3F10.5)') i, ">", xVec;
!			end do
!		end if
!
!		if (rang == 0) then
!			call DispCarvalhol(kMax  , "kMax")
!			call DispCarvalhol(kNStep, "kNStep")
!			write(*,*) "Permutation K (Nmc = 1)"
!			do i = 1, kNTotal
!			call get_Permutation(i, kMax, kNStep, kVecUnsigned);
!				write(*,*) i, ">", kVecUnsigned;
!			end do
!		end if
!
!		!Spectrum
!		write(*,*) "Spectrum (Nmc = 1)"
!		do i = 1, kNTotal
!			call get_Permutation(j, kMax, kNStep, kVecUnsigned);
!		    do m = 1, size(kSign,1)
!		    	kVec           = kVecUnsigned * kSign(m, :)
!				Sk             = get_SpectrumND(kVec, corrMod, corrL);
!			end do
!		end do
!
!		call DispCarvalhol(randField)

		!---------------------------------------------------

		if(allocated(deltaK))       deallocate(deltaK);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
		if(allocated(kVec))         deallocate(kVec);
		if(allocated(xVec))         deallocate(xVec)
		if(allocated(kMax))         deallocate(kMax);
		if(allocated(kNStep))       deallocate(kNStep);
		if(allocated(angleVec))     deallocate(angleVec);
		if(allocated(kSign))        deallocate(kSign);
		if(allocated(phiN))         deallocate(phiN);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
        if(allocated(xMaxGlob))     deallocate(xMaxGlob);
        if(allocated(xMinGlob))     deallocate(xMinGlob);
        if(allocated(kDelta))       deallocate(kDelta);


		!if(rang == 0) write(*,*) "";
        !if(rang == 0) write(*,*) "------------END randomFieldND (proc = ", rang, "---------";
		!if(rang == 0) write(*,*) "";
    end subroutine createRandomFieldStructured

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine set_XPoints(corrL, xMin, xMax, xPoints, nPoinsPerCorrL)
    	implicit none

    	!INPUT
    	double precision, dimension(:), intent(in) :: corrL, xMin, xMax;
    	integer                       , intent(in) :: nPoinsPerCorrL !Number of points by correlation length
    	!OUTPUT
    	double precision, dimension(:,:), allocatable, intent(OUT) :: xPoints;
    	!LOCAL VARIABLES
    	integer :: nDim, i, xNTotal;
    	integer , dimension(:) , allocatable :: xNStep;

    	nDim    = size(corrL)
    	allocate(xNStep(nDim))
		xNStep  = nPoinsPerCorrL*ceiling((xMax-xMin)/corrL);
		xNTotal = product(xNStep)
		allocate(xPoints(xNTotal, nDim))

		!call DispCarvalhol(xNStep,"xNStep")

		do i = 1, xNTotal
			call get_Permutation(i, xMax, xNStep, xPoints(i,:), xMin);
		end do

		deallocate(xNStep)

    end subroutine set_XPoints

end module randomFieldND
