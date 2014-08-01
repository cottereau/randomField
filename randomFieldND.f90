module randomFieldND

    use displayCarvalhol
    use spectra_RF
    use math_RF
	use mpi

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine createRandomFieldND(Nmc, corrMod, corrL, xMax, xPeriod, &
    						       kMax, xNStep, kNStep, randField);
        implicit none

        !INPUT
        integer,                          intent(in) :: Nmc;
        character (len=*),                intent(in) :: corrMod;
        double precision,   dimension(:), intent(in) :: corrL, xMax, xPeriod;

        !OUTPUT
        double precision, dimension(:),    allocatable, intent(out) :: kMax;
        integer,          dimension(:),    allocatable, intent(out) :: xNStep, kNStep;
        double precision, dimension(:, :), allocatable, intent(out) :: randField;

        !LOCAL VARIABLES
        integer :: i, j, k, m, nDim, kNS;
        double precision, dimension(:,:), allocatable :: kSign, phiN;
        double precision, dimension(:),   allocatable :: xVec, kVec, kVecUnsigned; !Allocated in function
        double precision, dimension(:),   allocatable :: deltaK, randLoc;
        double precision :: Sk;
        double precision :: pi = 3.1415926535898, zero = 0d0;
        integer          :: xNStepTotal, kNStepTotal;
        integer          :: testDim; !Only for tests
        integer          :: nb_procs, rang, code, xStart, xEnd, sizeUnif, sizeLoc;

		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
		call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)

		if(rang == 0) write(*,*) "";
        if(rang == 0) write(*,*) "------------START randomFieldND (proc = ", rang, ")---------";
		if(rang == 0) write(*,*) "";

		nDim = size(xMax,1);

		if(rang == 0) write(*,*) ">>>>>>>>> Variables initialization: kMax, xNStep and kNStep";
		!Allocating
		allocate(kMax       (nDim));
		allocate(xNStep     (nDim));
		allocate(kNStep     (nDim));
		!Initializing
		call set_kMaxND(corrMod, corrL, kMax)
		xNStep = ceiling(xMax/(pi/kMax))+1;
		kNStep = ceiling(kMax/(2.0d0*pi/xPeriod))+1;


		!ONLY FOR TESTS (Overwriting)----------------------

		testDim = 5;

		kMax  (1) = 2*pi*corrL(1);
	    xNStep(1) = testDim;
	    kNStep(1) = testDim;

		if(nDim > 1) then
			kMax  (2) = 2*pi*corrL(2);
		    xNStep(2) = testDim + 1;
			kNStep(2) = testDim + 1;
		end if

		if(nDim > 2) then
			kMax  (3) = 2*pi*corrL(3);
		    xNStep(3) = testDim + 2;
			kNStep(3) = testDim + 2;
		end if
		!---------------------------------------------------

		xNStepTotal = product(xNStep);
		kNStepTotal = product(kNStep);

		call get_sizes_MPI(xNStep, sizeLoc, sizeUnif, xStart, xEnd)

	    write(*,*) "Proc ", rang, " xStart = ", xStart, "   xEnd = ", xEnd, &
	    	       "sizeUnif"  , sizeUnif, "sizeLoc"   , sizeLoc

		!Random Field
		allocate(randField((sizeUnif),Nmc));
		allocate(randLoc   (sizeUnif));
		allocate(deltaK   (nDim));
		allocate(kSign    (2**(nDim-1), nDim));
		allocate(phiN     (size(kSign,1), kNStepTotal));
		allocate(kVecUnsigned(nDim));
		allocate(kVec(nDim));
		allocate(xVec(nDim));

		if(rang == 0) then
			write(*,*) "";
			write(*,*) ">>>>>>>>> Random Field Creation";
			write(*,*) "";
			write(*,*) "Number of dimensions = ", nDim;
			write(*,*) "Number of events     = ", Nmc;
			write(*,*) "Number x points      = ", xNStepTotal;
			write(*,*) "Number k points      = ", kNStepTotal;
			write(*,*) "";
		end if

		randField = 0;
		randLoc   = 0;
		deltaK    = 0;
	    deltaK    = (kMax)/(kNStep-1); !Defines deltaK
	    kNS       = kNStepTotal;
		call set_kSign(kSign)

		do k = 1, Nmc
			call random_number(phiN(:,1:kNS))
			do j = 1, kNS
				call get_Permutation(j, kMax, kNStep, kVecUnsigned);
			    do m = 1, size(kSign,1)
			    	kVec           = kVecUnsigned * kSign(m, :)
					Sk             = get_SpectrumND(kVec, corrMod, corrL);
					randLoc        = 0;
					randLoc(1:sizeLoc) = 2*pi*phiN(m, j) !TO VERIFY - Random part of the angle matrix for this k permutation
					do i = 1, sizeLoc
						call get_Permutation(i+(xStart-1), xMax, xNStep, xVec);
						randLoc(i) = randLoc(i) + dot_product(kVec, xVec); !Not-random part of the angle matrix
					end do
					randLoc(1:sizeLoc)     = sqrt(Sk*2*product(deltaK)) &
										     *cos(randLoc(1:sizeLoc))
					randField(1:sizeLoc,k) = randField(1:sizeLoc,k) + randLoc(1:sizeLoc)
				end do
			end do
			randField(1:sizeLoc,k) = product(corrL)*randField(1:sizeLoc,k)
			if(rang == 0) write(*,*) "Event ",k, "of", Nmc, "completed (counting only in proc 0)";
		end do

		randField(1:sizeLoc,:) = sqrt(2d0)*randField(1:sizeLoc,:)

		!Only printing------------------------------------------

!		call DispCarvalhol(xNStep, "xNStep");
!		call DispCarvalhol(kNStep, "kNStep");
!
!		write(*,*) "Permutation X (Nmc = 1)"
!		do i = 1, xNStepTotal
!			call get_Permutation(i, xMax, xNStep, xVec);
!			write(*,*) i, ">", xVec;
!		end do
!
!		write(*,*) "Permutation K (Nmc = 1)"
!		do i = 1, kNStepTotal
!		call get_Permutation(j, kMax, kNStep, kVecUnsigned);
!			write(*,*) i, ">", kVecUnsigned;
!		end do
!
!		!Spectrum
!		write(*,*) "Spectrum (Nmc = 1)"
!		do i = 1, kNStepTotal
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

		if(rang == 0) write(*,*) "";
        if(rang == 0) write(*,*) "------------END randomFieldND (proc = ", rang, "---------";
		if(rang == 0) write(*,*) "";
    end subroutine createRandomFieldND

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine gather_RF_MPI(randField, totalRandField, xNStepTotal, rang)
    	implicit none
    	!INPUT
    	double precision, dimension(:, :), intent(in) :: randField
    	integer                          , intent(in) :: rang, xNStepTotal

    	!OUTPUT
    	double precision, dimension(:, :), allocatable, intent(out) :: totalRandField;

    	!LOCAL VARIABLES
    	integer :: nb_procs, code, Nmc;
    	integer :: type_RF, type_Temp, dblBitSize, newExtent;
    	double precision, dimension(:, :), allocatable :: tempRandField;

		Nmc = size(randField, 2);
    	call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
    	if(rang == 0) then
			allocate (tempRandField(nb_procs*size(randField,1), size(randField,2)))
			tempRandField = 0.0d0
		end if

		!Assembling all the fields in one matrix";
		call MPI_TYPE_VECTOR(size(randField,2), size(randField,1), size(randField,1)*nb_procs, &
		                     MPI_DOUBLE_PRECISION, type_Temp, code) !type_Temp and code are outputs
		call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dblBitSize, code)
		newExtent = size(randField,1)*dblBitSize
		call MPI_TYPE_CREATE_RESIZED(type_Temp, 0, newExtent, type_RF, code) !Changing the starting point to the next "slice"
		call MPI_TYPE_COMMIT(type_RF, code)
		call MPI_GATHER(randField,     size(randField), MPI_DOUBLE_PRECISION, &
		                tempRandField, 1,               type_RF,               &
		                0,             MPI_COMM_WORLD,  code)
		call MPI_TYPE_FREE(type_RF, code)

		!Puting the unused spaces in the trash
		if(rang == 0) then
			allocate (totalRandField(xNStepTotal, Nmc))
			totalRandField = 0.0d0
			totalRandField = tempRandField(1:xNStepTotal,:)
			deallocate(tempRandField)
		end if

    end subroutine gather_RF_MPI

end module randomFieldND
