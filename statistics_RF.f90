module statistics_RF

	use displayCarvalhol
	use math_RF

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_StatisticsND(randField, xNStep, xMax, average, stdDeviation, averageCorrL)

		implicit none
		!INPUT
		double precision, dimension(:, :),              intent(in) :: randField, xMax;
		integer,          dimension(:, :),              intent(in) :: xNStep;
		!OUTPUT
		double precision, dimension(:), allocatable, intent(out) :: average, stdDeviation, &
																	averageCorrL;
		!LOCAL VARIABLES
		integer :: i, j, k, Nmc, nPermut, nPlanes, nDim, &
				   beg, step, end;
		double precision, dimension(:), allocatable :: Ra, Rb, Rc, R
		double precision :: corrL

		write(*,*) "";
		write(*,*) "------------START set_StatisticsND-----------------------------------------";
		write(*,*) "";

		nPermut = size(randField, 1);
		Nmc     = size(randField, 2);
		nDim    = size(xNStep   , 1)

		allocate(average     (nPermut));
		allocate(stdDeviation(nPermut));
		allocate(averageCorrL(nDim   ));

		write(*,*) ">>>>>>>>> Calculating Average and Stantard Deviation";
		do i = 1, nPermut
			average(i)      = calculateAverage(randField(i,:));
			stdDeviation(i) = calculateStdDeviation(randField(i,:), average(i));
		end do

		write(*,*) ">>>>>>>>> Calculating Correlation Length";
		call calculateAverageCorrL(randField, xMax, xNStep, average, stdDeviation, averageCorrL)

		write(*,*) "";
		write(*,*) "------------END set_StatisticsND-----------------------------------------";
		write(*,*) "";

	end subroutine

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	function calculateAverage(vector) result(average)

		implicit none
		!INPUT
		double precision, dimension(:), intent(in) :: vector;
		!OUTPUT
		double precision :: average;
		!LOCAL VARIABLES
		integer :: i;


		average = 0d0;

		do i = 1, size(vector)
			average = average + vector(i);
		end do

		average = average/dble(size(vector))

	end function

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	function calculateStdDeviation(vector, average) result(stdDeviation)
		implicit none
		!INPUT
		double precision, dimension(:), intent(in) :: vector;
		double precision, intent(in), optional :: average;
		!OUTPUT
		double precision :: stdDeviation;
		!LOCAL VARIABLES
		integer :: i;
		double precision :: avg;


		if(present(average)) then
			avg = average;
		else
			avg = calculateAverage(vector);
		endif

		stdDeviation = (calculateAverage(vector**2d0) - avg**2d0)**(0.5d0);

	end function

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine calculateAverageCorrL(randField, xMax, xNStep, average, stdDeviation, averageCorrL)
		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in) :: randField, xMax;
		double precision, dimension(:),    intent(in) :: average, stdDeviation;
		integer,          dimension(:, :), intent(in) :: xNStep;
		!OUTPUT
		double precision, dimension(:) :: averageCorrL;
		!LOCAL VARIABLES
		integer :: i, j, k, Nmc, nPermut, nPlanes, nDim, &
				   beg, step, end;
		double precision, dimension(:), allocatable :: Ra, Rb, Rc, R
		double precision :: corrL
		logical          :: possible

		nPermut  = size(randField, 1);
		Nmc      = size(randField, 2);
		nDim     = size(xNStep   , 1)
		possible = .true.
		averageCorrL = 0

		do i = 1, Nmc
			do j = 1, nDim
				if (xMax(j,i) /= xMax(j,1) .or. xNStep(j,i) /= xNStep(j,1)) then
					possible = .false.
					write(*,*) "Random field realizations are not compatible to calculate the correlation length"
					exit
				end if
			end do
			if (.not.possible) exit
		end do

		if (possible) then
			do i = 1, nDim
!				write(*,*) "Dimension ", i
				nPlanes = nPermut/xNStep(i,1)
				do j = 1, nPlanes
					call get_SequenceParam(i, j, xNStep(:,1), beg, step, end)

					Ra = 1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),&
											    randField(beg,:))
					Rb = average(beg+step:end:step)*average(beg)
					Rc = stdDeviation(beg+step:end:step)*stdDeviation(beg)
					R = (Ra - Rb) / Rc
					averageCorrL(i) = sum(R) * xMax(i,1)/xNStep(i,1) + averageCorrL(i)

!					call disp1D(Ra,"Ra")
!					call disp1D(Rb,"Rb")
!					call disp1D(Rc,"Rc")
!					call disp1D(R,  "R")
!					write(*,*) R

				end do
				averageCorrL(i) = averageCorrL(i) / dble(nPlanes)
			end do
		end if

	end subroutine

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine get_SequenceParam(axis, nPlane, nStep, beg, step, end)

        implicit none

        !INPUT
        integer,               intent(in)    :: axis, nPlane
        integer, dimension(:), intent(in) :: nStep;

        !OUTPUT
        integer, intent(out) :: beg, step, end;

        !LOCAL VARIABLES
        integer :: i, nDim, position, posPlane, gap;

		beg  = 0
		nDim = size(nStep)

		step = product(nStep(axis+1:));
		if (i == nDim) step = 1;
		gap = step*nStep(axis)
		beg = cyclicMod(nPlane, step) &
			  + int((nPlane-0.9)/step)*gap

		end = beg + step*(nStep(axis)-1)

    end

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    function get_Positions(vectorPos, fixedDim, nStep) result(posVec)
!
!        implicit none
!
!        !INPUT
!        integer,               intent(in) :: vectorPos, fixedDim;
!        integer, dimension(:), intent(in) :: nStep;
!        !OUTPUT
!        integer, dimension(:), allocatable :: posVec;
!        !LOCAL VARIABLES
!        integer :: i, j, pos;
!        integer :: seedStep, gapStep, nDim, total;
!
!		nDim     = size(nStep);
!		pos      = cyclicMod(vectorPos, nStep(fixedDim))
!
!		!Defining Seed Step
!	    seedStep = product(nStep(fixedDim+1:));
!	    if (fixedDim == nDim) seedStep = 1;
!
!	    !Defining Gap Step
!	    gapStep = seedStep * nStep(fixedDim)
!
!	    !Allocating vector
!		total = product(nStep)
!		allocate(posVec(total/nStep(fixedDim)))
!
!		do j = 1, total
!			i =   int((j-0.9)/seedStep)*gapStep &
!			    + (pos-1)*seedStep &
!				+ cyclicMod(j, seedStep)
!			posVec(j) = i;
!		end do
!    end

end module statistics_RF
