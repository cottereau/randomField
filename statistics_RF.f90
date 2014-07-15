module statistics_RF

	use displayCarvalhol
	use math_RF

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_StatisticsND(randField, xNStep, xMax, average, stdDeviation, averageCorrL)

		implicit none
		!INPUT
		double precision, dimension(:, :),           intent(in) :: randField;
		double precision, dimension(:),              intent(in) :: xMax;
		integer,          dimension(:),              intent(in) :: xNStep;
		!OUTPUT
		double precision, dimension(:), allocatable, intent(out) :: average, stdDeviation, &
																	averageCorrL;
		!LOCAL VARIABLES
		integer :: i, Nmc, nPermut,  nDim;

		write(*,*) "";
		write(*,*) "------------START set_StatisticsND-----------------------------------------";
		write(*,*) "";

		nPermut = size(randField, 1);
		Nmc     = size(randField, 2);
		nDim    = size(xNStep   , 1)

		allocate(average     (nPermut));
		allocate(stdDeviation(nPermut));
		allocate(averageCorrL(nDim   ));
		average = -1
		stdDeviation = -1
		averageCorrL = -1
		write(*,*) ">>>>>>>>> Calculating Average and Stantard Deviation";
		do i = 1, nPermut
			average(i)      = calculateAverage(randField(i,:));
			stdDeviation(i) = calculateStdDeviation(randField(i,:), average(i));
		end do

		if(Nmc > 1) then
			write(*,*) ">>>>>>>>> Calculating Correlation Length";
			call calculateAverageCorrL(randField, xMax, xNStep, average, stdDeviation, averageCorrL)
		else
			write(*,*) ">>>>>>>>> WARNING! Correlation Length coudn't be computed (only one event)";
		end if

		write(*,*) "";
		write(*,*) "------------END set_StatisticsND-----------------------------------------";
		write(*,*) "";

	end subroutine set_StatisticsND

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

	end function calculateAverage

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

	end function calculateStdDeviation

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine calculateAverageCorrL(randField, xMax, xNStep, average, stdDeviation, averageCorrL)
		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in) :: randField;
		double precision, dimension(:),    intent(in) :: xMax;
		double precision, dimension(:),    intent(in) :: average, stdDeviation;
		integer,          dimension(:),    intent(in) :: xNStep;
		!OUTPUT
		double precision, dimension(:),    intent(out) :: averageCorrL;
		!LOCAL VARIABLES
		integer :: i, j, k, Nmc, nPermut, nPlanes, nDim, &
				   beg, step, end;
		!double precision, dimension(:), allocatable :: Ra, Rb, Rc, R
		double precision :: corrL

		nPermut      = size(randField, 1);
		Nmc          = size(randField, 2);
		nDim         = size(xNStep   , 1);
		averageCorrL = 0;


		do i = 1, nDim
			nPlanes = nPermut/xNStep(i)
			do j = 1, nPlanes
				call get_SequenceParam(i, j, xNStep, beg, step, end)

				averageCorrL(i) = sum(((1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),randField(beg,:))) - &
									  (average(beg+step:end:step)*average(beg))) / &
									  (stdDeviation(beg+step:end:step)*stdDeviation(beg))) &
								       + averageCorrL(i)



!					Ra = 1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),&
!											    randField(beg,:))
!					Rb = average(beg+step:end:step)*average(beg)
!					Rc = stdDeviation(beg+step:end:step)*stdDeviation(beg)
!					R = (Ra - Rb) / Rc
!					averageCorrL(i) = sum(R) * xMax(i,1)/xNStep(i,1) + averageCorrL(i)
!
!					call dispCarvalhol(1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),&
!											    randField(beg,:)),"Ra")
!					call dispCarvalhol(average(beg+step:end:step)*average(beg),"Rb")
!					call dispCarvalhol(stdDeviation(beg+step:end:step)*stdDeviation(beg),"Rc")
!					call dispCarvalhol(((1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),randField(beg,:))) - &
!										  (average(beg+step:end:step)*average(beg))) / &
!										  (stdDeviation(beg+step:end:step)*stdDeviation(beg)),"R")
!					write(*,*) "R"
!					write(*,*)(((1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),randField(beg,:))) - &
!										  (average(beg+step:end:step)*average(beg))) / &
!										  (stdDeviation(beg+step:end:step)*stdDeviation(beg)))
!					write(*,*) "averageCorrL(",i,") = ", averageCorrL(i)

			end do
			averageCorrL(i) = (xMax(i)/xNStep(i))*averageCorrL(i) / dble(nPlanes)
		end do
		averageCorrL = 2*averageCorrL !Symmetry

	end subroutine calculateAverageCorrL

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

    end subroutine get_SequenceParam

end module statistics_RF
