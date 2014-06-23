module statistics_RF

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_Statistics1D(randField, xNStep, kNStep, average, stdDeviation, averageCorrL)

		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in) :: randField;
		integer,          dimension(:, :), intent(in) :: xNStep, kNStep;
		!OUTPUT
		double precision, dimension(:, :), intent(out) :: average, stdDeviation, averageCorrL;
		!LOCAL VARIABLES
		integer :: i;

		do i = 1, size(randField, 2)
			average(1,i)      = calculateAverage&
								(randField(1:xNStep(i, 1),i));
			stdDeviation(1,i) = calculateStdDeviation&
								(randField(1:xNStep(i, 1),i), average(1,i));
			averageCorrL(1,i) = calculateAverageCorrL&
								(randField(1:xNStep(i, 1),i), average(1,i), stdDeviation(1,i));
		end do
	end subroutine

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_StatisticsND(randField, xNStep, kNStep, average, stdDeviation, averageCorrL)

		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in) :: randField;
		integer,          dimension(:, :), intent(in) :: xNStep, kNStep;
		!OUTPUT
		double precision, dimension(:, :), intent(out) :: average, stdDeviation, averageCorrL;
		!LOCAL VARIABLES
		integer :: i, Nmc, totalNStep;

		Nmc = size(randField, 2);

		do i = 1, Nmc
			totalNStep        = product(xNStep(:,i))
			average(1,i)      = calculateAverage&
								(randField(1:totalNStep,i));
			stdDeviation(1,i) = calculateStdDeviation&
								(randField(1:totalNStep,i), average(1,i));
			averageCorrL(1,i) = calculateAverageCorrL&
								(randField(1:totalNStep,i), average(1,i), stdDeviation(1,i));
		end do
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
	function calculateAverageCorrL(vector, average, stdDeviation) result(averageCorrL)
		implicit none
		!INPUT
		double precision, dimension(:), intent(in) :: vector;
		double precision, intent(in), optional :: average, stdDeviation;
		!OUTPUT
		double precision :: averageCorrL;
		!LOCAL VARIABLES
		integer :: i, j, t;
		double precision :: avg, stdDev;

		if(present(average)) then
			avg = average;
		else
			avg = calculateAverage(vector);
		endif

		if(present(average)) then
			stdDev = stdDeviation;
		else
			stdDev = calculateStdDeviation(vector);
		endif

		averageCorrL = 0;
		t=0
		do i = 1, size(vector)
			do j = i+1, size(vector)
				averageCorrL = averageCorrL + ((vector(i)-avg)*(vector(j)-avg))
				t = t + 1
			end do
		end do

		averageCorrL = averageCorrL/(t*stdDev**2);

	end function
end module statistics_RF
