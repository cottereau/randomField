module statistics_RF

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_Statistics1D(randField, average, stdDeviation, averageCorrL)

		implicit none
		!INPUT
		double precision, dimension(:, :), allocatable, intent(in) :: randField;
		!OUTPUT
		double precision, dimension(:), allocatable, intent(out) :: average, stdDeviation, averageCorrL;
		!LOCAL VARIABLES
		integer :: i;

		!Allocation
		allocate(average(size(randField, 2)));
		allocate(stdDeviation(size(randField, 2)));
		allocate(averageCorrL(size(randField, 2)));

		do i = 1, size(randField, 2)
			average(i) = calculateAverage(randField(:,i));
			stdDeviation(i) = calculateStdDeviation(randField(:,i), average(i));
			averageCorrL(i) = calculateAverageCorrL(randField(:,i))
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

		write (*,*) "calculateAverage(vector**2d0) = ", calculateAverage(vector**2d0);
		write (*,*) "avg**2d0 = ", avg**2d0;

		stdDeviation = (calculateAverage(vector**2d0) - avg**2d0)**(0.5d0);

	end function

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	function calculateAverageCorrL(vector) result(averageCorrL)
		implicit none
		!INPUT
		double precision, dimension(:), intent(in) :: vector;
		!OUTPUT
		double precision :: averageCorrL;
		!LOCAL VARIABLES
		integer :: i;

		averageCorrL = 1234; !Not implemented

	end function
end module statistics_RF
