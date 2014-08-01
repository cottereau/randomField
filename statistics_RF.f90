module statistics_RF

	use displayCarvalhol
	use math_RF
	use mpi

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_StatisticsND(totalRandField, xNStep, xMax, &
					 			totalAverage, totalStdDeviation, totalAverageCorrL)

		implicit none
		!INPUT
		double precision, dimension(:, :),           intent(in) :: totalRandField;
		double precision, dimension(:),              intent(in) :: xMax;
		integer,          dimension(:),              intent(in) :: xNStep;
		!OUTPUT
		double precision, dimension(:), allocatable, intent(out) :: totalAverage, totalStdDeviation, &
																	totalAverageCorrL;
		!LOCAL VARIABLES
		integer :: i, Nmc, nPermut,  nDim;

		write(*,*) "";
		write(*,*) "------------START set_StatisticsND-----------------------------------------";
		write(*,*) "";

		nPermut = size(totalRandField, 1);
		Nmc     = size(totalRandField, 2);
		nDim    = size(xNStep   , 1)

		allocate(totalAverage(nPermut));
		allocate(totalStdDeviation(nPermut));
		allocate(totalAverageCorrL(nDim   ));
		totalAverage = -1
		totalStdDeviation = -1
		totalAverageCorrL = -1
		write(*,*) ">>>>>>>>> Calculating Average and Stantard Deviation";
		do i = 1, nPermut
			totalAverage(i)      = calculateAverage(totalRandField(i,:));
			totalStdDeviation(i) = calculateStdDeviation(totalRandField(i,:), totalAverage(i));
		end do

		if(Nmc > 1) then
			write(*,*) ">>>>>>>>> Calculating Correlation Length";
			call calculateAverageCorrL(totalRandField, xMax, xNStep, totalAverage, totalStdDeviation, totalAverageCorrL)
		else
			write(*,*) ">>>>>>>>> WARNING! Correlation Length coudn't be computed (only one event)";
		end if

		write(*,*) "";
		write(*,*) "------------END set_StatisticsND-----------------------------------------";
		write(*,*) "";

	end subroutine set_StatisticsND

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_StatisticsND_MPI(randField, xNStep, xMax, rang,                  &
									sumRF, sumRFsquare, ptAvg, ptStdDev, procCorrL, &
									totalSumRF, totalSumRFsquare,                   &
									evntAvg, evntStdDev,                            &
									globalAvg, globalStdDev, globalCorrL)

		implicit none
		!INPUT
		integer                          , intent(in) :: rang;
		integer         , dimension(:)   , intent(in) :: xNStep;
		double precision, dimension(:, :), intent(in) :: randField;
		double precision, dimension(:)   , intent(in) :: xMax;
		!OUTPUT
		double precision, dimension(:), allocatable, intent(out) :: sumRF, sumRFsquare, totalSumRF, totalSumRFsquare;
		double precision, dimension(:), allocatable, intent(out) :: evntAvg, evntStdDev, procCorrL;
    	double precision, dimension(:), allocatable, intent(out) :: ptAvg, ptStdDev, globalCorrL;
    	double precision                           , intent(out) :: globalAvg, globalStdDev;

		!LOCAL VARIABLES
		integer :: i, Nmc, nPermut,  nDim;
		integer :: sizeLoc;
		integer :: code, xNStepTotal;

		write(*,*) "";
		write(*,*) "------------START set_StatisticsNDmpi (proc = ", rang, ")-----------------------------------------";
		write(*,*) "";

		write (*,*) ">>>>get_sizes_MPI Proc = ", rang
		call get_sizes_MPI (xNStep, sizeLoc)
		nDim        = size(xNStep)
		xNStepTotal = product(xNStep)
		Nmc         = size(randField,2)

		write (*,*) ">>>>Allocating Proc = ", rang

		!Allocating
		allocate(sumRF(Nmc))
		allocate(sumRFsquare(Nmc))
		allocate(ptAvg(sizeLoc))
		allocate(ptStdDev(sizeLoc))
		allocate(procCorrL(nDim));
		if(rang == 0) then
			allocate (totalSumRF(Nmc))
			allocate (totalSumRFsquare(Nmc))
			totalSumRF = -1.0d0
			totalSumRFsquare = -1.0d0
			allocate (evntAvg(Nmc))
			allocate (evntStdDev(Nmc))
		end if


		!write (*,*) ">>>>sizeLoc = ", sizeLoc, " Proc = ", rang

		sumRF(:)       = sum( randField(1:sizeLoc, :)    , dim = 1)
		sumRFsquare(:) = sum((randField(1:sizeLoc, :))**2, dim = 1)

		!call DispCarvalhol(sumRF, "sumRF");
		!call DispCarvalhol(sumRFsquare, "sumRFsquare", "F15.8");

		!Calculating events statistics
		!write (*,*) ">>>>MPI_REDUCE Proc = ", rang

		call MPI_REDUCE(sumRF, totalSumRF, Nmc, MPI_DOUBLE_PRECISION, &
					    MPI_SUM, 0, MPI_COMM_WORLD, code)
		call MPI_REDUCE(sumRFsquare, totalSumRFsquare, Nmc, MPI_DOUBLE_PRECISION, &
					MPI_SUM, 0, MPI_COMM_WORLD, code)

		!write (*,*) ">>>>after MPI_REDUCE Proc = ", rang

		if(rang == 0) then
			evntAvg      = totalSumRF/dble(xNStepTotal);
			evntStdDev   = sqrt(totalSumRFsquare/dble(xNStepTotal) &
		             	   - (totalSumRF/dble(xNStepTotal))**2)
			globalAvg    = sum(totalSumRF)/dble(xNStepTotal*Nmc);
			globalStdDev = sqrt(sum(totalSumRFsquare)/dble(xNStepTotal*Nmc) &
		              	   - (sum(totalSumRF)/dble(xNStepTotal*Nmc))**2)
!
		    call DispCarvalhol(evntAvg, "evntAvg")
		    call DispCarvalhol(evntStdDev, "evntStdDev")
		    write(*,*) "globalAvg    = ", globalAvg
		    write(*,*) "globalStdDev = ", globalStdDev

!		write(*,*) ">>>>>>>>> Showing Statistics of each event";
!		write(*,*) "totalSumRF = ", totalSumRF;
!		write(*,*) "totalSumRFsquare = ", totalSumRFsquare;
!		write(*,*) " Average = ", totalSumRF/dble(xNStepTotal);
!		write(*,*) "Variance = ", totalSumRFsquare/dble(xNStepTotal) - (totalSumRF/dble(xNStepTotal))**2;
!		write(*,*) " Std Dev = ", sqrt(totalSumRFsquare/dble(xNStepTotal) - (totalSumRF/dble(xNStepTotal))**2)
!		write(*,*) ">>>>>>>>> Showing Global Statistics (all events)";
!		write(*,*) " Average = ", sum(totalSumRF)/dble(xNStepTotal*Nmc);
!		write(*,*) "Variance = ", sum(totalSumRFsquare)/dble(xNStepTotal*Nmc) - (sum(totalSumRF)/dble(xNStepTotal*Nmc))**2;
!		write(*,*) " Std Dev = ", sqrt(sum(totalSumRFsquare)/dble(xNStepTotal*Nmc) - (sum(totalSumRF)/dble(xNStepTotal*Nmc))**2)
		end if

		do i = 1, sizeLoc
			ptAvg(i)    = calculateAverage(randField(i,:));
			ptStdDev(i) = calculateStdDeviation(randField(i,:), ptAvg(i));
		end do

		!write (*,*) "Rang = ", rang
		!call DispCarvalhol(ptAvg, "ptAvg")
		!if(rang == 0) call DispCarvalhol(ptStdDev, "ptStdDev")

		call calculateAverageCorrL_MPI(randField, xMax, xNStep, ptAvg, ptStdDev, procCorrL)

		!write(*,*) "      sumRF = ", sumRF
		!write(*,*) "sumRFsquare = ", sumRFsquare

		!write(*,*) "";
		!write(*,*) "------------END set_StatisticsNDmpi (proc = ", rang, ")-----------------------------------------";
		!write(*,*) "";

	end subroutine set_StatisticsND_MPI

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
	subroutine calculateAverageCorrL(randField, xMax, xNStep, ptAvg, ptStdDev, averageCorrL)
		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in) :: randField;
		double precision, dimension(:),    intent(in) :: xMax;
		double precision, dimension(:),    intent(in) :: ptAvg, ptStdDev;
		integer,          dimension(:),    intent(in) :: xNStep;
		!OUTPUT
		double precision, dimension(:),    intent(out) :: averageCorrL;
		!LOCAL VARIABLES
		integer :: i, j, k, Nmc, nPermut, nPlanes, nDim, &
				   beg, step, end;
		!double precision, dimension(:), allocatable :: Ra, Rb, Rc, R
		double precision :: corrL

		nPermut      = product(xNStep);
		Nmc          = size(randField, 2);
		nDim         = size(xNStep);
		averageCorrL = 0;


		do i = 1, nDim
			nPlanes = nPermut/xNStep(i)
			do j = 1, nPlanes
				call get_SequenceParam(i, j, xNStep, beg, step, end)

				if(end > nPermut) end = nPermut

				averageCorrL(i) = sum(((1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),randField(beg,:))) &
				                      - (ptAvg(beg+step:end:step)    * ptAvg(beg)))                               &
									  / (ptStdDev(beg+step:end:step) * ptStdDev(beg)))                            &
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
	subroutine calculateAverageCorrL_MPI(randField, xMax, xNStep, ptAvg, ptStdDev, procCorrL)
		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in) :: randField;
		double precision, dimension(:),    intent(in) :: xMax;
		double precision, dimension(:),    intent(in) :: ptAvg, ptStdDev;
		integer,          dimension(:),    intent(in) :: xNStep;
		!OUTPUT
		double precision, dimension(:),    intent(out) :: procCorrL;
		!LOCAL VARIABLES
		integer :: i, j, k, Nmc, nPermut, nPlanes, nDim, rang, code, &
				   beg, step, end, supLim, infLim, spam;
		!double precision, dimension(:), allocatable :: Ra, Rb, Rc, R
		double precision :: corrL, countPlanes

		call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)

		nPermut     = product(xNStep);
		Nmc         = size(randField, 2);
		nDim        = size(xNStep);
		procCorrL   = 0.0d0;
		countPlanes = 0.0d0

		call get_sizes_MPI(xNStep, start = infLim, end = supLim)

		do i = 1, nDim
			nPlanes = nPermut/xNStep(i)
			countPlanes = 0.0d0
			do j = 1, nPlanes
				!Reasoning in global coordinates (the whole random field)
				call get_SequenceParam(i, j, xNStep, beg, step, end)
				spam = end - beg + 1
				beg = step*(infLim/step) + beg !lower cut limitation
				if(step*(infLim/step) /= infLim) beg = beg + step !Round up
				if(end > supLim) end = supLim !upper cut limitation
				if(beg + step < supLim) then
					!Transformate the coordinates in local again
					beg = beg - infLim + 1
					end = end - infLim + 1
					write(*,*) "Proc = ", rang, "Plane = ", j, "beg = ", beg, "end = ", end, "infLim = ", infLim, "supLim = ", supLim
					countPlanes  = countPlanes + 1.0d0
!					procCorrL(i) = sum(((1d0/dble(Nmc) * matmul(randField(beg+step:end:step, :),randField(beg,:))) &
!										- (ptAvg(beg+step:end:step)    * ptAvg(beg)))                                 &
!										/ (ptStdDev(beg+step:end:step) * ptStdDev(beg)))                           &
!									    + procCorrL(i)
				end if
			end do
			if (countPlanes > 0.0d0) procCorrL(i) = (xMax(i)/xNStep(i))*procCorrL(i) / countPlanes
		end do
		procCorrL = 2*procCorrL !Symmetry

	end subroutine calculateAverageCorrL_MPI

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
