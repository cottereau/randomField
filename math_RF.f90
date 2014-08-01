module math_RF

	use mpi
    !All logic and math routines

contains

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

    end function cyclicMod

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine get_Permutation(pos, qmax, nStep, pVec, qmin)

        implicit none

        !INPUT
        integer                       , intent(in)           :: pos;
        double precision, dimension(:), intent(in)           :: qmax;
        double precision, dimension(:), intent(in), optional :: qmin;
        integer,          dimension(:), intent(in)           :: nStep;
        !OUTPUT
        double precision, dimension(:), intent(out) :: pVec;
        !LOCAL VARIABLES
        integer :: i, j;
        integer :: seedStep, nDim;

		nDim = size(nStep);

		if (present(qmin)) then
			do j = 1, nDim
				seedStep = product(nStep(j+1:));
		 		if (j == nDim) seedStep = 1;
				i = cyclicMod(int((pos-0.9)/seedStep)+1, nStep(j))
				pVec(j) = (dble(i)-0.5d0)*(qmax(j)-qmin(j))/(nStep(j));
			end do
		else
			do j = 1, nDim
		  		seedStep = product(nStep(j+1:));
		 		if (j == nDim) seedStep = 1;
				i = cyclicMod(int((pos-0.9)/seedStep)+1, nStep(j))
				pVec(j) = (dble(i)-0.5d0)*(qmax(j))/(nStep(j)); !qmin = 0
			end do
		end if

    end subroutine get_Permutation

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine get_sizes_MPI(xNStep, sizeLoc, sizeUnif, start, end)

    	!INPUT
    	integer, dimension(:), intent(in) :: xNStep;
    	!OUTPUT
        integer, intent(out), optional :: sizeLoc, sizeUnif;
        integer, intent(out), optional :: start, end;

        !LOCAL VARIABLES
        integer :: rang, nb_procs, code;
        integer :: xStart, xEnd, xNStepTotal;

    	call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

		xNStepTotal = product(xNStep)

    	if(rang == nb_procs-1) then
			xStart = (nb_procs-1)*ceiling(dble(xNStepTotal)/dble(nb_procs)) + 1
			xEnd   = xNStepTotal
		else
			xStart = rang*ceiling(dble(xNStepTotal)/dble(nb_procs)) + 1
			xEnd   = (rang+1)*ceiling(dble(xNStepTotal)/dble(nb_procs))
		end if

		if(present(sizeLoc))  sizeLoc  = xEnd - xStart + 1 !Used to escape the not-used places of sizeUnif
		if(present(sizeUnif)) sizeUnif = ceiling(dble(xNStepTotal)/dble(nb_procs)) !Used to allow the MPI_GATHER afterwards
		if(present(start))    start    = xStart
		if(present(end))      end      = xEnd

	end subroutine get_sizes_MPI

end module math_RF
