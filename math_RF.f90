module math_RF

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

end module math_RF
