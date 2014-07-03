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

    end

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    function get_Permutation(pos, qmin, qmax, nStep) result(pVec)

        implicit none

        !INPUT
        integer                       , intent(in) :: pos;
        double precision, dimension(:), intent(in) :: qmin, qmax;
        integer,          dimension(:), intent(in) :: nStep;
        !OUTPUT
        double precision, dimension(:), allocatable :: pVec;
        !LOCAL VARIABLES
        integer :: i, j;
        integer :: seedStep, nDim;

		nDim = size(nStep);
		allocate(pVec(nDim))

		do j = 1, nDim
		    seedStep = product(nStep(j+1:));
		    if (j == nDim) seedStep = 1;
			i = cyclicMod(int((pos-0.9)/seedStep)+1, nStep(j))
			pVec(j) = (dble(i)-0.5d0)*(qmax(j)-qmin(j))/(nStep(j));
		end do
    end

end module math_RF
