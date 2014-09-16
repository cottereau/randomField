module spectra_RF
	use displayCarvalhol
	use math_RF
contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine set_kMaxND(corrMod, corrL, kMax);
    	!INPUT
        character (len=15),                intent(in) :: corrMod;
        double precision,   dimension(:),  intent(in) :: corrL;

        !OUTPUT
        double precision, dimension(:),   intent(out) :: kMax;

        !LOCAL VARIABLES
        double precision :: pi = 3.1415926535898
        integer          :: i

!        do i = 1, 100
!        	kMax = i/10.0 * corrL(:)
!        	write(*,*) "kMax = ", kMax
!        	write(*,*) "Spectrum = ", get_SpectrumND(kMax, corrMod, corrL)
!        	call DispCarvalhol (kMax, "kMax")
!        	call DispCarvalhol (get_SpectrumND(kMax, corrMod, corrL), "Spectrum")
!        end do

    	select case(corrMod)
   				case("gaussian")
				kMax(:) = 2*pi*corrL(:); !CRITERIA STILL TO BE TESTED
   		end select

    end subroutine set_kMaxND

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine set_kSign(kSign);
    	!INPUT and OUTPUT
        double precision,   dimension(:, :), intent(inout) :: kSign;

        !LOCAL VARIABLES
        integer:: i, j, pos, seedStep;

		nDim = size(kSign, 2);

		kSign(:,:) = 1d0

		do i = 2, nDim
			do j = 1, size(kSign,1)
			    seedStep = 2**(nDim-i);
				pos      = cyclicMod(int((j-0.9)/seedStep), 2)
				if (mod(pos,2) == 1) kSign(j,i) = -1d0
			end do
		end do
    end subroutine set_kSign

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    function get_SpectrumND(kVector, corrMod, corrL) result (Sk)
        ! Return Spectrum from a chosen correlation model
        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: kVector;
        double precision, dimension(:), intent(in) :: corrL;
        character (len=*),              intent(in) :: corrMod

        !OUTPUT
        double precision :: Sk;

        !LOCAL VARIABLES
        integer :: j, nDim;
        double precision, dimension(:), allocatable :: eta;
        double precision :: pi = 3.1415926535898

        Sk = 0;
        nDim = size(kVector)
        allocate (eta(nDim))

        select case(corrMod)
            case("gaussian" : "lognormal")
				Sk  = exp(-dot_product((kVector**2),(corrL**2))/(4.0d0)); !Amplitude part "product(corrL)" is external to the function
!            	eta = corrL/(2*sqrt(pi))
!				Sk  = exp(-dot_product((eta**2),(kVector**2))); !Amplitude part "product(corrL)" is external to the function
        end select

        deallocate (eta)
    end function get_SpectrumND

end module spectra_RF
