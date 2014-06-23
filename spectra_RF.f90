module spectra_RF
	use displayCarvalhol
contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    function make_Spectrum1D(kVector, corrMod, corrL) result (Sk)
        ! Return Spectrum from a chosen correlation model
        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: kVector;
        character (len=*), intent(in) :: corrMod
        double precision, intent(in) :: corrL;

        !OUTPUT
        double precision, dimension(:), allocatable :: Sk;

        !LOCAL VARIABLES
        integer :: j;
        double precision :: eta;
        double precision :: pi = 3.1415926535898

        allocate (Sk(size(kVector)));
        Sk = 0;

        select case(corrMod)

            case("gaussian")
				eta = corrL/sqrt(pi)
				Sk = ((1./(2.*pi)) * exp((-(eta**2)/4)*(kVector**2)));

        end select
    end

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine set_kMax1D(corrMod, corrL, kMax);
    	!INPUT
        character (len=15), dimension(:),    intent(in) :: corrMod;
        double precision,   dimension(:, :), intent(in) :: corrL;

        !OUTPUT
        double precision, dimension(:, :), intent(out) :: kMax;

        !LOCAL VARIABLES
        integer:: i;
        double precision :: pi = 3.1415926535898

        do i = 1, size(corrMod)
        	select case(corrMod(i))
   				case("gaussian")
					kMax(i,:) = 2*pi*corrL(i,:); !NOT REALLY IMPLEMENTED, CRITERIA TO BE DEFINED
   			end select
   		end do
    end

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine set_kMaxND(corrMod, corrL, kMax);
    	!INPUT
        character (len=15), dimension(:),    intent(in) :: corrMod;
        double precision,   dimension(:, :), intent(in) :: corrL;

        !OUTPUT
        double precision, dimension(:, :), intent(out) :: kMax;

        !LOCAL VARIABLES
        integer:: i;
        double precision :: pi = 3.1415926535898

        do i = 1, size(corrMod)
        	select case(corrMod(i))
   				case("gaussian")
					kMax(:, i) = 2*pi*corrL(:, i); !NOT REALLY IMPLEMENTED, CRITERIA TO BE DEFINED
   			end select
   		end do
    end

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

            case("gaussian")
            	eta = corrL/sqrt(pi)
            	if (nDim == 1) then
					Sk = (1./(2.*pi)) * exp(dot_product((-(eta**2)/4),(kVector**2)));
				else if (nDim == 2) then
					Sk = (1./(2.*pi)) * exp(dot_product((-(eta**2)/4),(kVector**2)));
				else if (nDim == 3) then
					Sk = (1./(2.*pi)) * exp(dot_product((-(eta**2)/4),(kVector**2)));
				else
					stop "No spectrum founded for this specifications"
				end if
        end select

        deallocate (eta)
    end

end module spectra_RF
