module spectra_RF

contains


    function make_Spectrum1D(kVector, corrMod, corrL) result (Sk)
        ! Return Spectrum from a chosen correlation model
        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: kVector;
        character (len=20), intent(in) :: corrMod
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
end module spectra_RF
