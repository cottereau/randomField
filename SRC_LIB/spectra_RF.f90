module spectra_RF
    use displayCarvalhol
    use math_RF
    use write_Log_File
    use constants_RF
    use type_RF

    implicit none
contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kMaxND(corrMod, kMax);
        implicit none
        !INPUT
        character (len=*),                intent(in) :: corrMod;
        !double precision,   dimension(:),  intent(in), optional :: corrL;

        !OUTPUT
        double precision, dimension(:),   intent(out) :: kMax;

        !LOCAL VARIABLES
        double precision :: pi = 3.1415926535898
        integer          :: i, nDim

        nDim = size(kMax)

        select case(trim(adjustL(corrMod)))
            case("gaussian")
                select case(nDim)
                    case(1)
                        kMax(:) = 6.457D0; !Value to cover 99% Spectra area
                    case(2)
                        kMax(:) = 7.035D0; !Value to cover 99% Spectra area
                    case(3)
                        kMax(:) = 7.355D0; !Value to cover 99% Spectra area
                end select
        end select

        write(get_fileId(),*) "kMax = ", kMax


    end subroutine set_kMaxND

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kPoints(RDF);
        implicit none

        !INPUT OUTPUT
        type(RF) :: RDF

        !LOCAL
        integer :: i
        double precision :: kAdjust    = 1.0D0 !"kNStep minimum" multiplier
        double precision :: periodMult = 1.1D0 !"range" multiplier
        double precision :: rAdjust    = 1.0D0 !"rNStep minimum" multiplier

        if(allocated(RDF%kPoints)) deallocate(RDF%kPoints)

        call set_kMaxND(RDF%corrMod, RDF%kMax) !Defining kMax according to corrMod

        if(RDF%independent) then
            RDF%kDelta(:) = 2.0D0*PI/(periodMult*(RDF%xMaxBound - RDF%xMinBound)) !Delta max in between two wave numbers to avoid periodicity
        else
            RDF%kDelta(:) = 2.0D0*PI/(periodMult*(RDF%xMaxGlob - RDF%xMinGlob)) !Delta max in between two wave numbers to avoid periodicity
            !kDelta(:,1) = 2.0D0*PI/(periodMult*(RDF%xMaxBound - RDF%xMinBound))
        end if

        select case (RDF%method)
            case(ISOTROPIC)

            case(SHINOZUKA)
                RDF%kNStep(:)   = 1 + kAdjust*(ceiling(RDF%kMax/RDF%kDelta(:))); !Number of points in k
                RDF%kDelta(:) = (RDF%kMax)/(RDF%kNStep-1); !Redefining kDelta after ceiling and adjust
                RDF%kNTotal = product(RDF%kNStep);

                write(get_fileId(),*) "RDF%kNStep = ", RDF%kNStep
                write(get_fileId(),*) "RDF%kNTotal = ", RDF%kNTotal

                allocate(RDF%kPoints(RDF%nDim, RDF%kNTotal))

                do i = 1, RDF%kNTotal
                    call get_Permutation(i, RDF%kMax, RDF%kNStep, RDF%kPoints(:, i), snapExtremes = .true.);
                end do

            case(RANDOMIZATION)
                RDF%kNStep(:)   = 1 + kAdjust*(ceiling(RDF%kMax/RDF%kDelta(:))); !Number of points in k
                RDF%kDelta(:) = (RDF%kMax)/(RDF%kNStep-1); !Redefining kDelta after ceiling and adjust
                RDF%kNTotal = product(RDF%kNStep);

                allocate(RDF%kPoints(RDF%nDim, RDF%kNTotal))
                call random_number(RDF%kPoints(:,:))
                do i = 1, RDF%nDim
                    RDF%kPoints(i,:) = RDF%kPoints(i,:) * RDF%kMax(i)
                end do

            case(FFT)
                RDF%kNStep(:) = 2*RDF%xNStep(:); !Number of points in k
                RDF%kNTotal   = product(RDF%kNStep);
                RDF%kMax(:)   = dble(RDF%kNStep(:) - 1) * RDF%kDelta(:)

        end select

    end subroutine set_kPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_SkVec(RDF);
        implicit none

        !OBS: corrL is supposed = 1. The complete formula is RDF%SkVec = product(corrL) * exp(-dot_product(kVector**2, corrL_effec**2)/(4.0d0*pi))

        !INPUT OUTPUT
        type(RF) :: RDF

        !LOCAL
        integer :: i, j, k
        logical :: cycleX, cycleY, cycleZ


        if(allocated(RDF%SkVec)) deallocate(RDF%SkVec)
        if(RDF%method /= FFT ) allocate(RDF%SkVec(RDF%kNTotal))

        select case(RDF%corrMod)

            case("gaussian")
                select case (RDF%method)
                    case(FFT)
                        if(RDF%nDim == 3 ) then

                            allocate(RDF%Sk3D(RDF%kNStep(1), RDF%kNStep(2), RDF%kNStep(3)))
                            RDF%Sk3D(:,:,:) = 0.0D0;

                            do k = 1, RDF%kNStep(3)
                                do j = 1, RDF%kNStep(2)
                                    do i = 1, RDF%kNStep(1)

                                        RDF%Sk3D(i,j,k) = exp(                                  &
                                            -(((i-1)*RDF%kDelta(1))**(2.0D0)  &
                                            + ((j-1)*RDF%kDelta(2))**(2.0D0)  &
                                            + ((k-1)*RDF%kDelta(3))**(2.0D0)) &
                                            /(4.0d0*pi))

                                        !TODO: Optimize to cycle when reaching Sk3D < TOLERANCE)

!                                        if(i==1) write(*,*) "RDF%Sk3D(i,j,k) = ",  RDF%Sk3D(i,j,k);
!                                        if(i==1) write(*,*) "sum k**2 = ",  ((i-1)*RDF%kDelta(1))**(2.0D0)  &
!                                                                          + ((j-1)*RDF%kDelta(2))**(2.0D0)  &
!                                                                          + ((k-1)*RDF%kDelta(3))**(2.0D0)
                                    end do

                                end do
                            end do
                        else
                            write(*,*) "ERROR!, FFT only is implemented in 3D)"
                            stop
                        end if

                    case default
                        RDF%SkVec = exp(-sum(RDF%kPoints**(2.0D0), 1)/(4.0d0*pi))

                end select
        end select

    end subroutine set_SkVec

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_rMax(corrMod, rMax, corrL);
        implicit none
        !INPUT
        character (len=*),                intent(in) :: corrMod;
        double precision,   dimension(:),  intent(in), optional :: corrL;

        !OUTPUT
        double precision, dimension(:),   intent(out) :: rMax;

        !LOCAL VARIABLES
        double precision :: pi = 3.1415926535898
        integer          :: i
        double precision, dimension(:), allocatable:: corrL_effec

        allocate(corrL_effec (size(rMax)))

        if (present(corrL)) corrL_effec = corrL
        if (.not. present(corrL)) corrL_effec = 1
        !        do i = 1, 100
        !            kMax = i/10.0 * corrL(:)
        !            write(*,*) "kMax = ", kMax
        !            write(*,*) "Spectrum = ", get_SpectrumND(kMax, corrMod, corrL)
        !            call DispCarvalhol (kMax, "kMax")
        !            call DispCarvalhol (get_SpectrumND(kMax, corrMod, corrL), "Spectrum")
        !        end do

        select case(trim(adjustL(corrMod)))
        case("gaussian")
            rMax(:) = 2*pi*corrL_effec(:); !CRITERIA STILL TO BE TESTED
        end select

        deallocate(corrL_effec)

    end subroutine set_rMax

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kSign(kSign);
        implicit none
        !INPUT and OUTPUT
        double precision,   dimension(:, :), intent(inout) :: kSign;

        !LOCAL VARIABLES
        integer:: i, j, pos, seedStep, nDim;

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

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kArray_rand(corrMod, kArray_rand, kMaxScal, nDim);
        implicit none
        !INPUT
        character (len=*), intent(in) :: corrMod;
        double precision, optional :: kMaxScal
        integer :: nDim

        !OUTPUT
        double precision, dimension(:),  intent(out) :: kArray_rand;

        !LOCAL VARIABLES
        integer :: i, j
        double precision ::  bound,mean,cdf_x,q,sd,x, uniRand
        integer :: st, which
        real :: av = 0.0
        real :: gennor
        real :: std = 1.0
        real :: normRand

!        mean = 0.0D0
!        sd = 1.0D0
!        test2 = gennor (real(mean), real(sd))
!        write(*,*) "gennor (0.0, 1.0) = ", test2

!        call random_number(kArray_rand)
!        kArray_rand = kMaxScal*kArray_rand

!TODO RANDOMIZATION
        which = 2 !Find x, for a given cdf_x

        !write(*,*) "corrMod = ", corrMod

        select case(trim(adjustL(corrMod)))
            case("gaussian")
                !write(*,*) "Gaussian kArray"
                mean = 0.0D0
                sd = SQRT_2PI
                !sd = 1.0D0/SQRT_2PI
                do i = 1, size(kArray_rand)
                    call random_number(uniRand)
                    cdf_x = (0.5D0 + uniRand/2.0D0)
                    q = 1-cdf_x
                    call cdfnor(which,cdf_x,q,kArray_rand(i),mean,sd,st,bound) !Normal CDF
                    kArray_rand = kArray_rand**(dble(nDim))
                end do

!                sd = SQRT_2PI
!                !sd = 1.0D0/SQRT_2PI
!                do i = 1, size(kArray_rand)
!                    cdf_x = 1.0D0
!                    do j = 1, nDim
!                        call random_number(uniRand)
!                        cdf_x = cdf_x*(0.5D0 + uniRand/2.0D0)
!                    end do
!                    q = 1-cdf_x
!                    call cdfnor(which,cdf_x,q,kArray_rand(i),mean,sd,st,bound) !Normal CDF
!                end do


                !write(*,*) "gennor (mean, sd) = ",gennor (0.5, 1.0)
                ! PROBLEM, how to put the seed
!                mean = 0.0D0
!                sd = 1.0D0/SQRT_2PI
!
!                kArray_rand(:) = 1
!                do i = 1, size(kArray_rand)
!                    do j = 1, nDim
!                        normRand = gennor (real(mean), real(sd))
!                        kArray_rand(i) = kArray_rand(i)*abs(dble(normRand))
!                    end do
!                end do
        end select



        !call reorder_vector(kArray_rand)

        !kArray_rand = kArray_rand**(dble(nDim))

        !call dispCarvalhol(kArray, "kArray")

    end subroutine set_kArray_rand

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kDelta_rand(kArray_rand, kDelta_rand);
        implicit none
        !INPUT
        double precision, dimension(:),  intent(in) :: kArray_rand;

        !OUTPUT
        double precision, dimension(:),  intent(out) :: kDelta_rand;

        !LOCAL VARIABLES
        integer :: i
        double precision :: inc

        kDelta_rand = 0.0D0

        !Correction for the first term
        inc = (kArray_rand(1) - 0.0D0)
        kDelta_rand(1) = kDelta_rand(1) + inc

        !Calculating distance between numbers (we supose they're ordered)
        do i = 1, size(kArray_rand)-1
            inc = (kArray_rand(i+1) - kArray_rand(i))/2.0D0
            kDelta_rand(i)   = kDelta_rand(i)   + inc
            kDelta_rand(i+1) = kDelta_rand(i+1) + inc
        end do

    end subroutine set_kDelta_rand

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function get_SpectrumND(kVector, corrMod, corrL) result (Sk)
        ! Return Spectrum from a chosen correlation model
        implicit none

        !INPUT
        double precision, dimension(:), intent(in) :: kVector;
        character (len=*),              intent(in) :: corrMod
        double precision, dimension(:), intent(in), optional :: corrL;

        !OUTPUT
        double precision :: Sk;

        !LOCAL VARIABLES
        integer :: j, nDim;
        double precision, dimension(:), allocatable :: eta;
        double precision :: pi = 3.1415926535898
        double precision, dimension(:), allocatable:: corrL_effec

        allocate(corrL_effec (size(kVector)))

        if (present(corrL)) corrL_effec = corrL
        if (.not. present(corrL)) corrL_effec = 1

        Sk = 0;
        nDim = size(kVector)

        select case(corrMod)
        case("gaussian")

            !REGIS
            !Sk  = exp(-dot_product((kVector**2),(corrL_effec**2))/(4.0d0)); !Amplitude part "product(corrL)" is external to the function

            !MEU
            Sk = exp(-dot_product(kVector**2, corrL_effec**2)/(4.0d0*pi)); !Amplitude part "product(corrL)" is external to the function
            !write(*,*) "Sk = ", Sk
            !write(*,*) "kVector = ", kVector

        end select

        deallocate(corrL_effec)

    end function get_SpectrumND

end module spectra_RF
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
