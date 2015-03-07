module type_RF

    use mpi

    type :: RF
        !MPI VARIABLES
        integer :: comm
        integer :: rang
        integer :: nb_procs
        integer :: Nmc, nDim;

        !GENERATION VARIABLES
            !nDim independent
        character (len=15) :: corrMod, margiFirst, meshType, meshMod;
        double precision   :: fieldAvg, fieldVar;
        integer            :: method !1 for Isotropic, 2 for Shinozuka
            !nDim dependent
        double precision, dimension(:)   , allocatable :: corrL;
        double precision, dimension(:, :), allocatable :: randField, xPoints;

    end type RF

    contains
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_RF(RF_a, nDim)
            type(RF) :: RF_a

            RF_a%nDim = nDim
            allocate(RF_a%corrL(nDim))

        end subroutine

end module type_RF
