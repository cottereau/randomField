module type_MESH

    use mpi

    type :: MESH
        !MPI VARIABLES
        integer :: comm
        integer :: rang
        integer :: nb_procs
        integer :: nDim;

        !MESH VARIABLES
            !nDim independent
        character (len=15) :: meshType, meshMod;
            !nDim dependent
        integer         , dimension(:)   , allocatable :: xNStep;
        double precision, dimension(:)   , allocatable :: xMax, xMin, xStep;

    end type MESH

    contains
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_MESH(MESH_a, nDim)
            type(MESH) :: MESH_a

            MESH_a%nDim = nDim
            allocate(MESH_a%xMax(nDim))
            allocate(MESH_a%xMin(nDim))
            allocate(MESH_a%xStep(nDim))
            allocate(MESH_a%xNStep(nDim))

        end subroutine

end module type_MESH
