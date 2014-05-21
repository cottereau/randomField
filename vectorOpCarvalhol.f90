module vectorOpCarvalhol

    !All vector operations

contains

    function createLinearProg(qmin, qmax, nSteps) result (lpVec)
        implicit none
        double precision, intent(in) :: qmin, qmax;
        integer, intent(in) :: nSteps;
        double precision :: delta;
        double precision, dimension(:), allocatable :: lpVec;
        integer :: i;
        allocate (lpVec(nSteps))

        lpVec = (/((i-1)*(qmax-qmin)/((DBLE(nSteps))-1.0), i=1, nSteps)/)

    end

    function permut2D(matrix,mRange, direction) result (pMatrix)

        implicit none
        integer, dimension(:), intent(in) :: mRange; !components to be permuted in each row/column
        double precision, dimension(:,:), allocatable, intent(in) :: matrix; !matrix to be permuted
        character (len=1), intent(in) :: direction; !Horizontal (H) or Vertical (V)
        integer :: i, j, k, l, m, dimNumber, mSlice, mElem;
        double precision, dimension(:,:), allocatable :: pMatrix;!permuted matrix


        dimNumber = size(mRange)

        m = 1;        !factorial
        do i = 1, dimNumber
            m = m*mRange(i)
        enddo

        mSlice = 1;        !Jump between Slices
        select case(direction)

            case("V")
                allocate (pMatrix(m, size(mRange)));
                do i = dimNumber, 1, -1 !Loop Over Dimensions

                    mElem = mSlice !Jump between Elements
                    mSlice = mSlice * mRange(i) !Jump between Slices (should be = 1 before first loop)
                    do l = 1,mRange(i) !Loop over k_l elements
                        do k = 1, m/mSlice !Loop Over Slices, for permutation
                            do j = 1, mElem !Loop Over Elements of each Slice
                                pMatrix(((k-1)*mSlice+(l-1)*mElem+j),i) = matrix(i, l);
                            enddo
                        enddo

                    enddo
                enddo

            case("H")
                allocate (pMatrix(size(mRange), m));
                do i = dimNumber, 1, -1 !Loop Over Dimensions

                    mElem = mSlice !Jump between Elements
                    mSlice = mSlice * mRange(i) !Jump between Slices (should be = 1 before first loop)

                    do l = 1,mRange(i) !Loop over k_l elements
                        do k = 1, m/mSlice !Loop Over Slices, for permutation
                            do j = 1, mElem !Loop Over Elements of each Slice
                                pMatrix(i, ((k-1)*mSlice+(l-1)*mElem+j)) = matrix(i, l);
                            enddo
                        enddo

                    enddo
                enddo

        end select
    end

    function multVecEl(vector) result (prodVecEl)

        !Takes a vector and returns its component product.
        !It's intended for integer vectors like (shape(vectorND)) or a vector of the number of steps in each dimension

        implicit none
        integer, dimension(:), intent(in) :: vector;
        integer :: prodVecEl, i;

        prodVecEl = 1;

        do i = 1, size(vector)
            if (vector(i).le.0) then
                write(*,*) "WARNING: In multVecEl your vector has a strange dimension in", i, "position."
            else
                prodVecEl = prodVecEl*vector(i);
            endif
        enddo
    end

    function averageValue(matrix) result(mean)
        implicit none
        double precision, dimension(:), intent(in) :: matrix;
        double precision :: mean, nElem;
        integer :: i;

        nElem = DBLE(size(matrix));
        mean = 0;

        do i = 1, size(matrix)
            mean = mean + matrix(i);
        enddo
        mean = mean/DBLE(size(matrix))
    end
end module vectorOpCarvalhol
