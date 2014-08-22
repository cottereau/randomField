module displayCarvalhol

	use mpi
    !All display routines
    interface DispCarvalhol
		module procedure Disp1Ddble,   &
		                 Disp2Ddble,   &
		                 Disp3Ddble,   &
		                 DispScalDble, &
		                 Disp1Dint,    &
		                 Disp2Dint,    &
		                 Disp1Dchar,   &
		                 Disp2Dchar
	end interface DispCarvalhol

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine Ordering_MPI_Start()

		implicit none

        !LOCAL VARIABLES
        integer, dimension( MPI_STATUS_SIZE ) :: statut
        integer :: rang, code, id = 15, rangCount = 0

		call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)

        if (rang /= 0) call MPI_RECV (rangCount,1, MPI_INTEGER, &
        						      rang-1, id, MPI_COMM_WORLD ,statut,code)
        write(*,*) "RANG ", rang, "-------"

    end subroutine Ordering_MPI_Start

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine Ordering_MPI_End()

    	implicit none
        !LOCAL VARIABLES
        integer :: rang, code, nb_procs, id = 15, rangCount

    	call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
    	call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

        if (rang /= nb_procs - 1) call MPI_SEND (rang + 1, 1, MPI_INTEGER , &
        			                             rang + 1, id, MPI_COMM_WORLD ,code)
    end subroutine Ordering_MPI_End

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine DispScalDble(scalar, title, format, nColumns, mpi)
        ! Displays Scalar

		implicit none

        !INPUT
        double precision,            intent(in) :: scalar
        character (len=*), optional, intent(in) :: title, format
        integer,           optional, intent(in) :: nColumns
        logical,           optional, intent(in) :: mpi

        !LOCAL VARIABLES
        double precision, dimension(:,:), allocatable :: matrix2d

		allocate(matrix2D(1,1));
		matrix2d = scalar;

		if(present(mpi)) then
			call Disp2Ddble(matrix2D, title, format, nColumns, mpi);
		else
			call Disp2Ddble(matrix2D, title, format, nColumns);
		end if

		deallocate(matrix2D);


    end subroutine DispScalDble

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp1Ddble(vector, title, format, nColumns, mpi)
        ! Displays 1D Matrix (Vector)

		implicit none

        !INPUT
        double precision,     dimension(:),          intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical, optional, intent(in) :: mpi

        !LOCAL VARIABLES
        double precision, dimension(:,:), allocatable :: matrix2d

		allocate(matrix2D(size(vector),1));
		matrix2d(:,1) = vector;

		if(present(mpi)) then
			call Disp2Ddble(matrix2D, title, format, nColumns, mpi);
		else
			call Disp2Ddble(matrix2D, title, format, nColumns);
		end if

		deallocate(matrix2D);

    end subroutine Disp1Ddble

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp1Dint(vector, title, format, nColumns, mpi)
        ! Displays 1D Matrix (Vector)

		implicit none

        !INPUT
        integer,            dimension(:),            intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical, optional, intent(in) :: mpi

        !LOCAL VARIABLES
        integer, dimension(:,:), allocatable :: matrix2D

		allocate(matrix2D(size(vector),1));
		matrix2d(:,1) = vector;

		if(present(mpi)) then
			call Disp2Dint(matrix2D, title, format, nColumns, mpi);
		else
			call Disp2Dint(matrix2D, title, format, nColumns);
		end if

		deallocate(matrix2D);

    end subroutine Disp1Dint

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp1Dchar(vector, title, format, nColumns, mpi)
        ! Displays 1D Matrix (Vector)

		implicit none

        !INPUT
        character(len=*),   dimension(:),            intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical, optional, intent(in) :: mpi

        !LOCAL VARIABLES
        character(len=30), dimension(:,:), allocatable :: matrix2D

		allocate(matrix2D(size(vector),1));
		matrix2d(:,1) = vector;

		if(present(mpi)) then
			call Disp2Dchar(matrix2D, title, format, nColumns, mpi);
		else
			call Disp2Dchar(matrix2D, title, format, nColumns);
		end if
		deallocate(matrix2D);

    end subroutine Disp1Dchar

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp2Ddble(matrix2D, title, format, nColumns, mpi)
        ! Displays 2D Matrix, "div" columns at a time

		implicit none

		!INPUT
        double precision, dimension(:, :),           intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical, optional, intent(in) :: mpi

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
		character (len=40) :: doubleFmt;
		character (len=10) :: tempFmt;

		if(present(mpi)) then
			if(mpi) call Ordering_MPI_Start()
		end if

		if(present(format))                     tempFmt = format;
		if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
		if(.not.present(format))    tempFmt = "F10.5";
		if(.not.present(nColumns))  div = 15;
		write(doubleFmt, fmt="(I3, A)") div, tempFmt;

		write(*,*) ""
		if(present(title)) write(*,*) "/////// ", title, " ///////";

		tol = 1E-10;

        write(*,*) ""
        do k=1, size(matrix2D,2)/div
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
            write(*,*)
            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
            enddo
        enddo

        if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
            write(*,*)

            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
            enddo
        end if

        write(*,*)

        if(present(mpi)) then
			if(mpi) call Ordering_MPI_End()
		end if

    end subroutine Disp2Ddble

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp2Dint(matrix2D, title, format, nColumns, mpi)
        ! Displays 2D Matrix, "div" columns at a time

		implicit none

		!INPUT
        integer,           dimension(:, :),          intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical, optional, intent(in) :: mpi

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
		character (len=40) :: doubleFmt;
		character (len=10) :: tempFmt;

		if(present(mpi)) then
			if(mpi) call Ordering_MPI_Start()
		end if

		if(present(format))                     tempFmt = format;
		if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
		if(.not.present(format))    tempFmt = "I6";
		if(.not.present(nColumns))  div = 5;
		write(doubleFmt, fmt="(I3, A)") div, tempFmt;

		write(*,*) ""
		if(present(title)) write(*,*) "/////// ", title, " ///////";

		tol = 1E-10;

        write(*,*) ""
        do k=1, size(matrix2D,2)/div
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
            write(*,*)
            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
            enddo
        enddo

        if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
            write(*,*)

            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
            enddo
        end if

        write(*,*)

        if(present(mpi)) then
			if(mpi) call Ordering_MPI_End()
		end if

    end subroutine Disp2Dint

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp2Dchar(matrix2D, title, format, nColumns, mpi)
        ! Displays 2D Matrix, "div" columns at a time

		implicit none

		!INPUT
        character (len=*), dimension(:, :),          intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical, optional, intent(in) :: mpi

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
		character (len=40) :: charFmt;
		character (len=10) :: tempFmt;

		if(present(mpi)) then
			if(mpi) call Ordering_MPI_Start()
		end if

		if(present(format))                     tempFmt = format;
		if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
		if(.not.present(format))    tempFmt = "A";
		if(.not.present(nColumns))  div = 5;
		write(charFmt, fmt="(I3, A)") div, tempFmt;

		write(*,*) ""
		if(present(title)) write(*,*) "/////// ", title, " ///////";

		tol = 1E-10;

        write(*,*) ""
        do k=1, size(matrix2D,2)/div
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
            write(*,*)
            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//charFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
            enddo
        enddo

        if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
            write(*,*)
            write(*,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
            write(*,*)

            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,fmt="(I4, A, ("//charFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
            enddo
        end if

        write(*,*)

        if(present(mpi)) then
			if(mpi) call Ordering_MPI_End()
		end if

    end subroutine Disp2Dchar

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine Disp3Ddble(matrix, title, format, nColumns, mpi)
        ! Displays 1D Matrix (Vector)

		implicit none

        !INPUT
        double precision , dimension(:,:,:),          intent(in) :: matrix
        character (len=*),                  optional, intent(in) :: title, format
        integer,                            optional, intent(in) :: nColumns
        logical, optional, intent(in) :: mpi

        !LOCAL VARIABLES
        integer :: i
		do i = 1, size(matrix, 3)
			write(*,*) "---------Slice ", i, "----------"
			if(present(mpi)) then
				call Disp2Ddble(matrix(:,:,i), title, format, nColumns, mpi);
			else
				call Disp2Ddble(matrix(:,:,i), title, format, nColumns);
			end if
		end do

    end subroutine Disp3Ddble

end module displayCarvalhol
