module displayCarvalhol

    !All display routines

contains

    subroutine Disp1D(matrix1D, name)
        ! Displays 1D Matrix (Vector)

        !INPUT
        character (len=*), optional, intent(in) :: name
        double precision, dimension(:), intent(in) :: matrix1D

        !LOCAL VARIABLES
        integer :: j

		if(present(name)) write(*,*) name;

        write(*,*)
        do j=1,size(matrix1D)
            write(*,*) j,"->",matrix1D(j)
        enddo
        write(*,*)
    end subroutine Disp1D

    subroutine Disp2D(matrix2D, name, nColumns, tolerance)
        ! Displays 2D Matrix, "div" columns at a time

		!INPUT
        double precision, dimension(:, :), intent(in) :: matrix2D
        character (len=*), optional, intent(in) :: name
        integer, optional, intent(in) :: nColumns
        double precision, optional, intent(in) :: tolerance;

        !LOCAL VARIABLES
        integer :: k, j, div;
        double precision :: tol;


		if(present(name)) write(*,*) name;
		if(present(nColumns).and.nColumns.gt.0) div = nColumns;
		if(present(tolerance).and.tolerance.gt.0d0) tol = tolerance;

		if(.not.present(nColumns)) div = 4;
		if(.not.present(tolerance)) tol = 1E-10;

        write(*,*)
        do k=1, size(matrix2D,2)/div
            write(*,*)
            write(*,*) "Columns", (k-1)*div+1, " to ", k*div ;
            write(*,*)
            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,*) j,"->",matrix2D(j,(k-1)*div+1:k*div)
            enddo
        enddo

        if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
            write(*,*)
            write(*,*) "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
            write(*,*)

            do j= lbound(matrix2D,1), ubound(matrix2D,1)
                write(*,*) j,"->",matrix2D(j,(k-1)*div+1:)
            enddo
        end if

        write(*,*)
    end subroutine Disp2D

end module displayCarvalhol
