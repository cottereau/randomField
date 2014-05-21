module displayCarvalhol

    !All display routines

contains

    subroutine Disp1D(matrix1D)
        ! Displays 1D Matrix (Vector)
        integer :: j
        double precision, dimension(:), intent(in) :: matrix1D

        write(*,*)
        do j=1,size(matrix1D)
            write(*,*) j,"->",matrix1D(j)
        enddo
        write(*,*)
    end subroutine Disp1D

    subroutine Disp2D(matrix2D)
        ! Displays 2D Matrix, "div" columns at a time
        integer :: k, j, div=4;
        double precision :: tol = 1E-10;
        double precision, dimension(:, :), intent(in) :: matrix2D

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
