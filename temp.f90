program temp 
    
    implicit none

    integer :: N, i, j, k, l, c, r
    
    integer, allocatable :: grid(:,:), updated(:,:) 
    
    N = 4
    call srand(3) 
    allocate(grid(N+2,N+2))
    allocate(updated(N+2,N+2))
    
    ! call fill_rand(grid, N) 
    
    grid(3, 2) = 1 
    grid(4, 3) = 1
    grid(2, 4) = 1 
    grid(3, 4) = 1 
    grid(4, 4) = 1
    

    do i = 1, 3
        print *, "iteration: ",i
        print *
        call wrap(grid, N) 
        call print_grid(grid, N)
    
        call evo(grid, updated, N)

        ! call wrap(updated, N)
        ! call print_grid(updated, N) 
        grid = updated
    enddo

    deallocate(grid)
    deallocate(updated)


end program temp 

subroutine fill_rand(A, N) 
    implicit none 
    
    integer :: N, i, j, r
    integer :: A(N+2, N+2) 
    
    do j = 2, N+1 
        do i = 2, N+1 
            call rnum(r)
            A(i, j) = r
        enddo 
    enddo
end subroutine fill_rand


subroutine evo(A, B, N)
    implicit none 

    integer :: i, j, k, l, c, N
    integer :: A(N+2, N+2) 
    integer :: B(N+2, N+2) 

    do j = 2, N+1
        do i = 2, N+1 
            c = -A(i,j)
            do k = j-1, j+1 
                do l = i-1, i+1 
                    c = c + A(l,k) 
                enddo
            enddo
            if(c == 3) then 
                B(i,j) = 1
            else if(c == 2) then 
                B(i,j) = A(i,j)
            else 
                B(i,j) = 0
            endif
        enddo
    enddo 


end subroutine evo


subroutine rnum(r) 
    implicit none 

    integer :: r
    r = floor(2 * rand()) 

end subroutine rnum


subroutine print_grid(A, N) 
    implicit none 
    integer :: N, i
    integer :: A(N+2,N+2)
    integer :: c
    c = 0
    do i = 2, N+1 
        print *, A(i, 2:N+1)
        c = c + SUM(A(i, 2:N+1))
    enddo
    print *
    print *, "total alive", c 
    print * 
    ! do i = 1, N+2  
    !    print *, A(i,:)
    ! enddo  
    print *
    print *

end subroutine print_grid

subroutine wrap(A, N) 
    implicit none 
    integer :: N
    integer :: A(N+2, N+2) 

    A(1, :) = A(N+1, :) 
    A(N+2, :) = A(2,:) 
    A(:, 1) = A(:, N+1) 
    A(:, N+2) = A(:,2) 
    
    A(1, 1) = A(N+1,N+1) 
    A(N+2,N+2) = A(2,2) 
    A(N+2, 1) = A(2, N+1) 
    A(1, N+2) = A(N+1, 2)

end subroutine wrap 

