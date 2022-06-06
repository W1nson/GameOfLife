program gol_col

include 'mpif.h' 


integer :: ierr, myid, numprocs, request
integer :: leftid, rightid, tag 
integer :: status(MPI_STATUS_SIZE)
integer :: colid
integer, parameter:: N = 4
integer, parameter:: M = 4
integer :: cols 
integer :: left(N), right(N) 
integer :: grid(N, N)
integer, allocatable :: subgrid(:,:), updated(:,:), mycols(:,:) 

double precision :: t1, t2, time 

call MPI_INIT(ierr) 
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr) 

tag = 1234

! check if you can split up the columns evenly with the procs 
! N = 10, procs = 5 
! each gets two columns 

if (N < numprocs) then 
    cols = 1 
else  
    cols = N / numprocs
    if ((mod(N, numprocs) .ne. 0) .and. (mod(N, numprocs) <   (N / 2.0))) then 
        cols = cols + 1 
    endif 
endif 


allocate(mycols(M, cols)) 

if (myid .eq. 0) then 
    ! do j = 1,N 
    !     do i = 1,N 
    !         grid(i, j) = 0 
    !     enddo 
    ! enddo
    ! grid(2, 1) = 1
    ! grid(3, 2) = 1
    ! grid(1, 3) = 1
    ! grid(2, 3) = 1
    ! grid(3, 3) = 1
    
    call fill_rand(grid, M, N) 
        
    ! call print_grid(grid, M, N)


endif


t1 = MPI_Wtime() 


! scatter each col to each processor
call MPI_SCATTER(grid, cols * N, MPI_INTEGER, mycols, cols *  N, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
! scatter each row to each processor 
! call MPI_SCATTERV(grid, N, M, MPI_INTEGER, myrows, N, MPI_INTEGER, 

! scatter 
colid = myid*cols + 1
rowid = myid + 1

if (colid <= N) then 
     
    ! print *, 'processor ', myid, "col: ", colid
    if(myid .eq. numprocs-1) then
        cols = N - colid + 1
        ! mycols = mycols(:, 1: N-colid+1)
        ! call print_grid(mycols, M, cols) 
    else  
        ! call print_grid(mycols, M, cols)
    endif
endif

allocate(subgrid(M+2, cols+2)) 
allocate(updated(M+2, cols+2))
subgrid = 0 
updated = 0 
! print *, mycol
! print *

! send the column to the left and right 
! check if the col is edge 

if (myid .eq. 0) then 
    leftid  = numprocs - 1
    rightid = myid + 1 
else if (myid .eq. numprocs-1) then 
    leftid = myid -1 
    rightid = 0 
else 
    leftid = myid -1 
    rightid = myid + 1
endif 


! iterations  
do iter = 1, 1000 
   
    ! send the right most column to the right proc
    call MPI_ISEND(mycols(:, cols), M, MPI_INTEGER, rightid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! receive the right most colum from the left proc
    call MPI_IRECV(left, M, MPI_INTEGER, leftid, tag, MPI_COMM_WORLD,request, ierr)
    
    ! waiting for each proc to complete
    call MPI_WAIT(request, status, ierr)
    
    ! send the left most column to the right proc
    call MPI_ISEND(mycols(:, 1), M, MPI_INTEGER, leftid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! receive the left most column from the right proc 
    call MPI_IRECV(right, M, MPI_INTEGER, rightid, tag, MPI_COMM_WORLD,request, ierr)
    
    ! waiting for each proc to complete
    call MPI_WAIT(request, status, ierr)
   
    
    ! mesh the left and right column with the center columns 
    do i = 2, N+1 
        subgrid(i, :) = [left(i-1), mycols(i-1, :), right(i-1)]
    enddo 
    
    ! print *, 'processor ', myid, "col: ", colid

    ! wrapping top and bottom 
    subgrid(1, :) = subgrid(N+1,:) 
    subgrid(N+2, :) = subgrid(2,:)
    ! call print_grid(subgrid, M+2,  cols+2)
    
    ! evolution updates
    do i = 2, cols+1 
        call evo_col(subgrid, updated, M, cols, i) 
        ! call print_grid(updated, M+2, cols+2)
        mycols(:, i-1) = updated(2:N+1, i) 
    enddo 
    ! call print_grid(mycols, M, cols)
enddo 

! gather all the columns from other procs back to proc 0 matrix 
call MPI_GATHER(mycols, cols*N, MPI_INTEGER, grid, cols*N, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 

t2 = MPI_Wtime()

if(myid .eq. 0) then 
    print *, "iteration: ", iter  
    print *,  "time: ", (t2-t1)
    ! call print_grid(grid, M, N) 
endif 

deallocate(mycols)
deallocate(subgrid) 
deallocate(updated) 

call MPI_FINALIZE(ierr) 

stop
end program gol_col 









! subroutine insert(m, b, col)
!     integer, allocatable, intent(inout) :: m(:, :)
!     integer, intent(in) :: b(size(m, 1)), col
!     integer, allocatable :: temp(:, :)
!     integer :: rows, cols
!     rows = size(m, 1)
!     cols = size(m, 2)
!     allocate(temp(rows, cols + 1))
!     temp(:, 1:col) = m(:, 1:col)
!     temp(:, col) = b
!     temp(:, col + 1:cols + 1) = m(:,col:cols)
!     call move_alloc(temp, m)
! end subroutine insert
! function  wrap_top_bot(A, N) result(B)
!     integer :: N  
!     integer :: A(N, 3) 
!     integer :: B(N+2, 3) 
!     
!     do i = 2, N+1 
!         B(i, :) = A(i-1,:)
!     enddo 
!     B(1,:) = A(N, :) 
!     B(N+2, :) = A(1,:)
! 
! endfunction 

subroutine evo_col(A, B, M, N, col) 
    implicit none 
    
    integer :: col, M, N, c, k, l, row
    integer :: A(M+2, N+2) 
    integer :: B(M+2, N+2) 

     
    do row = 2, M+1 
        c = -A(row, col)
        do k = col-1, col+1
            do l = row-1, row+1 
                c = c + A(l,k) 
            enddo 
        enddo 
        if(c == 3) then 
            B(row, col) = 1
        else if(c==2) then 
            B(row, col) = A(row, col) 
        else 
            B(row, col) = 0
        endif    
    enddo

end subroutine evo_col


subroutine print_grid(A, M, N) 
    implicit none 
    integer :: N, M, i 
    integer :: A(M, N) 
    integer :: c 
    c = 0 

    do i = 1, M
        print *, A(i, :) 
        c = c+SUM(A(i,:))
    enddo 
    print * 
    print *, "totoal alive", c 
    print *
    print *

end subroutine print_grid


subroutine fill_rand(A, M, N)
    implicit none

    integer :: M, N, i, j, r
    integer :: A(M, N)

    do j = 1, N
        do i = 1, M
            call rnum(r)
            A(i, j) = r
        enddo
    enddo
end subroutine fill_rand


subroutine rnum(r)
    implicit none

    integer :: r
    r = floor(2 * rand())

end subroutine rnum
