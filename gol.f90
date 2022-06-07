program gol

include 'mpif.h' 

integer :: ierr, myid, numprocs, request
integer :: leftid, rightid, topid, botid, diagid, topleftid, toprightid, botleftid,botrightid, tag 
integer :: colid, rowid 
integer :: status(MPI_STATUS_SIZE)
integer, parameter:: M = 20
integer, parameter:: N = 20
integer :: cols, rows, col, row 
integer :: tl, tr, bl, br 
integer, allocatable :: left(:), right(:), top(:), bot(:)
integer, allocatable :: displace(:)
integer :: grid(M, N), temp(M,N)
integer :: mycorners(4), corners(4)
integer, allocatable :: subgrid(:,:), updated(:,:)

integer :: gridprocs

! change the number inside to how many outputs you want to see throughout the iterations 
integer, parameter :: outs = 4
integer :: output_iter(outs) 
integer :: iter_out(outs, M, N)
integer :: c 

integer, allocatable:: mysubmat(:,:)
integer :: submat, resizedtype
integer :: sizes(2), blockDim(2), starts(2) 
integer(kind=MPI_ADDRESS_KIND) :: extent
integer :: intsize

integer, allocatable:: counts(:) 
double precision :: t1, t2, time 

call MPI_INIT(ierr) 
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr) 



output_iter = (/0, 20, 40, 80/)


! check if you can split up the columns and rows evenly with the procs 
! M = 4, N = 4, procs = 4 
! each gets two columns 

if(floor(sqrt(real(numprocs)))**2 - numprocs .ne. 0) then
    if(myid .eq. 0)  then
        print *, "Your number of proccesors doesn't make a square grid, change number of nodes"
    endif
    call MPI_FINALIZE(ierr)
    stop
endif 

gridprocs = sqrt(real(numprocs))
cols = N / sqrt(real(numprocs)) 
rows = M / sqrt(real(numprocs))

if(cols * sqrt(real(numprocs)) .ne. N) then 
    if(myid .eq. 0) then 
        print *, "Your grid doesn't divide evenly, change size of your grid" 
    endif 
    call MPI_FINALIZE(ierr) 
    stop
endif 
if(rows * sqrt(real(numprocs)) .ne. M) then 
    if(myid .eq. 0) then 
        print *, "Your grid doesn't divide evenly, change size of your grid" 
    endif 
    call MPI_FINALIZE(ierr) 
    stop
endif

c = 1
! initialize data 
if (myid .eq. 0) then


    print *, rows, cols

    do j = 1,N 
        do i = 1,N 
            grid(i, j) = 0 
        enddo 
    enddo
    grid(2, 1) = 1
    grid(3, 2) = 1
    grid(1, 3) = 1
    grid(2, 3) = 1
    grid(3, 3) = 1

    ! call fill_rand(grid, M, N) 
        
    call print_grid(grid, M, N)
    if(output_iter(c) == 0) then 
        iter_out(c,:,:) = grid
        c = c+1
    endif
        
endif



! 
! if (N < numprocs) then 
!     cols = 1
! 
! else  
!     cols = N / numprocs
!     if ((mod(N, numprocs) .ne. 0) .and. (mod(N, numprocs) <   (N / 2.0))) then 
!         cols = cols + 1 
!     endif 
! endif 


sizes = (/M, N/) 
blockDim = (/rows, cols/)
starts = (/0, 0/)

! create submat
call MPI_Type_create_subarray(2, sizes, blockDim, starts, MPI_ORDER_FORTRAN, MPI_INTEGER,submat ,ierr)

call MPI_Type_size(MPI_INTEGER, intsize, ierr) 
extent = rows * intsize 
call MPI_Type_create_resized(submat, 0, extent, resizedtype, ierr) 

call MPI_Type_commit(resizedtype, ierr)

! counts is how many submatrix needs to send
allocate(counts(0:numprocs-1))
allocate(displace(numprocs))

if (myid .eq. 0) then 
    
    counts = 1 
    forall( col=1:gridprocs, row=1:gridprocs)
        displace(1+(row-1)+gridprocs*(col-1)) = (row-1) + rows*gridprocs*(col-1)
    endforall
    ! setup displacement   
    ! displace = (/0, 1, 8, 10/)


    print *, "counts" 
    print *, counts
    print *  
    print *, "displace" 
    print *, displace
    print * 
endif
allocate(mysubmat(rows, cols))
! scatter each row to each processor 
call MPI_SCATTERV(grid, counts, displace, resizedtype, mysubmat, rows*cols, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr);

! print *, "myid", myid
! print *
! call print_grid(mysubmat, rows, cols) 

allocate(subgrid(rows+2, cols+2)) 
allocate(updated(rows+2, cols+2))
subgrid = 0 
updated = 0 

allocate(left(rows)) 
allocate(right(rows))
allocate(top(cols))
allocate(bot(cols))


t1 = MPI_Wtime() 
! iterations  
do iter = 1, 80 
     
    tag = 1234 
    ! find the left and right id 
    ! send the column to the left and right 
    ! check if the col is edge
    if (myid + gridprocs >= gridprocs**2) then 
        leftid  = mod(myid + gridprocs, gridprocs)
        rightid = mod(myid + gridprocs, gridprocs)
    else 
        leftid  = myid + gridprocs
        rightid = myid + gridprocs
    endif 

    ! print *, "myid", myid, "left", leftid, "right", rightid 
    ! print *, "my left", mysubmat(:,1) 
    ! print *, "my right", mysubmat(:,cols)

    ! send the right most column to the right proc
    call MPI_ISEND(mysubmat(:, cols), rows, MPI_INTEGER, rightid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! receive the right most colum from the left proc
    call MPI_IRECV(left, rows, MPI_INTEGER, leftid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! waiting for each proc to complete
    call MPI_WAIT(request, status, ierr)
    
    ! send the left most column to the right proc
    call MPI_ISEND(mysubmat(:, 1), rows, MPI_INTEGER, leftid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! receive the left most column from the right proc 
    call MPI_IRECV(right, rows, MPI_INTEGER, rightid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! waiting for each proc to complete
    call MPI_WAIT(request, status, ierr)
    
    ! print *, "left: ", left 
    ! print *, "right: ", right
   
    ! only works for 4 procs 2x2 grid 
    if (mod(myid, gridprocs) == 0) then 
        topid  = myid + 1
        botid = myid + 1
    else 
        topid  = myid - 1
        botid = myid - 1
    endif 

    ! print * 
    ! print *, "myid", myid, "top", topid, "bot", botid 

    ! print *, "my top", mysubmat(1,:) 
    ! print *, "my bot", mysubmat(rows,:)

    ! send the bot most row to the bot proc
    call MPI_ISEND(mysubmat(rows, :), cols, MPI_INTEGER, botid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! receive the bot  most colum from the top proc
    call MPI_IRECV(top, cols, MPI_INTEGER, topid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! waiting for each proc to complete
    call MPI_WAIT(request, status, ierr)
    
    ! send the left most column to the right proc
    call MPI_ISEND(mysubmat(1, :), cols, MPI_INTEGER, topid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! receive the left most column from the right proc 
    call MPI_IRECV(bot, cols, MPI_INTEGER, botid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! waiting for each proc to complete
    call MPI_WAIT(request, status, ierr)

    ! print *, "top: ", top 
    ! print *, "bot: ", bot


    ! sending the diagnals 
    ! 2x2 gridprocs only!!
    if(mod(myid, gridprocs) == 0) then 
        diagid = myid-1 
        if(diagid < 0) then 
            diagid = numprocs-1
        endif
    else 
        diagid = myid+1 
        if(diagid >= numprocs) then 
            diagid = 0
        endif
    endif
    
    
    ! send all the corners of diagnals 
    mycorners(1) = mysubmat(1, 1) !tl 
    mycorners(2) = mysubmat(1, cols) !tr
    mycorners(3) = mysubmat(rows, 1) !bl
    mycorners(4) = mysubmat(rows, cols) !br
    
    call MPI_ISEND(mycorners, 4, MPI_INTEGER, diagid, tag, MPI_COMM_WORLD, request, ierr)
    
    call MPI_IRECV(corners, 4, MPI_INTEGER, diagid, tag, MPI_COMM_WORLD, request, ierr)
    
    call MPI_WAIT(request, status, ierr)
    
    
    tl = corners(4) 
    tr = corners(3) 
    bl = corners(2) 
    br = corners(1)
    

    ! mesh the left and right column with the center columns
    
    do j = 2, cols+2-1 
        do i = 2, rows+2-1 
            subgrid(i, j) = mysubmat(i-1, j-1) 
        enddo 
    enddo
    subgrid(2:rows+2-1, 1) = left 
    subgrid(2:rows+2-1, cols+2) = right 
    subgrid(1, 2:cols+2-1) = top 
    subgrid(rows+2, 2:cols+2-1) = bot 
   
    subgrid(1,1) = tl 
    subgrid(1,cols+2) = tr
    subgrid(rows+2, 1) = bl
    subgrid(rows+2, cols+2) = br
     

    ! print *, myid, "subgrid"
    ! call print_grid(subgrid, rows+2, cols+2) 

    ! evolution update 
    do j = 2, cols+1 
        do i = 2, rows+1
            call evo(subgrid, updated, rows+2, cols+2, i, j) 
            mysubmat(i-1, j-1) = updated(i, j) 
        enddo 
    enddo

    if(myid .eq. 0) then 
        print *, iter, output_iter(c)
    endif
    !if(iter == output_iter(c)) then 
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        ! print *, "gather"
    call MPI_Gatherv(mysubmat, cols**2, MPI_INTEGER, grid, counts, displace, resizedtype, 0, MPI_COMM_WORLD, ierr) 
    if(iter == output_iter(c)) then 
        iter_out(c, :,:) = grid
        !if(myid .eq. 0) then 
        !     print *, iter
            ! call print_grid(grid, M, N)
            !print *, "load"
            ! iter_out(c,:,:) = temp
        ! endif
        c = c+1
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
enddo 


! gather everything back 

call MPI_Gatherv(mysubmat, cols**2, MPI_INTEGER, grid, counts, displace, resizedtype, 0, MPI_COMM_WORLD, ierr) 
t2 = MPI_Wtime()

if(myid .eq. 0) then 
    print *, "iteration: ", iter-1 
    print *,  "time: ", (t2-t1)
    do i = 1, outs
        print *, "iteration:", output_iter(i)  
        call print_grid(iter_out(i,:,:), M, N)
    enddo 
    ! call print_grid(grid, M, N)
endif 

deallocate(subgrid) 
deallocate(updated)
deallocate(left) 
deallocate(right)
deallocate(top)
deallocate(bot)
call MPI_TYPE_FREE(resizedtype, ierr)
call MPI_FINALIZE(ierr) 

stop
end program gol 


subroutine sendCol(myid, numprocs, mysubmat, numcols, numrows, left, right, status, request, ierr)

    integer :: myid, numprocs, numcols, numrows, stat_size 
    integer :: leftid, rightid 
    integer :: tag  
    integer :: mysubmat(numrows, numcols)
    integer, intent(out) :: left(numrows), right(numrows)
    integer :: ierr, request
    ! integer :: status(stat_size)
    

    tag = 1234 
    ! find the left and right id 
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


   ! send the right most column to the right proc
    call MPI_ISEND(mysubmat(:, numcols), numrows, MPI_INTEGER, rightid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! receive the right most colum from the left proc
    call MPI_IRECV(left, numrows, MPI_INTEGER, leftid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! waiting for each proc to complete
    call MPI_WAIT(request, MPI_STATUS_IGNORE, ierr)
    
    ! send the left most column to the right proc
    call MPI_ISEND(mysubmat(:, 1), numrows, MPI_INTEGER, leftid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! receive the left most column from the right proc 
    call MPI_IRECV(right, numrows, MPI_INTEGER, rightid, tag, MPI_COMM_WORLD, request, ierr)
    
    ! waiting for each proc to complete
    call MPI_WAIT(request, MPI_STATUS_IGNORE, ierr)
end subroutine sendCol





subroutine evo(A, B, M, N, row, col) 
    implicit none 
    integer :: col, row, M, N, c, l, k
    integer :: A(M, N) 
    integer :: B(M, N) 
    c = -A(row, col) 
    do k = col-1, col+1 
        do l = row-1, row+1
            c = c + A(l,k) 
        enddo 
    enddo 
    if(c == 3) then 
        B(row, col) = 1
    else if (c==2) then 
        B(row, col) =  A(row, col) 
    else 
        B(row, col) = 0 
    endif 
end subroutine evo


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
    integer :: N, M, i,j 
    integer :: A(M, N) 
    integer :: c 
    c = 0 
     
    do j = 1, N
        do i = 1, M
            write (*, '(i2))', ADVANCE="NO"), A(i, j)
        enddo
        print *
    enddo 
    do i = 1, M
        c = c+SUM(A(i,:))
    enddo 
    print * 
    print *, "total alive", c 
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
