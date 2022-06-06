program test

implicit none 

integer, dimension (5) :: b = (/ 21, 22, 23, 24, 25 /)
integer, dimension (5, 1) :: c, d, e 
integer :: i  

integer, dimension (7, 3):: subgrid, temp 

c = reshape( b, (/ 5, 1/) )
d = c
e = c 

do i = 1, 5
    print *, d(i, :) 
enddo 

do i = 2, 6 
    subgrid(i,:) = [c(i-1,:), d(i-1,:), e(i-1,:)] 
enddo 

do i = 1, 7 
    print *, subgrid(i, :)
enddo 


temp = wrap_top_bot(subgrid, 7) 

do i = 1, 7
    print *, temp(i,:)
enddo



stop 
end

function  wrap_top_bot(A, N) result(B)
    implicit none
    integer :: N, i
    integer :: A(N+2, 3)
    integer :: B(N+2, 3)
    B(1,:) = A(N-1, :)
    B(N+2, :) = A(2,:)
endfunction
