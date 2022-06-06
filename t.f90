program t 


integer :: temp(4,4) 

temp(1,1) = 1
temp(1,2) = 2 
temp(3,2) = 40

do i = 1, 4 
    print *, temp(i,:) 
enddo


end program t 
