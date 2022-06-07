
# MPI GameOfLife 

## PCAM process: 

* Partitioning: 
	
	I have decomposed the computation for each task is to tell whether a cell will be alive or not. 

* Communication: 

	Since the each cell needs its 8 nearing neighbors to compute, there will have 8 communication if each node takes each cell. 

* Agglomeration: 

	For column decomposition, I join the cells into columns, which creates a pencil shape for each node to compute. This leads to that each node will only need to communite to 2 neighboring nodes to get the left and right neighbors. 

	For 2D decomposition, I join the cells into small blocks, which creates combines both rows and columns of cells into small blocks. This can lead to that each node has to communicate with 8 different directions to get the neighboring cells. 

* Mapping: 

	For column decomposition, if columns can be divided evenly among the number of nodes, then each node will get the same amount of columns. If they are uneven, then some nodes will have one or two more than the others. 

	For 2D decomposition, if blocks can be divided evenly among the number of nodes, then each node will have same size of the block. If not, then some of the node will have a smaller block to compute. 


## Detail: 

My approach is to assume the grid of cells is square and the number of proccess has to be able to split into 2x2 or 3x3 grid in order to split the cells into blocks. 
I have a check to see if it is able to split, if not then it will terminate and ask for new setup.  

In `gol.f90`, I have just implement 4 nodes and 2x2 grid of nodes, to compute with the 20x20 grid. 


## Game of Life with 2D decomposition: 
* gol.f90: is the 2D  decomposition of Game of Life 

To compile and run the program
```
mpif90 -o gol gol.f90
```

boot up the compute node and go to the directory where you compile. 
```
mpirun -np <# of processes>(4) gol
```

## Game of Life with column decomposition: 
* gol_col.f90: is a column decomposition of the Game of Life

To compile and run the program
```
mpif90 -o gol_col gol_col.f90
```

boot up the compute node and go to the directory where you compile. 
```
mpirun -np <# of processes>(4) gol_col
```

## Game of Life with serial: 
* temp.f90: is a serial way of doing Game of Life 

To compile and run the program
```
gfortran -o temp temp.f90
```

boot up the compute node and go to the directory where you compile. 
```
./temp
```



Output of this porgram will use 20x20 grid and runs 80 iterations which will return the same output as initial location. 

If you want to change to different iteration, you can find the `do loop` to change the bounds of the loop to run 20, 40, 60, 80 iterations 

