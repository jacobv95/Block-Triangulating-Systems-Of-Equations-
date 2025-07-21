Example of block triangulating a system of equations.

This represents the system of equations as a bipartite graph. 

It then uses Hopcraft-Kart to find a maximum matching of the bipartite graph.

Finally, Tarjans algorithm is used to find strongly connected components.

The strongly connected components are then returned as blocks. 

The block include equation and variable indexes.

The blocks can be solve in succession. This reduces the number of variables, that has to be solved for at any given time.
