CC=gcc -g



clust: main.o pdb_structure.o math_func.o trajectories.o cluster.o
	$(CC) main.o pdb_structure.o math_func.o trajectories.o cluster.o -o clust -lm

main.o: main.c
	$(CC) -c main.c -o main.o

pdb_structure.o: pdb_structure.c
	$(CC) -c pdb_structure.c -o pdb_structure.o

math_func.o: math_func.c
	$(CC) -c math_func.c -o math_func.o

trajectories.o: trajectories.c
	$(CC) -c trajectories.c -o trajectories.o

cluster.o: cluster.c
	$(CC) -c cluster.c -o cluster.o -lm

clean:
	rm *.o
	rm clust

