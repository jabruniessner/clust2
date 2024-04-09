# The clust2 program
# Copyright (C) 2024 Jakob Niessner.
# Contact: jabruniessner@gmail.com
#
#
# This program was written at the Heidelberg Institute for theoretical studies,
# Schloß-Wolfsbrunnenweg 35, 69118 Heidelberg, Germany
# Under the supervision of Prof. Dr. Rebecca C. Wade
# Contact: rebecca.wade@h-its.org
#
# Permission to use, copy, modify, and distribute this software and its
# documentation with or without modifications and for any purpose and
# without fee is hereby granted, provided that any copyright notices
# appear in all copies and that both those copyright notices and this
# permission notice appear in supporting documentation, and that the
# names of the contributors or copyright holders not be used in
# advertising or publicity pertaining to distribution of the software
# without specific prior permission.
#
# THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
# OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
# OR PERFORMANCE OF THIS SOFTWARE.


CC=gcc -g

all: clust cut

cut: cut.o pdb_structure.o math_func.o trajectories.o cluster.o shared_functions.o
	$(CC) cut.o pdb_structure.o math_func.o trajectories.o cluster.o shared_functions.o -o cut -lm

clust: main.o pdb_structure.o math_func.o trajectories.o cluster.o shared_functions.o
	$(CC) main.o pdb_structure.o math_func.o trajectories.o cluster.o shared_functions.o -o clust -lm

cut.o: cut.c
	$(CC) -c cut.c -o cut.o

main.o: main.c
	$(CC) -c main.c -o main.o

shared_functions.o: shared_functions.c
	$(CC) -c shared_functions.c -o shared_functions.o

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

