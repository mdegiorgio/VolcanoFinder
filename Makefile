#COPTS=-g -pg
COPTS = -O3 -Wall
#COPTS = -O0 -Wall

VolcanoFinder: VolcanoFinder.o freq.o factorials.o bfgs.o sort.o my_rand.o
	gcc -o VolcanoFinder VolcanoFinder.o $(COPTS) freq.o factorials.o bfgs.o sort.o my_rand.o -lm

VolcanoFinder.o: VolcanoFinder.c VolcanoFinder.h
	gcc -c VolcanoFinder.c $(COPTS)

freq.o: freq.c freq.h
	gcc -c freq.c $(COPTS)

factorials.o: factorials.c factorials.h
	gcc -c factorials.c $(COPTS)

bfgs.o: bfgs.c bfgs.h
	gcc -c bfgs.c $(COPTS)

sort.o: sort.c sort.h
	gcc -c sort.c $(COPTS)

my_rand.o: my_rand.c my_rand.h
	gcc -c my_rand.c $(COPTS)
