objects = globalvars.o randmaps.o
nmin?=12
nmax?=18
ns?=1024
nt?=4
DEFINES = -DnMin=$(nmin) -DnMax=$(nmax) -DFG_SAMPLE_N=$(ns) -DThreadsN=$(nt)
CC = mpiicpc
CFLAGS = -xHost -O3 -lippcp -lippcore  $(DEFINES)

FG : FG.cpp globalvars.cpp randmaps.cpp
	$(CC) $(CFLAGS) -c globalvars.cpp -o globalvars.o
	$(CC) $(CFLAGS) -c randmaps.cpp globalvars.o -o randmaps.o
	$(CC) $(CFLAGS) FG.cpp randmaps.o globalvars.o -o FG
	make clean


.PHONY:clean

clean:
	rm $(objects)