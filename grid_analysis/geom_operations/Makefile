.PHONY: all clean 

all: getpot.so getcoulene.so getpot_parallel.so getcoulene_parallel.so

getpot.so: getpot.f90
	f2py -c -m getpot getpot.f90 

getcoulene.so: getcoulene.f90
	f2py -c -m getcoulene getcoulene.f90

getpot_parallel.so: getpot_parallel.f90
	f2py --opt='-O3 -funroll-loops' -lgomp --f90flags='-fopenmp' -c -m getpot_parallel getpot_parallel.f90

getcoulene_parallel.so: getcoulene_parallel.f90
	f2py --opt='-O3 -funroll-loops' -lgomp --f90flags='-fopenmp' -c -m getcoulene_parallel getcoulene_parallel.f90

clean:
		rm -rf getpot.so
		rm -rf getcoulene.so
		rm -rf getpot_parallel.so
		rm -rf getcoulene_parallel.so
