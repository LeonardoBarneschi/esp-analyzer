.PHONY: all clean 

all: getmap.so

getmap.so: getmap.f90
		f2py -c -m getmap getmap.f90 
			
alternative:
		f2py --f90exec=/usr/bin/gfortran -c -m getmap getmap.f90

clean:
		rm -rf getmap.so
