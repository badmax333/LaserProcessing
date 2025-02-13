all: clean main.exe
	mkdir output
	./main.exe
	rm -f main.exe




### Create the simulation program executable file.
main.exe: main.cpp
	clang main.cpp -lfftw3 -lm -o main.exe



### Delete the whole "output" directory and the simulation program executable file.
clean: 
	rm -f -R output
	rm -f main.exe

