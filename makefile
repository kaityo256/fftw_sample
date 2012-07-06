all: a.out

a.out: main.cc
	g++ -O3 $< -lfftw3

acf.dat: a.out
	./a.out > $@

graph: acf.png

acf.png: acf.dat
	gnuplot acf.plt

clean:
	rm -f a.out acf.dat acf.png
