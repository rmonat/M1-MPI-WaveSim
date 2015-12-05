CC = mpicc
GCC = gcc
SCC = smpicc -pipe -march=native
CFLAGS = -Wall -Wextra -O3 -std=c99
DFLAGS = -D DEBUG -g
RM = rm -rf
SRC = main.c automaton.c args.c step0.c step123.c step4.c
DEPS = automaton.h args.h step0.h step123.h step4.h
OBJ = main.o automaton.o args.o step0.o step123.o step4.o
OBJ_DEBUG = main.o.debug automaton.o.debug args.o.debug step0.o.debug step123.o.debug step4.o.debug

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

%.o.debug: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(DFLAGS)


main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

main.debug: $(OBJ_DEBUG)
	$(CC) -o $@ $^ $(CFLAGS) $(DFLAGS)

main.sim: $(OBJ)
	$(SCC) -o $@ $^ $(CFLAGS)


young.o: young.c
	$(GCC) -c -o $@ $< $(CFLAGS)

young: young.o
	$(GCC) -o $@ $^ $(CFLAGS)


latex:
	pdflatex report.tex

clean:
	$(RM) *.o *.aux *.log *.out *.o.debug monat/ monat.tar.gz

dist: clean latex
	mkdir -p monat/
	cp $(SRC) monat/
	cp $(DEPS) monat/
	cp report.pdf monat/
	cp Makefile monat/
	cp LICENSE monat/
	cp -r results monat/
	tar -cvzf monat.tar.gz monat

demo: main
	@echo "\t\033[0;31msample 1 without walls\033[0m"
	mpirun -n 4 ./main -step 1 -iteration 100000 -dt 0.005 -grid 2 2 -i ./toolbox/sample_type1.in -lastdump res.dump
	./toolbox/viewer.py -i res.dump -f show -s 128 256
	@echo "\t\033[0;31msample 1 now with walls\033[0m"
	mpirun -n 4 ./main -step 2 -iteration 100000 -dt 0.005 -grid 2 2 -i ./toolbox/sample_type1.in -lastdump res_walls.dump
	./toolbox/viewer.py -i res_walls.dump -f show -s 128 256
