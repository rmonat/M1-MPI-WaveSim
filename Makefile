CC = gcc
CFLAGS = -Wall -Wextra -O3
DFLAGS = -D DEBUG -g
RM = rm -rf
SRC = main.c automaton.c args.c step0.c
DEPS = automaton.h args.h step0.h
OBJ = main.o automaton.o args.o step0.o
OBJ_DEBUG = main.o.debug automaton.o.debug args.o.debug step0.o.debug

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

%.o.debug: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(DFLAGS)

main.debug: $(OBJ_DEBUG)
	$(CC) -o $@ $^ $(CFLAGS) $(DFLAGS)

latex:
	pdflatex report.tex

clean:
	$(RM) *.o *.aux *.log *.out *.o.debug
