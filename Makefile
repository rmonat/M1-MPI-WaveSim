CC = gcc
CFLAGS = -Wall -Wextra
DFLAGS = -D DEBUG -g
RM = rm -rf
SRC = main.c automaton.c
DEPS = automaton.h
OBJ = main.o automaton.o
OBJ_DEBUG = main.o.debug automaton.o.debug

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
	$(RM) *.o *.aux *.log *.out
