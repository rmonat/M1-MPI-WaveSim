CC = gcc
CFLAGS = -Wall -Wextra
DFLAGS = -D DEBUG -g
RM = rm -rf
SRC = main.c
OBJ = main.o
OBJ_DEBUG = main.o.debug

%.o: %.c
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
