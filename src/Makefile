
CC = gcc
CFLAGS = -O3 -o
LDFLAGS = -lm

objects = psipred psipass2 solvpred cov21stats deepmetapsicov_makepredmap alnstats psicov

all: $(objects)

psipred:	sspred_avpred.c ssdefs.h sspred_net.h
	$(CC) $(CFLAGS) $@ $< $(LDFLAGS)

psipass2:	sspred_hmulti.c ssdefs.h sspred_net2.h
	$(CC) $(CFLAGS) $@ $< $(LDFLAGS)

solvpred:	solvpred_avpred.c solvation_net.h
	$(CC) $(CFLAGS) $@ $< $(LDFLAGS)

cov21stats:	cov21stats.c
	$(CC) $(CFLAGS) $@ $< $(LDFLAGS)

deepmetapsicov_makepredmap:	deepmetapsicov_makepredmap.c
	$(CC) $(CFLAGS) $@ $< $(LDFLAGS)

alnstats:	alnstats.c
	$(CC) $(CFLAGS) $@ $< $(LDFLAGS)

psicov:	psicov25.c
	$(CC) $(CFLAGS) $@ $< $(LDFLAGS) -fopenmp

install:	all
	/bin/mkdir -p ../bin
	/bin/mv $(objects) ../bin

clean:
	/bin/rm -rf ../bin
	/bin/rm -f $(objects)
