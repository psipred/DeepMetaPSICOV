######################################################################
#	Makefile for MetaPSICOV                                      #
######################################################################

#### CC		Command to execute C compiler
#### CFLAGS	Flags to pass to C compiler.

INCDIR		= .
CC		= cc
CPP             = g++
CFLAGS		= -O
LIBS		= -lm
GD              = -lgd

all:		plotmap

plotmap:	plotmap.cpp
		$(CPP) $(CFLAGS) plotmap.cpp $(LIBS) $(GD) -o plotmap
