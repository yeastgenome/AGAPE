include ../Flags.mk

ARCH ?= $(shell uname -m)
MAIN = ..
DEST= ../$(MAIN)/bin
KLIB = $(MAIN)/ucsc/lib/$(ARCH)
KINC = $(MAIN)/ucsc/inc

CLIB = $(KLIB)/jkweb.a -lm
CFLAGS += -I. -I$(KINC)

homologs: 
	$(CC) $(CFLAGS) \
	util.c util_i.c regions.c util_exons.c repeats.c readChain.c \
	$(CLIB) \
	-o $(DEST)/homologs 
