# Defines compiler and linker settings for CHAP Makefiles

CC = gcc
CFLAGS = -Wall -Wextra
CFLAGS += -g -ggdb
LFLAGS = -lm

# suppress new warning re: "variable 'foo' set but not used" in gcc 4.6 and above
ifeq "$(CC)" "gcc"
  MAJOR := $(shell gcc -dumpversion | cut -d '.' -f1)
  MINOR := $(shell gcc -dumpversion | cut -d '.' -f2)
  NEWGCC := $(shell [ $(MAJOR) -gt 4 ] || \
                   ([ $(MAJOR) -eq 4 ] && [ $(MINOR) -ge 6 ]); echo $$?)
else
  NEWGCC = 1    # false, in shell-speak
endif
ifeq "$(NEWGCC)" "0"    # true
  CFLAGS += -Wno-unused-but-set-variable
endif
