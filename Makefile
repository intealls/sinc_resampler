CC = gcc
WARNINGS = -Wall -Wextra -pedantic -Wno-unused-function
CFLAGS = -g -O3 $(WARNINGS) -march=native -DRS_BLOCKSIZE=8 -DTYPE_IN=float -DTYPE_OUT=float -funroll-loops -ftree-vectorize
ODIR = bin
NAME = $(ODIR)/test_resamp

all: release

release: directories executable

debug: CFLAGS = -g3 -O0 $(WARNINGS)
debug: NAME := $(NAME)_dbg
debug: directories executable

directories:
	mkdir -p $(ODIR)

executable:
	$(CC) $(CFLAGS) sinc_resampler.c -lm -o $(NAME)

clean:
	rm -rf $(NAME) $(ODIR) data.csv *.s
