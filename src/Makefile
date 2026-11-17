# Build the FEAL-8 key-recovery tools.
#
#   make            CPU tools (attack stages 1-3, test-data generators, analysis programs)
#   make gpu        CUDA tools (needs nvcc; set ARCH, e.g. ARCH=sm_120 for RTX 5090)
#   make check      quick self-checks (test vector, key verification, CBC-MAC of the data)
#   make clean

CC      = gcc
CFLAGS  ?= -O3 -march=native
OMP     ?= -fopenmp
NVCC    ?= nvcc
ARCH    ?= sm_86
NVFLAGS ?= -O3 -arch=$(ARCH)
PYTHON  ?= python3
DATA    ?= ../data

CPU_MAIN = attackC2 stage3 mkct mkdt mkct_seed
CPU_ANAL = analysis/attackC analysis/corr analysis/null analysis/null2 \
           analysis/refstat analysis/bias analysis/mlc
GPU_ALL  = attackC2_gpu bs_gpu analysis/attackC_gpu

all: $(CPU_MAIN) $(CPU_ANAL)
gpu: $(GPU_ALL)

attackC2: attackC2.c feal8.h approx8.h
	$(CC) $(CFLAGS) $(OMP) -o $@ $< -lm
stage3: stage3.c feal8.h
	$(CC) $(CFLAGS) $(OMP) -o $@ $< -lm
mkct mkdt mkct_seed: %: %.c feal8.h
	$(CC) $(CFLAGS) -o $@ $< -lm

analysis/%: analysis/%.c feal8.h approx8.h
	$(CC) $(CFLAGS) $(OMP) -I. -o $@ $< -lm

attackC2_gpu: attackC2_gpu.cu feal8.h approx8.h
	$(NVCC) $(NVFLAGS) -o $@ $<
bs_gpu: bs_gpu.cu feal8.h
	$(NVCC) $(NVFLAGS) -o $@ $<
analysis/attackC_gpu: analysis/attackC_gpu.cu feal8.h
	$(NVCC) $(NVFLAGS) -I. -o $@ $<

check:
	$(PYTHON) verify_key.py $$(cat $(DATA)/key.txt) $(DATA)/plaintext.bin $(DATA)/ciphertext.bin
	$(PYTHON) check_mac.py $(DATA)

clean:
	rm -f $(CPU_MAIN) $(CPU_ANAL) $(GPU_ALL)

.PHONY: all gpu check clean
