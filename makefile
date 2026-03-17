# ---- STATIC ELEMENTS FOR CONFIGURATION ----
CC = gcc
PyC = sage
TARGET = pmns

# ---- CONFIGURABLE PARAMETERS ----
NTEST ?= 1000
BITS ?= 128
EXTENSION ?= 2
OPT ?= -O3
TYPE ?= 0
METHOD ?= 1

# ---- SCRIPTS PATH ----
C_GEN = pmns_generator/orchestrator.py
SRC_GEN = generated_code/reduction.c

# ---- RULES ----
all: $(TARGET)

$(TARGET): $(SRC_GEN)
	$(CC) $(OPT) $(SRC_GEN) -o $(TARGET)

$(SRC_GEN): $(C_GEN)
	$(PyC) $(C_GEN) -ntest $(NTEST) -nbits $(BITS) -k $(EXTENSION) -Etype $(TYPE) -method $(METHOD)

clean:
	rm -f $(SRC_GEN) $(TARGET)

.PHONY: all clean
