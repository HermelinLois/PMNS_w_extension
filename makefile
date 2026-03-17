# ---- STATIC ELEMENTS FOR CONFIGURATION ----
CC = gcc
PyC = sage
TARGET = pmns

# ---- CONFIGURABLE PARAMETERS ----
NTEST ?= 1000
NBITS ?= 128
K ?= 2
OPT ?= -O3
ETYPE ?= 0
METHOD ?= 0

# ---- SCRIPTS PATH ----
C_GEN = pmns_generator/orchestrator.py
SRC_GEN = generated_code/reduction.c

# ---- RULES ----
all: $(TARGET)

$(TARGET): $(SRC_GEN)
	$(CC) $(OPT) $(SRC_GEN) -o $(TARGET)

$(SRC_GEN): $(C_GEN)
	$(PyC) $(C_GEN) -ntest $(NTEST) -nbits $(NBITS) -k $(K) -Etype $(ETYPE) -method $(METHOD)

clean:
	rm -f $(SRC_GEN) $(TARGET)

.PHONY: all clean $(SRC_GEN)
