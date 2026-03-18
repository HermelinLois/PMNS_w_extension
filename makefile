# ---- STATIC ELEMENTS FOR CONFIGURATION ----
CC = gcc
PyC = sage
OUTPUT_DIR = generated_code
TARGET = pmns

# ---- CONFIGURABLE PARAMETERS ----
NTEST ?= 100
NBITS ?= 128
K ?= 2
OPT ?= -O3 -funroll-loops
ETYPE ?= 0
METHOD ?= 0

# ---- SCRIPTS PATH ----
C_GEN = pmns_generator/orchestrator.py
SRC_GEN = $(OUTPUT_DIR)/reduction.c

# ---- RULES ----

all: $(TARGET)

$(TARGET): $(SRC_GEN)
	$(CC) $(OPT) $(SRC_GEN) -o $(TARGET)

$(SRC_GEN): $(C_GEN)
	$(PyC) $(C_GEN) -ntest $(NTEST) -nbits $(NBITS) -k $(K) -Etype $(ETYPE) -method $(METHOD) -name $(OUTPUT_DIR)

clean:
	rm -f $(SRC_GEN) $(TARGET)
	rm -r $(OUTPUT_DIR)

.PHONY: all clean $(SRC_GEN)
