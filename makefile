# ---- STATIC PARAMETERS ----
CC = gcc
PyC = sage
OUTPUT_DIR = generated_code
TARGET = pmns

# ---- PARAMETERS ----
NTEST ?= 100				# number of elements randomly generated
NBITS ?= 128				# bit size of prime used
K ?= 2						# extension degree used
OPT ?= -O3 -funroll-loops	# compilation option
ETYPE ?= 0					# type of external polynomial construction used
METHOD ?= 0					# reduction method used

# ---- PATH TO TARGETS ----
C_GEN = pmns_generator/orchestrator.py
SRC_GEN = $(OUTPUT_DIR)/reduction.c

# ---- RULES ----
all: $(TARGET)

$(TARGET): $(SRC_GEN)
	@$(CC) $(OPT) $(SRC_GEN) -o $(TARGET)

$(SRC_GEN): $(C_GEN)
	@$(PyC) $(C_GEN) -ntest $(NTEST) -nbits $(NBITS) -k $(K) -Etype $(ETYPE) -method $(METHOD) -name $(OUTPUT_DIR)

clean:
	@rm -rf $(OUTPUT_DIR)
	@rm -f $(TARGET)

.PHONY: all clean