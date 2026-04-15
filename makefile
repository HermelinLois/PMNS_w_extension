# ---- STATIC PARAMETERS ----
CC = gcc
PyC = sage
OUTPUT_DIR = generated_pmns
TARGET = test_reduction

# ---- PARAMETERS ----
NTEST ?= 10
NBITS ?= 128
K ?= 2
OPT ?= -O3 -funroll-loops
ETYPE ?= 0
METHOD ?= 0

# ---- PATH TO TARGETS ----
C_GEN = pmns_generator/orchestrator.py
SRC_GEN = $(OUTPUT_DIR)/code/reduction.c

# ---- RULES ----
all: $(TARGET)

show-config:
	@printf "=================================================================\n"
	@printf "|%18s%s%16s|\n" "" "PMNS GENERATION CONFIGURATION" ""
	@printf "|%15s%s%13s|\n" "" "[FORMAT] DECSRIPTION (NAME) : VALUE" ""
	@printf "=================================================================\n"
	@printf "| %-61s |\n" "PRIME BIT SIZE (NBITS) : $(NBITS)"
	@printf "| %-61s |\n" "EXTENSION DEGREE (K) : $(K)"
	@printf "| %-61s |\n" "EXTERNAL REDUCTION USED (ETYPE) : $(ETYPE)"
	@printf "| %-61s |\n" "REDUCTION METHOD USED (METHOD) : $(METHOD)"
	@printf "| %-61s |\n" "NUMBER OF TEST (NTEST) : $(NTEST)"
	@printf "|---------------------------------------------------------------|\n"
	@printf "| %-61s |\n" "COMPILATION OPTION (OPT) : $(OPT)"
	@printf "| %-61s |\n" "CODE GENERATED IN : $(OUTPUT_DIR)"
	@printf "=================================================================\n"

$(TARGET): $(SRC_GEN)
	@$(CC) $(OPT) $(SRC_GEN) -o $(TARGET)

$(SRC_GEN): show-config $(C_GEN)
	@$(PyC) $(C_GEN) -ntest $(NTEST) -nbits $(NBITS) -k $(K) -Etype $(ETYPE) -method $(METHOD) -name $(OUTPUT_DIR)

clean:
	@rm -rf $(OUTPUT_DIR)
	@rm -f $(TARGET)

.PHONY: all clean