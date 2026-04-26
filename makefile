SAGE ?= sage
CC ?= gcc
CFLAGS ?= -O3 -Wall -Wextra -std=gnu11
GMP_LIBS ?= -lgmp
GENERATOR ?= pmns_generator/orchestrator.py
OUTPUT_DIR = pmns_exec

# Generator parameters
NTEST ?= 100
NBITS ?= 128
K ?= 2
ETYPE ?= 0

GENERATED_FILES = \
	$(OUTPUT_DIR)/config_pmns \
	$(OUTPUT_DIR)/params/pmns_params.h \
	$(OUTPUT_DIR)/params/reductions_params.h \
	$(OUTPUT_DIR)/params/conversions_params.h \
	$(OUTPUT_DIR)/tests/conversions_values.h \
	$(OUTPUT_DIR)/tests/reductions_values.h

TEST_REDUCTION_SRC = $(OUTPUT_DIR)/tests/test_reduction.c
REDUCTIONS_SRC = $(OUTPUT_DIR)/codes/reductions.c
TEST_REDUCTION_BIN = test_reduction
TEST_CONVERSION_SRC = $(OUTPUT_DIR)/tests/test_conversion.c
CONVERSIONS_SRC = $(OUTPUT_DIR)/codes/conversions.c
TEST_CONVERSION_BIN = test_conversion

.PHONY: all help show-config check-sage generate test-reduction test-conversion clean

all: generate

help:
	@echo "Targets:"
	@echo "  make generate         -> run pmns generator"
	@echo "  make test-reduction   -> generate then compile reduction test"
	@echo "  make test-conversion  -> generate then compile conversion test"
	@echo "  make clean            -> remove generated files"
	@echo "  make show-config      -> print current generation config"
	@echo ""
	@echo "Parameters (override with VAR=value):"
	@echo "  NTEST=$(NTEST) NBITS=$(NBITS) K=$(K) ETYPE=$(ETYPE) WITH_VERIF=$(WITH_VERIF)"
	@echo ""
	@echo "Example:"
	@echo "  make generate NTEST=500 NBITS=256 K=3 ETYPE=1 WITH_VERIF=False"

show-config:
	@printf "===============================================================\n"
	@printf "|%20s%s%18s|\n" "" "PMNS GENERATOR CONFIG" ""
	@printf "===============================================================\n"
	@printf "| %-61s |\n" "SAGE CMD : $(SAGE)"
	@printf "| %-61s |\n" "GENERATOR : $(GENERATOR)"
	@printf "| %-61s |\n" "OUTPUT DIR : $(OUTPUT_DIR)"
	@printf "|-------------------------------------------------------------|\n"
	@printf "| %-61s |\n" "NTEST : $(NTEST)"
	@printf "| %-61s |\n" "NBITS : $(NBITS)"
	@printf "| %-61s |\n" "K : $(K)"
	@printf "| %-61s |\n" "ETYPE : $(ETYPE)"
	@printf "| %-61s |\n" "WITH_VERIF : $(WITH_VERIF)"
	@printf "===============================================================\n"

check-sage:
	@command -v $(SAGE) >/dev/null 2>&1 || { echo "Error: '$(SAGE)' not found in PATH"; exit 127; }

generate: check-sage show-config
	@mkdir -p $(OUTPUT_DIR)
	@$(SAGE) $(GENERATOR) \
		-ntest $(NTEST) \
		-nbits $(NBITS) \
		-k $(K) \
		-Etype $(ETYPE)
	@echo "Generation complete."

$(TEST_REDUCTION_BIN): generate $(TEST_REDUCTION_SRC) $(REDUCTIONS_SRC)
	@$(CC) $(CFLAGS) \
		-DN_TEST=$(NTEST) \
		-I$(OUTPUT_DIR)/tests \
		-I$(OUTPUT_DIR)/codes \
		$(TEST_REDUCTION_SRC) $(REDUCTIONS_SRC) \
		$(GMP_LIBS) \
		-o $(TEST_REDUCTION_BIN)
	@echo "Built $(TEST_REDUCTION_BIN)."

test-reduction: $(TEST_REDUCTION_BIN)

$(TEST_CONVERSION_BIN): generate $(TEST_CONVERSION_SRC) $(CONVERSIONS_SRC) $(REDUCTIONS_SRC)
	@$(CC) $(CFLAGS) \
		-DN_TEST=$(NTEST) \
		-I$(OUTPUT_DIR)/tests \
		-I$(OUTPUT_DIR)/codes \
		$(TEST_CONVERSION_SRC) $(CONVERSIONS_SRC) $(REDUCTIONS_SRC) \
		$(GMP_LIBS) \
		-o $(TEST_CONVERSION_BIN)
	@echo "Built $(TEST_CONVERSION_BIN)."

test-conversion: $(TEST_CONVERSION_BIN)

clean:
	@rm -f $(GENERATED_FILES) $(TEST_REDUCTION_BIN) $(TEST_CONVERSION_BIN)
	@echo "Generated files removed (no directory deleted)."









# # ---- STATIC PARAMETERS ----
# CC = gcc
# PyI = sage
# OUTPUT_DIR = pmns_exec
# TEST_DIR = $(OUTPUT_DIR) / tests

# GENERATED_ELEMENTS = $(OUTPUT_DIR) $(TARGET_red) $(TARGET_conv)
# GENERATED_FILES = $(OUTPUT_DIR)/config_pmns $(OUTPUT_DIR)/params $(OUTPUT_DIR)/tests/*.h 

# # ---- PARAMETERS ----
# NTESTS ?= 100
# NBITS ?= 128
# K ?= 2
# OPT ?= -O3 -funroll-loops
# ETYPE ?= 0

# ifeq LOAD
#     LOAD_FLAG = --load
# else
#     LOAD_FLAG =
# endif

# # ---- PATH TO TARGETS ----
# C_GEN = pmns_generator/orchestrator.py
# TARGET_TEST_red = $(TEST_DIR) / test_reduction.c
# TARGET_TEST_conv = $(TEST_DIR) / test_conversion.c

# # ---- RULES ----
# all: $(TARGET_red) $(TARGET_conv)

# show-config:
# 	@printf "=================================================================\n"
# 	@printf "|%18s%s%16s|\n" "" "PMNS GENERATION CONFIGURATION" ""
# 	@printf "|%15s%s%13s|\n" "" "[FORMAT] DECSRIPTION (NAME) : VALUE" ""
# 	@printf "=================================================================\n"
# 	@printf "| %-61s |\n" "PRIME BIT SIZE (NBITS) : $(NBITS)"
# 	@printf "| %-61s |\n" "EXTENSION DEGREE (K) : $(K)"
# 	@printf "| %-61s |\n" "EXTERNAL REDUCTION USED (ETYPE) : $(ETYPE)"
# 	@printf "| %-61s |\n" "NUMBER OF TEST (NTESTS) : $(NTESTS)"
# 	@printf "|---------------------------------------------------------------|\n"
# 	@printf "| %-61s |\n" "COMPILATION OPTION (OPT) : $(OPT)"
# 	@printf "| %-61s |\n" "EXTRA FLAGS (LOAD) : $(LOAD)"
# 	@printf "=================================================================\n"


# $(TARGET_conv): $(SRC_GEN_conv)
# 	@ $(CC) $(OPT) $(SRC_GEN_conv) -o $(TARGET_conv) -lgmp

# $(TARGET_red): $(SRC_GEN_red)
# 	@ $(CC) $(OPT) $(SRC_GEN_red) -o $(TARGET_red)

# $(SRC_GEN_conv) $(SRC_GEN_red) &: show-config $(C_GEN)
# 	@ $(PyI) $(C_GEN) -ntests $(NTESTS) -nbits $(NBITS) -k $(K) -Etype $(ETYPE) $(LOAD_FLAG)

# clean:
# 	@ rm -rf $(GENERATED_ELEMENTS)

# .PHONY: all clean