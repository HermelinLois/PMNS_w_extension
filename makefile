# ---- STATIC PARAMETERS ----
CC = gcc
PyI = sage
OUTPUT_DIR = pmns_exec

C_GEN = pmns_generator/orchestrator.py
TEST_DIR = $(OUTPUT_DIR)/tests
CODE_DIR = $(OUTPUT_DIR)/codes

GENERATED_FILES = $(OUTPUT_DIR)/config_pmns $(OUTPUT_DIR)/params $(OUTPUT_DIR)/tests/*.h 

SRC_red = $(CODE_DIR)/reductions.c
SRC_conv = $(CODE_DIR)/conversions.c

TEST_SRC_red = $(TEST_DIR)/test_reductions.c
TEST_SRC_conv = $(TEST_DIR)/test_conversions.c

TARGET_TEST_red = test_reductions
TARGET_TEST_conv = test_conversions

# ---- PARAMETERS ----
NTESTS ?= 100
NBITS ?= 128
K ?= 2
OPT ?= -O3 -funroll-loops
ETYPE ?= 0


GMP_INC ?= $(CONDA_PREFIX)/include
GMP_LIB ?= $(CONDA_PREFIX)/lib

CFLAGS = -Wall -Wextra -std=gnu11 $(OPT) -I$(GMP_INC)
LDFLAGS = -L$(GMP_LIB) -Wl,-rpath,$(GMP_LIB) -lgmp

ifdef LOAD
    LOAD_FLAG = --load
	SHOW_LOAD = TRUE
else
    LOAD_FLAG =
	SHOW_LOAD = FALSE
endif

# ---- RULES ----
all: $(TARGET_TEST_red) $(TARGET_TEST_conv)

show-config:
	@printf "=================================================================\n"
	@printf "|%18s%s%16s|\n" "" "PMNS GENERATION CONFIGURATION" ""
	@printf "|%15s%s%13s|\n" "" "[FORMAT] DECSRIPTION (NAME) : VALUE" ""
	@printf "=================================================================\n"
	@printf "| %-61s |\n" "PRIME BIT SIZE (NBITS) : $(NBITS)"
	@printf "| %-61s |\n" "EXTENSION DEGREE (K) : $(K)"
	@printf "| %-61s |\n" "EXTERNAL REDUCTION USED (ETYPE) : $(ETYPE)"
	@printf "| %-61s |\n" "NUMBER OF TESTS (NTESTS) : $(NTESTS)"
	@printf "|---------------------------------------------------------------|\n"
	@printf "| %-61s |\n" "COMPILATION OPTION (OPT) : $(OPT)"
	@printf "| %-61s |\n" "LOAD PRECOMPUTE PMNS (LOAD) : $(SHOW_LOAD)"
	@printf "=================================================================\n"

$(TARGET_TEST_red): $(TEST_SRC_red) $(SRC_red)
	$(CC) $(CFLAGS) -I$(TEST_DIR) -I$(CODE_DIR) $(TEST_SRC_red) $(SRC_red) $(LDFLAGS) -o $@

$(TARGET_TEST_conv): $(TEST_SRC_conv) $(SRC_conv) $(SRC_red)
	$(CC) $(CFLAGS) -I$(TEST_DIR) -I$(CODE_DIR) $(TEST_SRC_conv) $(SRC_conv) $(SRC_red) $(LDFLAGS) -o $@


generate: show-config $(C_GEN)
	@ $(PyI) $(C_GEN) -ntests $(NTESTS) -nbits $(NBITS) -k $(K) -Etype $(ETYPE) $(LOAD_FLAG)

$(SRC_red) $(SRC_conv): generate

clean:
	@ rm -rf $(GENERATED_FILES)
	@ rm -f $(TARGET_TEST_red) $(TARGET_TEST_conv)

.PHONY: generate