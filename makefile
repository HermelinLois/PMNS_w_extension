# ---- STATIC PARAMETERS ----
CC = gcc
PyC = sage
OUTPUT_DIR = generated_pmns
TARGET_red = test_reduction
TARGET_conv = test_conversion
GENERATED_ELEMENTS = $(OUTPUT_DIR) $(TARGET_red) $(TARGET_conv)
GITIGNORE = .gitignore

# ---- PARAMETERS ----
NTEST ?= 10
NBITS ?= 128
K ?= 2
OPT ?= -O3 -funroll-loops
ETYPE ?= 0
METHOD ?= 0

# ---- PATH TO TARGETS ----
C_GEN = pmns_generator/orchestrator.py
SRC_GEN_red = $(OUTPUT_DIR)/code/reduction.c
SRC_GEN_conv = $(OUTPUT_DIR)/code/conversion.c

# ---- TAGS ----
START_TAG = <-- GENERATED-ELEMENTS-START -->
END_TAG = <-- GENERATED-ELEMENTS-END -->

# ---- RULES ----
all: update-ignore $(TARGET_red) $(TARGET_conv)

update-ignore:
	@touch $(GITIGNORE)
	@sed -i '/$(START_TAG)/,/$(END_TAG)/d' $(GITIGNORE)
	@echo "$(START_TAG)" >> $(GITIGNORE)
	@for file in $(GENERATED_ELEMENTS); do echo $$file >> $(GITIGNORE); done
	@echo "$(END_TAG)" >> $(GITIGNORE)

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


$(TARGET_conv): $(SRC_GEN_conv)
	@ $(CC) $(OPT) $(SRC_GEN_conv) -o $(TARGET_conv) -lgmp

$(TARGET_red): $(SRC_GEN_red)
	@ $(CC) $(OPT) $(SRC_GEN_red) -o $(TARGET_red)

$(SRC_GEN_conv) $(SRC_GEN_red) &: show-config $(C_GEN)
	@ $(PyC) $(C_GEN) -ntest $(NTEST) -nbits $(NBITS) -k $(K) -Etype $(ETYPE) -method $(METHOD) -name $(OUTPUT_DIR)

clean:
	@ sed -i '/$(START_TAG)/,/$(END_TAG)/d' $(GITIGNORE) 
	@ rm -rf $(GENERATED_ELEMENTS)

.PHONY: all clean