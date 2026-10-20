# Compiler settings
FC = gfortran
FFLAGS = -O2 -Wall

# Directories
SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod
BIN_DIR = bin

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.f90)
OBJS = $(SRCS:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)

# Executable name
EXE = $(BIN_DIR)/sim

# Make sure the order is correct for module dependencies
ORDERED_SRCS = $(SRC_DIR)/parameters.f90 \
               $(SRC_DIR)/particles.f90 \
               $(SRC_DIR)/forces.f90 \
               $(SRC_DIR)/integrator.f90 \
               $(SRC_DIR)/io.f90 \
               $(SRC_DIR)/main.f90

ORDERED_OBJS = $(ORDERED_SRCS:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)

.PHONY: all clean directories

all: directories $(EXE)

# Create necessary directories
directories:
	@mkdir $(OBJ_DIR)
	@mkdir $(MOD_DIR)
	@mkdir $(BIN_DIR)

# Link the executable
$(EXE): $(ORDERED_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -J$(MOD_DIR) -c $< -o $@

# Clean build files
clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(BIN_DIR)
