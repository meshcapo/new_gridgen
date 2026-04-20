#
# Gridexec
#
# @version 0.1
.PHONY: code clean test print-%
.DEFAULT_GOAL := code

# Rule to print variables from Makefile
print-% : ; @echo $* = $($*)

# Directories
IDIR = include
ODIR = obj
SDIR = source
EDIR = bin
MKDIR_P := mkdir -p
TEST_DIR = tests

# Executable
EXEC = $(EDIR)/gridexec
TESTEXEC = $(EDIR)/tests_runner

# Compiler and flags
CC = clang
CFLAGS = -I$(IDIR) -g -O0
CRITERION_CFLAGS := $(shell pkg-config --cflags criterion 2>/dev/null || echo)

# Linker flags and libraries
LDFLAGS = -g
LDLIBS  = -lgsl -lgslcblas -lm
CRITERION_LIBS := $(shell pkg-config --libs criterion 2>/dev/null || echo -lcriterion -lpthread)

# Header files
_DEPS = globvar.h types.h errors.h aux_functions.h file_operations.h geometry.h jacobian_ders.h newton.h grid_solvers.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# Source files for tests
TEST_SRC = $(wildcard $(TEST_DIR)/*.c)

# Object files
_OBJ = errors.o aux_functions.o file_operations.o geometry.o jacobian_ders.o newton.o grid_solvers.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
TEST_OBJ = $(patsubst $(TEST_DIR)/%.c, $(TEST_DIR)/%.o, $(TEST_SRC))

# Objects for linking tests: exclude main.o so Criterion's main is used
OBJ_NOMAIN = $(filter-out $(ODIR)/main.o, $(OBJ))

$(EDIR):
	$(MKDIR_P) $(EDIR)
$(ODIR):
	$(MKDIR_P) $(ODIR)
$(ODIR)/%.o: $(SDIR)/%.c $(DEPS) | $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS)
$(TEST_DIR)/%.o: $(TEST_DIR)/%.c $(DEPS) | $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS) $(CRITERION_CFLAGS)
$(EXEC): $(OBJ) | $(EDIR)
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)
$(TESTEXEC): $(OBJ_NOMAIN) $(TEST_OBJ) | $(EDIR)
	$(CC) -o $@ $^ $(LDFLAGS) $(CRITERION_LIBS) $(LDLIBS)

# 'make code'
code: $(EXEC)

# 'make test'
test: $(TESTEXEC)

# 'make clean'
clean:
	rm -f $(EXEC) $(ODIR)/*.o $(TEST_DIR)/*.o
	rm -rf $(EDIR) $(ODIR)

# end
