INCLUDE:= -Iinclude -I../include -I../RLlib/include -I../libconfig-1.4.9/lib -I../lapacke -I../arpack

LIBS:= -lblas -larpack  -lm -llapack -llapacke -L../RLlib -lRLlib -L../libconfig-1.4.9/lib -Wl,-Bstatic -lconfig++ -Wl,-Bdynamic 

SRC:= $(wildcard src/*.cc)
OBJ:= $(SRC:src/%.cc=bin/%.o)

MAINS:= RunTests

MAINSOBJ:= $(MAINS:%=bin/%.o)

NON_TESTABLE_BIN:= ../bin/Compute.o

TESTABLE_BIN:= $(filter-out $(NON_TESTABLE_BIN), $(wildcard ../bin/*.o))


CCFLAGS:=  -pthread -Wall -fPIC -g -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP -pedantic -Wno-long-long -O3

CC:= g++

.PHONY: clean update_headers all makerun

.NOTPARALLEL: update_headers include/ExistingTests.cff include/AllHeaders.cff

all: bin update_headers $(MAINS) temp

bin:
	@echo Making directory for binaries...
	@mkdir -p bin

temp:
	@echo Making temporary file directory...
	@mkdir -p temp

makerun: all
	@echo Running tests...
	@./RunTests

update_headers: 
	@./testlist_create.sh
	@./allheaders_create.sh

test: 
	@echo Compiling main test object...
	@$(MAKE) -C test makerun

include/ExistingTests.cff: testlist_create.sh
	@echo Creating list of existing tests...
	@./testlist_create.sh
	@touch include/ExistingTests.cff

include/AllHeaders.cff: allheaders_create.sh
	@echo Creating list of all other headers...
	@./allheaders_create.sh
	@touch include/AllHeaders.cff

$(MAINS): %: $(filter-out $(MAINSOBJ), $(OBJ)) include/ExistingTests.cff include/AllHeaders.cff bin/%.o $(TESTABLE_BIN)
	@echo Compiling $@
	@$(CC) $(CCFLAGS) $(INCLUDE) $(filter-out include/ExistingTests.cff include/AllHeaders.cff, $^) $(LIBS) -o $@  


bin/%.o: src/%.cc include/%.hh include/ExistingTests.cff include/AllHeaders.cff
	@echo Compiling $@
	@$(CC) $(CCFLAGS) $(INCLUDE) -c $< -o $@


clean:
	@echo Cleaning up test directories...
	@rm -rf $(MAINS)
	@rm -rf bin/*
	@rm -rf include/ExistingTests*.cff
	@rm -rf include/AllHeaders*.cff
	@rm -rf temp/*
	@find . -type f -name "*.o" -exec rm -f {} \;
	@find . -type f -name "*.a" -exec rm -f {} \;
	@find . -type f -name "*.d" -exec rm -f {} \;
	@find . -type f -name "*.of" -exec rm -f {} \;
	@find . -type f -name "*~" -exec rm -f {} \;
	@find . -type f -name "\#*\#" -exec rm -f {} \;
