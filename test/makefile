INCLUDE:= -Iinclude -I../include

LIBS:= $(shell root-config --ldflags --libs --glibs)
LIBS += -L../lib -lG4Dictionary

SRC:= $(wildcard src/*.cc)
OBJ:= $(SRC:src/%.cc=bin/%.o)

MAINS:= RunTests

MAINSOBJ:= $(MAINS:%=bin/%.o)

CCFLAGS:=  -pthread -Wall -DExcludeInternalDependencies=1 -fPIC

G4OBJ:= ../lib/libG4Dictionary.a

ROOTCCFLAGS:=  -Wl,--no-as-needed $(shell root-config --cflags) -g

CC:= g++

.PHONY: clean

all: bin $(MAINS) $(DICTIONARY) $(G4OBJ)

bin:
	@echo Making directory for binaries
	@mkdir -p bin

makerun: all
	@echo Running tests...
	@./RunTests

test: $(G4OBJ)
	@echo Compiling main test object...
	@$(MAKE) -C test makerun

$(MAINS): %: $(filter-out $(MAINSOBJ), $(OBJ)) bin/%.o
	@echo Compiling $@
	@$(CC) $(CCFLAGS) $(INCLUDE) $(ROOTCCFLAGS) $^ $(LIBS) -o $@  


bin/%.o: src/%.cc include/%.hh $(G4OBJ)
	@echo Compiling $@
	@$(CC) $(CCFLAGS) $(INCLUDE) $(ROOTCCFLAGS) -c $< $(LIBS) -o $@


clean:
	@echo Cleaning up test directories...
	@rm -rf $(MAINS)
	@rm -rf bin/*
	@find . -type f -name "*.o" -exec rm -f {} \;
	@find . -type f -name "*.a" -exec rm -f {} \;
	@find . -type f -name "*.d" -exec rm -f {} \;
	@find . -type f -name "*.of" -exec rm -f {} \;
	@find . -type f -name "*~" -exec rm -f {} \;
	@find . -type f -name "\#*\#" -exec rm -f {} \;
	@find . -type f -name "*AutoDict*" -exec rm -f {} \;