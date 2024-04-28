#IDIR = include
IDIR=src
LDIR=lib
SDIR=src
ODIR=obj

.SUFFIXES:
.SUFFIXES: .cpp .o

CC=g++ --std=c++14
CFLAGS=-I$(IDIR) -Wuninitialized -Weffc++ --pedantic-errors

## files needed to compile stellar models
#  dependencies
_STARTYPES = \
	STARS/Star \
		STARS/Polytrope \
		STARS/ChandrasekharWD++ \
		STARS/MESA \
		STARS/SimpleWD \
		STARS/Isopycnic 
_STARDEPS = constants.h 

STARDEPS = $(patsubst %, $(IDIR)/%.h, $(_STARTYPES)) $(patsubst %, $(IDIR)/%, $(_STARDEPS))
STARSRC  = $(patsubst %, $(IDIR)/%.cpp, $(_STARTYPES))

## files needed to compile mode drivers
#  dependencies
_DRVTYPES = \
	MODES/NonradialModeDriver \
	MODES/CowlingModeDriver
_DRVDEPS = \
	constants.h\
	STARS/Star.h\
	MODES/ModeDriver.h 
DRVDEPS = $(patsubst %, $(IDIR)/%.h, $(_DRVTYPES)) $(patsubst %, $(IDIR)/%, $(_DRVDEPS))
DRVSRC  = $(patsubst %, $(SDIR)/%.cpp, $(_DRVTYPES))

## files needed to compile the mode object
#  dependencies
_MODEDEPS = constants.h\
  STARS/Star.h MODES/ModeDriver.h MODES/Mode.h
MODEDEPS = $(patsubst %, $(IDIR)/%, $(_MODEDEPS))
#  source
_MODESRC = MODES/Mode.cpp
MODESRC  = $(patsubst %, $(SDIR)/%, $(_MODESRC))


## files needed to compile main program
#  dependencies
_MAINDEPS = constants.h ThrainMain.h ThrainIO.h ThrainUnits.h
MAINDEPS = $(patsubst %, $(IDIR)/%, $(_MAINDEPS)) $(STARDEPS) $(MODEDEPS) $(DRVDEPS)
#  soure
_MAINSRC = ThrainMain.cpp ThrainIO.cpp ThrainUnits.cpp ThrainStellar.cpp ThrainMode.cpp
MAINSRC  = $(patsubst %, $(SDIR)/%, $(_MAINSRC))


## prepare object names
STAROBJ = $(patsubst %, $(ODIR)/%.o, $(_STARTYPES))
DRVOBJ  = $(patsubst %, $(ODIR)/%.o, $(_DRVTYPES))
MODEOBJ = $(patsubst %.cpp, $(ODIR)/%.o, $(_MODESRC))
MAINOBJ = $(patsubst %.cpp, $(ODIR)/%.o, $(_MAINSRC))


## Main rule for THRAIN program
thrain:  $(MAINOBJ) $(MODEOBJ) $(STAROBJ) $(DRVOBJ) # |library
	$(CC) -o $@ $^ $(CFLAGS) $(LDIR)/mylib.a -lm

## Rules for each subsection -- only update if their dependencies change
$(STAROBJ): $(ODIR)/%.o: $(SDIR)/%.cpp $(STARDEPS) |obj/STARS
	$(CC) -c -o $@ $< $(CFLAGS)

$(DRVOBJ): $(ODIR)/%.o: $(SDIR)/%.cpp $(DRVDEPS) |obj/MODES
	$(CC) -c -o $@ $< $(CFLAGS)

$(MODEOBJ): $(ODIR)/%.o: $(SDIR)/%.cpp $(MODEDEPS) |obj/MODES
	$(CC) -c -o $@ $< $(CFLAGS)

$(MAINOBJ): $(ODIR)/%.o: $(SDIR)/%.cpp  $(MAINDEPS) |obj/
	$(CC) -c -o $@ $< $(CFLAGS)


mode: $(MAINOBJ) $(MODEOBJ) $(STAROBJ) $(DRVOBJ)
	touch $(SDIR)/MODES/Mode.h
	rm $(ODIR)/MODES/Mode.o
	$(CC) -c -o $(ODIR)/MODES/Mode.o $(IDIR)/MODES/Mode.cpp $(CFLAGS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDIR)/mylib.a

# these will create the necessary directories
# see solution 4 here: 
#    https://www.cmcrossroads.com/article/making-directories-gnu-make
obj/:
	mkdir -p obj
obj/STARS:
	mkdir -p obj/STARS
obj/MODES:
	mkdir -p obj/MODES

tests: thrain tests/*.h
	cxxtestgen --error-printer -o tests/tests.cpp tests/*.h
#	this line makes cxxtest print to stderr so that stdout can be captured
#	sed 's/CxxTest::ErrorPrinter tmp;/CxxTest::ErrorPrinter tmp(std::cerr);/' \
#		tests/tests.cpp > changed.cpp && mv changed.cpp tests/tests.cpp
	$(CC) -o tests/tests.out \
		$(ODIR)/ThrainUnits.o $(ODIR)/ThrainMode.o $(ODIR)/ThrainIO.o $(ODIR)/ThrainStellar.o \
		$(MODEOBJ) $(STAROBJ) $(DRVOBJ) tests/tests.cpp $(CFLAGS) $(LDIR)/mylib.a

cppcheck:
	cppcheck lib/ --error-exitcode=1 --std=c++14
	cppcheck src/ --error-exitcode=1 --std=c++14

.PHONY: clean pull library

library:
	rm -f lib/*.o
	rm -f lib/*.a
	$(MAKE) -C lib --makefile=makelib library
	rm -f lib/*.o

## this command is used on my local machine to handle centralized versioning
pull:
	$(MAKE) -f pull

clean:
	rm -f $(ODIR)/*.o $(ODIR)/STARS/*.o $(ODIR)/MODES/*.o

