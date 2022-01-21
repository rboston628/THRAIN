#IDIR = include
IDIR=src
LDIR=lib
SDIR=src
ODIR=obj

.SUFFIXES:
.SUFFIXES: .cpp .o

CC=g++ -std=c++11
CFLAGS=-I$(IDIR) -Wuninitialized -Weffc++

#VPATH = src include lib

## files needed to compile stellar models
#  dependencies
_STARDEPS = constants.h\
		STARS/Star.h \
		STARS/Polytrope.h \
		STARS/ChandrasekharWD++.h \
		STARS/MESA.h\
		STARS/SimpleWD.h
STARDEPS = $(patsubst %, $(IDIR)/%, $(_STARDEPS))
#  source
_STARSRC = STARS/Star.cpp \
		STARS/Polytrope.cpp\
		STARS/ChandrasekharWD++.cpp\
		STARS/MESA.cpp\
		STARS/SimpleWD.cpp
STARSRC  = $(patsubst %, $(SDIR)/%, $(_STARSRC))

## files needed to compile mode drivers
#  dependencies
_DRVDEPS = constants.h\
  STARS/Star.h\
  MODES/Mode.h MODES/Mode.cpp\
  MODES/ModeDriver.h MODES/NonradialModeDriver.h MODES/CowlingModeDriver.h
DRVDEPS = $(patsubst %, $(IDIR)/%, $(_DRVDEPS))
#  source
_DRVSRC  = MODES/NonradialModeDriver.cpp MODES/CowlingModeDriver.cpp
DRVSRC   = $(patsubst %, $(SDIR)/%, $(_DRVSRC))

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
_MAINDEPS = constants.h ThrainMain.h ThrainIO.h\
  STARS/Star.h STARS/Polytrope.h STARS/ChandrasekharWD++.h\
  MODES/Mode.h\
  MODES/ModeDriver.h MODES/NonradialModeDriver.h
MAINDEPS = $(patsubst %, $(IDIR)/%, $(_MAINDEPS))
#  soure
_MAINSRC = ThrainMain.cpp ThrainIO.cpp ThrainStellar.cpp ThrainMode.cpp ThrainUnits.cpp
MAINSRC  = $(patsubst %, $(SDIR)/%, $(_MAINSRC))


## prepare object names
STAROBJ = $(patsubst %.cpp,$(ODIR)/%.o, $(_STARSRC))
DRVOBJ  = $(patsubst %.cpp,$(ODIR)/%.o, $(_DRVSRC))
MODEOBJ = $(patsubst %.cpp,$(ODIR)/%.o, $(_MODESRC))
MAINOBJ = $(patsubst %.cpp,$(ODIR)/%.o, $(_MAINSRC))


## Main rule for THRAIN program
thrain:  $(MAINOBJ) $(MODEOBJ) $(STAROBJ) $(DRVOBJ) |library
	$(CC) -o $@ $^ $(CFLAGS) $(LDIR)/mylib.a

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

library: 
	$(CC) -c -o lib/Splinor.o lib/Splinor.cpp
	$(CC) -c -o lib/chandra.o lib/chandra.cpp
	$(CC) -c -o lib/stellar.o lib/stellar.cpp
	ar rcs lib/mylib.a lib/Splinor.o lib/chandra.o lib/stellar.o

obj/:
	mkdir -p obj
obj/STARS:
	mkdir -p obj/STARS
obj/MODES:
	mkdir -p obj/MODES


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(ODIR)/STARS/*.o $(ODIR)/MODES/*.o

