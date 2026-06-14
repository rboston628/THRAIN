#IDIR = include
IDIR=src
LDIR=lib
SDIR=src
ODIR=obj

.SUFFIXES:
.SUFFIXES: .cpp .o

CXX      ?= g++
CXXSTD   ?= -std=c++14
CXXFLAGS +=
CPPFLAGS += -I$(IDIR)
LDFLAGS  +=
LDLIBS   += -lm
CFLAGS=-I$(IDIR) -Wuninitialized -Weffc++ --pedantic-errors

# detect if using MSVC
IS_MSVC := $(filter cl,$(notdir $(CXX)))
ifeq ($(IS_MSVC),cl)
# MSVC flags
	CXXFLAGS += /std:c++14 /W4 /EHsc
else
# GCC/Clang flags
	CXXFLAGS += -std=c++14 -Wuninitialized -Weffc++ --pedantic-errors
endif

# Auto dependency files
CPPFLAGS += -MMD -MP

# If we're in a conda-like environment (pixi/conda/mamba), allow it to provide
# headers (e.g., cxxtest) without hard-coding paths.
ifneq ($(CONDA_PREFIX),)
	CFLAGS += -I$(CONDA_PREFIX)/include
endif

## files needed to compile stellar models
#  dependencies
_STARTYPES = \
	STARS/Star \
		STARS/Polytrope \
		STARS/ChandrasekharWD++ \
		STARS/MESA \
		STARS/SimpleWD
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

# collect dep files
DEPS = $(STAROBJ:.o=.d) $(DRVOBJ:.o=.d) $(MODEOBJ:.o=.d) $(MAINOBJ:.o=.d)

## Main rule for THRAIN program
thrain: $(MAINOBJ) $(MODEOBJ) $(STAROBJ) $(DRVOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDIR)/mylib.a $(LDLIBS)

$(ODIR)/%.o: $(SDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

mode: $(MAINOBJ) $(MODEOBJ) $(STAROBJ) $(DRVOBJ)
	$(MAKE) -B $(ODIR)/MODES/Mode.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDIR)/mylib.a

# include autogen dependencies if present
-include $(DEPS)

# -------------------------------------------------------------
# TESTS - mized doctest + cxxtest
# -------------------------------------------------------------

TDIR = tests

# Tools for CxxTest generation
CXXTESTGEN ?= cxxtestgen
SED        ?= sed

# Generated CxxTest runner
CXXTEST_CPP  := $(TDIR)/cxxtests.cpp
CXXTEST_OUT  := $(TDIR)/cxxtest.out
DOCTEST_OUT  := $(TDIR)/doctest.out

# Discover tests by extension
DOCTEST_SRCS := $(filter-out $(CXXTEST_CPP),$(wildcard $(TDIR)/*.cpp))
CXXTEST_HDRS := $(wildcard $(TDIR)/*.h)

# define test dependencies
# these are "mock" classes that are only useful for testing
_TESTSTARS = test_stars/Isopycnic test_stars/DummyStar
_TESTSTARDEPS = constants.h STARS/Star.h STARS/Star.cpp
_TESTMODES = test_modes/SineModeDriver test_modes/DummyMode
_TESTMODEDEPS = constants.h MODES/ModeDriver.h
TESTSTARDEPS = $(patsubst %, $(TDIR)/%.h, $(_TESTSTARS)) $(patsubst %, $(IDIR)/%, $(_TESTSTARDEPS))
TESTMODEDEPS = $(patsubst %, $(TDIR)/%.h, $(_TESTMODES)) $(patsubst %, $(IDIR)/%, $(_TESTMODEDEPS))
TESTSTARSRC = $(patsubst %, $(TDIR)/%.cpp, $(_TESTSTARS))
TESTMODESRC = $(patsubst %, $(TDIR)/%.cpp, $(_TESTMODES))
TESTSTAROBJ = $(patsubst %, $(ODIR)/%.o, $(_TESTSTARS))
TESTMODEOBJ = $(patsubst %, $(ODIR)/%.o, $(_TESTMODES))

TEST_CORE_OBJS := \
	$(ODIR)/ThrainUnits.o $(ODIR)/ThrainMode.o $(ODIR)/ThrainIO.o $(ODIR)/ThrainStellar.o \
	$(MODEOBJ) $(STAROBJ) $(DRVOBJ) \
	$(TESTMODEOBJ) $(TESTSTAROBJ)

# Build both during migration
tests: tests-doctest tests-cxxtest

# DOCTEST BUILD
tests-doctest: thrain $(DOCTEST_SRCS) $(TESTSTAROBJ) $(TESTMODEOBJ)
	$(CXX) $(LDFLAGS) -o $(DOCTEST_OUT) \
		$(TEST_CORE_OBJS) \
		$(DOCTEST_SRCS) -I$(TDIR) $(CPPFLAGS) $(CXXFLAGS) \
		$(LDIR)/mylib.a $(LDLIBS)

# CXXTEST BUILD
tests-cxxtest: thrain $(TESTSTAROBJ) $(TESTMODEOBJ)
	$(CXXTESTGEN) --error-printer -o $(CXXTEST_CPP) $(CXXTEST_HDRS)
#	this line makes cxxtest print to stderr so that stdout can be captured
	sed 's/CxxTest::ErrorPrinter tmp;/CxxTest::ErrorPrinter tmp(std::cerr);/' \
		$(CXXTEST_CPP) > changed.cpp && mv changed.cpp $(CXXTEST_CPP)
	$(CXX) $(LDFLAGS) -o $(CXXTEST_OUT) \
		$(TEST_CORE_OBJS) \
		$(CXXTEST_CPP) $(CPPFLAGS) $(CXXFLAGS) \
		$(LDIR)/mylib.a $(LDLIBS)

#TESTSRC := $(TDIR)/mode_test.cpp $(TDIR)/rootfind_test.cpp $(TDIR)/basic.cpp $(TDIR)/string_test.cpp $(TDIR)/logger_test.cpp $(TDIR)/fullcalc_test.cpp $(TDIR)/test_main.cpp $(TDIR)/matrix_test.cpp $(TDIR)/thrainunits_test.cpp
#TESTSRC := $(wildcard $(TDIR)/*.cpp)

# tests: thrain $(TESTSRC) $(TESTSTAROBJ) $(TESTMODEOBJ)
# 	$(CXX) $(LDFLAGS) -o $(TDIR)/tests.out $(TEST_CORE_OBJS)
# 		$(TESTSRC) -I$(TDIR) $(CPPFLAGS) $(CXXFLAGS) \
# 		$(LDIR)/mylib.a $(LDLIBS)

$(TESTSTAROBJ): $(ODIR)/%.o: $(TDIR)/%.cpp $(TESTSTARDEPS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(TDIR) -c -o $@ $<

$(TESTMODEOBJ): $(ODIR)/%.o: $(TDIR)/%.cpp $(TESTMODEDEPS)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(TDIR) -c -o $@ $<

cppcheck:
	mkdir -p cppcheck
	cppcheck --xml --xml-version=2 \
		--language=c++ \
		--std=c++14 \
		--enable=warning,performance,portability,information \
		--inconclusive \
		--inline-suppr \
		--force \
		--library=std.cfg \
		--suppress=missingIncludeSystem \
		--suppress=virtualCallInConstructor \
		--error-exitcode=1  \
		-I$(IDIR) -I$(SDIR) -I$(LDIR) \
		src/ lib/ \
		2> cppcheck/cppcheck_report.xml

.PHONY: clean cleantests pull library cppcheck

library:
	rm -f lib/*.o
	rm -f lib/*.a
	$(MAKE) -C lib --makefile=makelib library
	rm -f lib/*.o

## this command is used on my local machine to handle centralized versioning
pull:
	$(MAKE) -f pull

clean: cleantests
	rm -f $(ODIR)/*.o $(ODIR)/STARS/*.o $(ODIR)/MODES/*.o

cleantests:
	rm -f $(ODIR)/tests/*.o
	rm -f $(ODIR)/test_stars/*.o $(ODIR)/test_modes/*.o
