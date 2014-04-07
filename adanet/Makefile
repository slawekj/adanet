CC=g++
MEX=mex
CFLAGS=-Wall -fPIC
MKDIR=mkdir -p
TARDIR=decisionStump
OBJDIR=obj
SRCDIR=src

MEXSRC=$(SRCDIR)/adaboost.cpp
SOURCES=$(SRCDIR)/train_decision_stump.cpp
OBJECTS=$(OBJDIR)/train_decision_stump.o
MEXFILE=$(TARDIR)/adaboost.mexa64

all: $(SOURCES) $(MEXFILE)

$(MEXFILE): $(OBJECTS)
	$(MKDIR) $(TARDIR)
	$(MEX) $(MEXSRC) $(OBJECTS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(MKDIR) $(OBJDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR) $(TARDIR)
