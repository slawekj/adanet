CC=g++
MEX=mex
CFLAGS=-Wall -fPIC
MKDIR=mkdir
TARDIR=decisionStump
OBJDIR=obj
SRCDIR=src

MEXSRC=$(SRCDIR)/adaboost.cpp
SOURCES=$(SRCDIR)/train_decision_stump.cpp
OBJECTS=$(OBJDIR)/train_decision_stump.o
MEXFILE=$(TARDIR)/adaboost.mexa64

all: $(SOURCES) $(MEXFILE)

$(MEXFILE): $(OBJECTS) $(TARDIR)
	$(MEX) $(MEXSRC) $(OBJECTS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJDIR):
	$(MKDIR) $(OBJDIR)

$(TARDIR):
	$(MKDIR) $(TARDIR)

clean:
	rm -rf $(OBJDIR) $(TARDIR)
