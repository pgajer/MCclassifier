#############################################################################
#
# Makefile for building: trimAlign
#
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
FLAGS         = -O2 # -g -O0
CFLAGS        = $(FLAGS) -Wall
CXXFLAGS      = $(FLAGS) -Wall
INCPATH       = -Isrc
LINK          = g++
LIBS          = # -lm # -L../../lib -lkmerstats
AR            = ar cq
RANLIB        = ranlib -s
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = $(COPY_FILE)
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = $(COPY_FILE)
DEL_FILE      = rm -f
SYMLINK       = ln -sf
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Files

SRCDIR   = .
BUILDDIR = $(SRCDIR)/.build
BINDIR   = ../bin

create-build-dir := $(shell $(CHK_DIR_EXISTS) $(BUILDDIR) || $(MKDIR) $(BUILDDIR))
create-bin-dir := $(shell $(CHK_DIR_EXISTS) $(BINDIR) || $(MKDIR) $(BINDIR))

OBJECTS = $(BUILDDIR)/trimAlign.o    \
          $(BUILDDIR)/CUtilities.o   \
          $(BUILDDIR)/IOUtilities.o

####### Build rules

trimAlign: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) -o $(BINDIR)/trimAlign

$(BUILDDIR)/trimAlign.o: $(SRCDIR)/msa.h $(SRCDIR)/trimAlign.c $(SRCDIR)/IOUtilities.h \
	$(SRCDIR)/IOUtilities.c $(SRCDIR)/CUtilities.h $(SRCDIR)/CUtilities.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o $(BUILDDIR)/trimAlign.o $(SRCDIR)/trimAlign.c

$(BUILDDIR)/CUtilities.o: $(SRCDIR)/CUtilities.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o $(BUILDDIR)/CUtilities.o $(SRCDIR)/CUtilities.c

$(BUILDDIR)/IOUtilities.o: $(SRCDIR)/IOUtilities.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o $(BUILDDIR)/IOUtilities.o $(SRCDIR)/IOUtilities.c

clean:
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core
