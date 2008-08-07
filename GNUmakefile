#===============================================================================
# File and Version Information:
#      $Id: GNUmakefile,v 1.7 2008-08-07 01:13:21 adye Exp $
#
# Description:
#      Makefile for the RooUnfold package
#
# Instructions:
#      o Make sure the ROOTSYS environment variable is set and points
#        to your ROOT release, and that $ROOTSYS/bin is in your PATH.
#
#      o run 'make <target>'
#        - Default target makes shared library (libRooUnfold.so) for use
#          in ROOT.
#        - Add NOROOFIT=1 to build test programs (RooUnfoldTest*)
#          without RooFit (this is default if RooFit is not available).
#        - Add SHARED=1 to link test executables with shared library
#          (libRooUnfold.so). Otherwise links with static library
#          (libRooUnfold.a).
#        - Add ROOTBUILD=debug for debug version.
#
# Build targets:
#      shlib   - make libRooUnfold.so (default target)
#      include - make dependency files (*.d)
#      lib     - make libRooUnfold.a
#      bin     - make lib and example programs
#      commands- show commands to make each type of target
#      clean   - delete all intermediate and final build objects
#
# Author List:
#      Tim Adye <T.J.Adye@rl.ac.uk>
#
# Copyright Information:
#      Copyleft (C) 2007 Rutherford Appleton Laboratory
#
#===============================================================================

# === ROOT setup ===============================================================
-include $(ROOTSYS)/test/Makefile.arch
ifeq ($(ROOTCONFIG),)
ROOTCONFIG    = $(ROOTSYS)/bin/root-config
endif
ifeq ($(ARCH),)
# === This section is just in case ROOT's test/Makefile.arch didn't work =======
out := $(shell echo "$(ROOTSYS)/test/Makefile.arch not found - trying a basic Linux config" >&2)
ARCH          =   $(shell $(ROOTCONFIG) --arch)
ifeq ($(ARCH),)
out := $(shell echo "$(ROOTCONFIG) did not work - assume standard locations below $(ROOTSYS)" >&2)
ARCH          =   $(shell uname | tr '[A-Z]' '[a-z]')
ROOTLIBS      = -L$(ROOTSYS)/lib -lCore -lCint -lHist -lGraf -lGpad -lPostscript -lMatrix -ldl
ROOTINCLUDES  = -I$(ROOTSYS)/include
CXXFLAGS      =   $(ROOTINCLUDES)
NOROOFIT      = 1
else
ROOTLIBS      =   $(shell $(ROOTCONFIG) --libs)
ROOTINCLUDES  = -I$(shell $(ROOTCONFIG) --incdir)
CXXFLAGS      =   $(shell $(ROOTCONFIG) --cflags)
endif
CXX           = g++
CXXFLAGS     += -Wall -fPIC
LD            = g++
LDFLAGS       =
SOFLAGS       = -shared
ObjSuf        = o
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"
ifneq ($(findstring debug,$(ROOTBUILD)),)
CXXFLAGS     += -g
LDFLAGS      += -g
endif
else
ROOTINCLUDES  = -I$(shell $(ROOTCONFIG) --incdir)
endif

# === RooUnfold directories and options ========================================

CC            = $(CXX)
ifeq ($(CC),g++)
CCFLAGS       = $(filter-out -Woverloaded-virtual,$(CXXFLAGS)) -Wno-deprecated -Wno-parentheses -Wno-sign-compare
MFLAGS        = -MM -Wno-deprecated
else
CCFLAGS       = $(CXXFLAGS)
endif
SRCDIR        = $(CURDIR)/src/
WORKDIR       = $(CURDIR)/tmp/$(ARCH)/
LIBDIR        = $(CURDIR)/
SHLIBDIR      = $(CURDIR)/
EXEDIR        = $(CURDIR)/
EXESRC        = $(CURDIR)/examples/
INCLUDES      = -I$(SRCDIR)

# === Internal configuration ===================================================

PACKAGE       = RooUnfold
OBJDIR        = $(WORKDIR)obj/
DEPDIR        = $(WORKDIR)dep/

ifeq ($(NOROOFIT),)
ifneq ($(shell $(ROOTCONFIG) --has-roofit),yes)
out := $(shell echo "This version of ROOT does not support RooFit. We will build the test programs without it." >&2)
NOROOFIT      = 1
endif
endif

ifneq ($(NOROOFIT),)
CPPFLAGS     += -DNOROOFIT
else
ifneq ($(wildcard $(ROOTSYS)/lib/libRooFitCore.$(DllSuf)),)
ROOFITLIBS   += -lRooFit -lRooFitCore -lThread -lMinuit -lHtml
else
ROOFITLIBS   += -lRooFit -lMinuit -lHtml
endif
endif

MAIN          = $(notdir $(wildcard $(EXESRC)*.cxx))
MAINEXE       = $(addprefix $(EXEDIR),$(patsubst %.cxx,%$(ExeSuf),$(MAIN)))
ROOTSYS      ?= ERROR_RootSysIsNotDefined
HLIST         = $(filter-out $(SRCDIR)$(PACKAGE)_LinkDef.h,$(wildcard $(SRCDIR)*.h)) $(SRCDIR)$(PACKAGE)_LinkDef.h
CINTFILE      = $(WORKDIR)$(PACKAGE)Cint.cxx
CINTOBJ       = $(OBJDIR)$(PACKAGE)Cint.o
LIBFILE       = $(LIBDIR)lib$(PACKAGE).a
SHLIBFILE     = $(SHLIBDIR)lib$(PACKAGE).$(DllSuf)

ifneq ($(SHARED),)
LIBS          = -L$(SHLIBDIR)
LINKLIB       = $(SHLIBFILE)
LINKLIBOPT    = -l$(PACKAGE)
else
LIBS          = -L$(LIBDIR)
LINKLIB       = $(LIBFILE)
LINKLIBOPT    = -Wl,-static -l$(PACKAGE) -Wl,-Bdynamic
endif

# List of all object files to build
OLIST         = $(addprefix $(OBJDIR),$(patsubst %.cxx,%.o,$(notdir $(wildcard $(SRCDIR)*.cxx))))

ifeq ($(MFLAGS),)

# Can't make dependency files, so make every compilation dependent on all headers.
HDEP          = $(HLIST)

else

# List of all dependency file to make
DLIST         = $(addprefix $(DEPDIR),$(patsubst %.cxx,%.d,$(notdir $(wildcard $(SRCDIR)*.cxx $(EXESRC)*.cxx))))

ifeq ($(NOROOFIT),)
# List of programs that use RooFit. Should only be those in $(EXESRC)
ROOFITCLIENTS = $(patsubst $(DEPDIR)%.d,$(OBJDIR)%.o,$(shell fgrep -l '/RooGlobalFunc.h ' $(DLIST) 2>/dev/null))
endif
endif

# === Implicit rules ===========================================================

# Implicit rule making all dependency Makefiles included at the end of this makefile
$(DEPDIR)%.d : $(SRCDIR)%.cxx
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@rm -f $@
	@set -e; \
	 $(CC) $(MFLAGS) $(CPPFLAGS) $(INCLUDES) $(ROOTINCLUDES) $< \
	 | sed 's,\($(notdir $*)\.o\) *:,$(OBJDIR)\1 $@ :,g' > $@; \
	 [ -s $@ ] || rm -f $@

$(DEPDIR)%.d : $(EXESRC)%.cxx
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@rm -f $@
	@set -e; \
	 $(CC) $(MFLAGS) $(CPPFLAGS) $(INCLUDES) $(ROOTINCLUDES) $< \
	 | sed 's,\($(notdir $*)\.o\) *:,$(OBJDIR)\1 $@ :,g' > $@; \
	 [ -s $@ ] || rm -f $@

# Implicit rule to compile all classes
$(OBJDIR)%.o : $(SRCDIR)%.cxx $(HDEP)
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) $(CPPFLAGS) -c $< -o $(OBJDIR)$(notdir $@) $(INCLUDES)

# Implicit rule to compile main program
$(OBJDIR)%.o : $(EXESRC)%.cxx $(HDEP)
	@echo "Compiling main program $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) $(CPPFLAGS) -c $< -o $(OBJDIR)$(notdir $@) $(INCLUDES)


# === Explicit rules ===========================================================

default : shlib

# Rule to make ROOTCINT output file
$(CINTOBJ) : $(HLIST)
	@mkdir -p $(WORKDIR)
	@mkdir -p $(OBJDIR)
	@echo "Running rootcint for $(SRCDIR)$(PACKAGE)_LinkDef.h"
	@cd $(SRC) ; $(ROOTSYS)/bin/rootcint -f $(CINTFILE) -c -p $(INCLUDES) $(HLIST)
	@echo "Compiling $(CINTFILE)"
	@$(CC) $(CCFLAGS) $(CPPFLAGS) -c $(CINTFILE) -o $(CINTOBJ) $(INCLUDES)

# Rule to combine objects into a library
$(LIBFILE) : $(OLIST) $(CINTOBJ)
	@echo "Making $(LIBFILE)"
	@mkdir -p $(LIBDIR)
	@rm -f $(LIBFILE)
	@ar q $(LIBFILE) $(OLIST) $(CINTOBJ)
	@ranlib $(LIBFILE)

# Rule to combine objects into a shared library
$(SHLIBFILE) : $(OLIST) $(CINTOBJ)
	@echo "Making $(SHLIBFILE)"
	@mkdir -p $(SHLIBDIR)
	@rm -f $(SHLIBFILE)
	@$(LD) $(SOFLAGS) $(LDFLAGS) $(OLIST) $(CINTOBJ) $(OutPutOpt)$(SHLIBFILE) $(ROOTLIBS)

$(MAINEXE) : $(EXEDIR)%$(ExeSuf) : $(OBJDIR)%.o $(LINKLIB)
	@echo "Making executable $@"
	@mkdir -p $(EXEDIR)
	@$(LD) $(LDFLAGS) $< $(OutPutOpt)$@ $(LIBS) $(LINKLIBOPT) $(ROOTLIBS) $(if $(findstring $<,$(ROOFITCLIENTS)),$(ROOFITLIBS))

# Useful build targets
include: $(DLIST)
lib: $(LIBFILE)
shlib: $(SHLIBFILE)
bin: $(MAINEXE)

commands :
	@echo "Make $(DEPDIR)%.d:	$(CC) $(MFLAGS) $(CPPFLAGS) $(INCLUDES) $(ROOTINCLUDES) $(SRCDIR)%.cxx | sed 's,\(%\.o\) *:,$(OBJDIR)\1 $(DEPDIR)%.d :,g' > $(DEPDIR)%.d"
	@echo
	@echo "Compile $(SRCDIR)%.cxx:	$(CC) $(CCFLAGS) $(CPPFLAGS) -c $(SRCDIR)%.cxx -o $(OBJDIR)%.o $(INCLUDES)"
	@echo
	@echo "Make $(SHLIBFILE):	$(LD) $(SOFLAGS) $(LDFLAGS) *.o $(OutPutOpt)$(SHLIBFILE) $(ROOTLIBS)"
	@echo
	@echo "Make executable $(EXEDIR)RooUnfoldTest$(ExeSuf):	$(LD) $(LDFLAGS) $(OBJDIR)RooUnfoldTest.o $(OutPutOpt)$(EXEDIR)RooUnfoldTest$(ExeSuf) $(LIBS) $(LINKLIBOPT) $(ROOTLIBS) $(ROOFITLIBS)"

clean :
	rm -f $(DLIST)
	rm -f $(CINTFILE) $(basename $(CINTFILE)).h
	rm -f $(OLIST) $(CINTOBJ)
	rm -f $(LIBFILE)
	rm -f $(SHLIBFILE)
	rm -f $(addprefix $(OBJDIR),$(patsubst %.cxx,%.o,$(MAIN)))
	rm -f $(MAINEXE)

.PHONY : include shlib lib bin default clean

ifneq ($(DLIST),)
-include $(DLIST)
endif
