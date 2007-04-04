#===============================================================================
# File and Version Information:
#      $Id: GNUmakefile,v 1.1.1.1 2007-04-04 21:43:16 adye Exp $
#
# Description:
#      Makefile for the RooUnfold package
#
# Instructions:
#      o Review 'external configuration' section below
#        to match systems compilers setup
#
#      o Make sure the ROOTSYS environment variable is set and points
#        to your ROOT release.
#
#      o run 'make <target>'
#        Add NOROOFIT=1 to build example programs (RooUnfoldTest*)
#        in ROOT 4 or earlier.
#
# Build targets:
#      bin   - make shlib and example programs (default target)
#      lib   - make libRooUnfold.a
#      shlib - make libRooUnfold.so
#      clean - delete all intermediate and final build objects
#
# Author List:
#      Tim Adye <T.J.Adye@rl.ac.uk>
#
# Copyright Information:
#      Copyleft (C) 2007 Rutherford Appleton Laboratory
#
#===============================================================================

# --- External configuration ---------------------------------
include $(ROOTSYS)/test/Makefile.arch
CC       = $(CXX)
CCFLAGS  = -Wno-deprecated $(filter-out -Wall -Woverloaded-virtual,$(CXXFLAGS))
MFLAGS   = -MM -Wno-deprecated
SRCDIR   = $(CURDIR)/src/
WORKDIR  = $(CURDIR)/tmp/$(ARCH)/
LIBDIR   = $(CURDIR)/lib/$(ARCH)/
EXEDIR   = $(CURDIR)/
EXESRC   = $(CURDIR)/examples/
INCLUDES = -I$(SRCDIR)
# -------------------------------------------------------------

# Internal configuration
PACKAGE=RooUnfold
OBJDIR=$(WORKDIR)obj/
DEPDIR=$(WORKDIR)dep/
LIBS=-L$(LIBDIR)
ROOTLIBS += -lMinuit -lHtml
ifneq ($(EXPLLINKLIBS),)
EXPLLINKLIBS := $(LIBS) $(ROOTLIBS)
endif

ifneq ($(NOROOFIT),)
CPPFLAGS += -DNOROOFIT
else
LIBS += -lRooFit
endif

MAIN      = $(notdir $(wildcard $(EXESRC)*.cxx))
MAINEXE   = $(addprefix $(EXEDIR),$(patsubst %.cxx,%$(ExeSuf),$(MAIN)))
INCLUDES += -I$(ROOTSYS)/include
ROOTSYS  ?= ERROR_RootSysIsNotDefined
HLIST     = $(filter-out $(SRCDIR)$(PACKAGE)_LinkDef.h,$(wildcard $(SRCDIR)*.h)) $(SRCDIR)$(PACKAGE)_LinkDef.h
CINTFILE  = $(WORKDIR)$(PACKAGE)Cint.cxx
CINTOBJ   = $(OBJDIR)$(PACKAGE)Cint.o
LIBFILE   = $(LIBDIR)lib$(PACKAGE).a
SHLIBFILE = $(LIBDIR)lib$(PACKAGE).$(DllSuf)

default : shlib

# List of all object files to build
OLIST=$(addprefix $(OBJDIR),$(patsubst %.cxx,%.o,$(notdir $(wildcard $(SRCDIR)*.cxx))))

# List of all dependency file to make
DLIST=$(addprefix $(DEPDIR),$(patsubst %.cxx,%.d,$(notdir $(wildcard $(SRCDIR)*.cxx $(EXESRC)*.cxx))))

# Implicit rule making all dependency Makefiles included at the end of this makefile
$(DEPDIR)%.d : $(SRCDIR)%.cxx
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@rm -f $@
	@set -e; \
	 $(CC) $(MFLAGS) $(CPPFLAGS) $(INCLUDES) $< \
	 | sed 's,\($(notdir $*)\.o\) *:,$(OBJDIR)\1 $@ :,g' > $@; \
	 [ -s $@ ] || rm -f $@

$(DEPDIR)%.d : $(EXESRC)%.cxx
	@echo "Making $@"
	@mkdir -p $(DEPDIR)
	@rm -f $@
	@set -e; \
	 $(CC) $(MFLAGS) $(CPPFLAGS) $(INCLUDES) $< \
	 | sed 's,\($(notdir $*)\.o\) *:,$(OBJDIR)\1 $@ :,g' > $@; \
	 [ -s $@ ] || rm -f $@

# Implicit rule to compile all classes
$(OBJDIR)%.o : $(SRCDIR)%.cxx
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) $(CPPFLAGS) -c $< -o $(OBJDIR)$(notdir $@) $(INCLUDES)

# Implicit rule to compile main program
$(OBJDIR)%.o : $(EXESRC)%.cxx
	@echo "Compiling main program $<"
	@mkdir -p $(OBJDIR)
	@$(CC) $(CCFLAGS) $(CPPFLAGS) -c $< -o $(OBJDIR)$(notdir $@) $(INCLUDES)


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
	@mkdir -p $(LIBDIR)
	@rm -f $(SHLIBFILE)
	@$(LD) $(SOFLAGS) $(LDFLAGS) $(OLIST) $(CINTOBJ) $(OutPutOpt) $(SHLIBFILE) $(EXPLLINKLIBS)

$(MAINEXE) : $(EXEDIR)%$(ExeSuf) : $(OBJDIR)%.o $(SHLIBFILE)
	@echo "Making executable $@"
	@mkdir -p $(EXEDIR)
	@$(LD) $(LDFLAGS) $< $(OutPutOpt) $@ $(LIBS) -l$(PACKAGE) $(ROOTLIBS)

# Useful build targets
include: $(DLIST)
lib: $(LIBFILE)
shlib: $(SHLIBFILE)
bin: $(MAINEXE)

clean :
	rm -f $(DLIST)
	rm -f $(CINTFILE) $(basename $(CINTFILE)).h
	rm -f $(OLIST) $(CINTOBJ)
	rm -f $(LIBFILE)
	rm -f $(SHLIBFILE)
	rm -f $(addprefix $(OBJDIR),$(patsubst %.cxx,%.o,$(MAIN)))
	rm -f $(MAINEXE)

.PHONY : include shlib lib bin default clean

-include $(DLIST)
