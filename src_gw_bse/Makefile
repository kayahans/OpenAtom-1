include config.mk

C_FILES=$(wildcard *.C)
CI_FILES=$(wildcard *.ci)
MODULES=$(CI_FILES:.ci=.decl.h)
OBJS=$(addprefix $(BUILDDIR)/, $(C_FILES:.C=.o))

default: all
all: $(BUILDDIR) $(BUILDDIR)/libmodulemain.a pzheev

.PHONY : $(DIRS)

$(BUILDDIR):
	mkdir $(BUILDDIR)

$(BUILDDIR)/libmodulemain.a: make_ci $(OBJS) $(DIRS)
	$(CHARMC) $(LFLAGS) $(OBJS) $(BUILDDIR)/*.o -o $@
#	mv charmrun $(BUILDDIR)

LIBS=build/libmodulemain.a
pzheev: diagonalizer/pzheev.cpp $(LIBS)
	$(CHARMC) -c++ mpicxx $(CFLAGS) -c ./diagonalizer/pzheev.cpp -o ./diagonalizer/pzheev.o
	$(CHARMC) -ld++ mpicxx -mpi -nomain-module -o gw_bse ./diagonalizer/pzheev.o -L./ -L./build/  $(LFLAGS) -module CkMulticast  -module CkLoop -module main

make_ci: $(MODULES)
	for dir in $(DIRS); do make -C $${dir} make_ci; done

$(DIRS):
	make -C $@

%.decl.h: %.ci
	$(CHARMC) $<

$(BUILDDIR)/%.o: %.C %.h $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.C $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@


clean:
	rm -rf $(BUILDDIR) *.decl.h *.def.h
	for dir in $(DIRS); do make -C $${dir} clean; done
