
# By default (if you just run "make"), AGAPE is configured to keep its scripts
# in the main "package directory" (where the distribution
# was unpacked), while compiled binaries files are located
# in the "bin" subdirectories, respectively.  If you want to
# install it elsewhere, edit the following lines to specify the desired
# locations, and then run "make install".
#
INSTALL = install
TARGETS = build_utils

build : clean bindir $(TARGETS)
	@echo '=== Build complete ==='

bindir :
	@echo ---$@---
	mkdir -p bin

build_lastz :
	@echo ---$@---
	cd lastz/src; make lastz
	cp lastz/src/lastz ../bin/lastz

build_homologs :
	@echo ---$@---
	cd homologs; make

build_utils :
	@echo ---$@---
	cd utils; make

clean :
	@echo ---$@---
	cd lastz/src; make clean
	rm -rf bin
