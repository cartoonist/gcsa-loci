##  Directories:
# External libraries directory.
EXTDIR       := ext
# Source files directory.
SRCDIR       := src
# Binary files directory.
BINDIR       := ./bin
# The location where the binary files should be installed.
PREFIX       ?= ~/.local

##  Sources:
SOURCES  := $(wildcard ${SRCDIR}/*.cc)
SOURCES  += $(wildcard ${SRCDIR}/*.h)

##  Header-only libraries:
HEADERONLY_LIBS = spdlog kseq++

# Specifying phony targets.
.PHONY: all tags install install-debug clean distclean

## Functions:
define echotitle
	@echo
	@tput setaf 2
	@echo "$1"
	@tput sgr0
	@echo
endef

define uninstall_headers
	rm -rfv $(addprefix ${SRCDIR}/, ${HEADERONLY_LIBS})
endef

all: init
	$(call echotitle,"Building sources...")
	@make -C ${SRCDIR}

init-header:
	$(call echotitle,"Copying header-only libraries...")

init: init-header $(addprefix ${SRCDIR}/, ${HEADERONLY_LIBS})

${SRCDIR}/spdlog:
	@cp -rv ${EXTDIR}/spdlog/include/spdlog $@

${SRCDIR}/kseq++:
	@cp -rv ${EXTDIR}/kseq++/src $@

tags:
	$(call echotitle,"Updating ctags...")
	ctags ${SOURCES}

install: all
	$(call echotitle,"Installing binaries...")
	@install -v ${BINDIR}/gcsa_loci ${PREFIX}/bin

clean:
	@make -C ${SRCDIR} $@

distclean:
	$(call uninstall_headers)
	@make -C ${SRCDIR} $@
