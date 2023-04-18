# Usage:
# make -f print-variable.mk -f <target makefile> _print-<target variable>

_print-%:
	@echo "---PRINTING $*---"
	@echo ${$*}
