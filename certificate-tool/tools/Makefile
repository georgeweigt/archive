.PHONY: default

# https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html

%: %.c
	$(CC) -Wall -O0 -o $@ $<

default:
	for FILE in $$(basename -a -s .c *.c) ; do make $$FILE ; done
