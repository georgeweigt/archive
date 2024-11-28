.PHONY: default clean

DIRS =	feynman-hibbs \
	fitzpatrick \
	griffiths-schroeter \
	jordan \
	lancaster-blundell \
	miller \
	susskind-friedman

default:
	for FILE in $(DIRS) ; do make -C $$FILE ; done

clean:
	for FILE in $(DIRS) ; do make -C $$FILE clean ; done
