.PHONY: default clean

DIRS =	feynman-hibbs \
	fitzpatrick \
	griffiths-schroeter \
	jordan \
	lancaster-blundell \
	miller \
	susskind-friedman

default:
	for DIR in $(DIRS) ; do make -C $$DIR ; done

clean:
	for DIR in $(DIRS) ; do make -C $$DIR clean ; done
