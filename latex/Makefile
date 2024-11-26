.PHONY: default clean

# https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html

%.pdf: %.tex
	pdflatex $<
	if [ -e $*.toc ] ; then pdflatex $< ; pdflatex $< ; fi

default:
	for file in $$(basename -a -s .tex *.tex) ; do make $$file.pdf ; done

clean:
	rm -f *.aux *.log *.synctex.gz *.toc
