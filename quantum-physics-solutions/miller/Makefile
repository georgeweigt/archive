.PHONY: default clean

%.pdf: %.tex preamble
	pdflatex $<

%.html: %.txt ../text1 ../text2
	cat ../text1 $< ../text2 > $@

default:
	for file in $$(basename -a -s .tex *.tex) ; do make $$file.pdf ; done
	for file in $$(basename -a -s .txt *.txt) ; do make $$file.html ; done

clean:
	rm -f *.aux *.log *.synctex.gz
