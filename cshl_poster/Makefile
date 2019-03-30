.PHONY : clean

all:
	pdflatex poster
	biber poster
	pdflatex poster
	
clean:
	-rm a0header.ps poster.aux  poster.log  poster.pdf
