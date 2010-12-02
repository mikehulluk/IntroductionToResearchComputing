

#
all: build/document.pdf
	
clean:
	rm -f build/document.*

build/document.tex : document.tex
	cp document.tex build/document.tex

build/document.pdf: build/document.tex 
	make -C build

