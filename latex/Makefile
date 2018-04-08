all: build/thesis.pdf build/thesis-compressed.pdf


TeXOptions = -xelatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build \
			 --synctex=1
                                                                                
build/thesis.pdf: FORCE | build
	latexmk $(TeXOptions) thesis.tex
	
FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build

build/thesis-compressed.pdf: build/thesis.pdf
	gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dNOPAUSE -dQUIET -dBATCH -sOutputFile=build/thesis-compressed.pdf build/thesis.pdf

