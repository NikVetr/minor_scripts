#!/bin/bash
filename="$1"
basename="$(basename "$filename")"
dirname="$(dirname "$filename")"
cd "$dirname"
/Library/TeX/texbin/pdflatex --synctex=1 --interaction=batchmode "$basename"
/Library/TeX/texbin/bibtex "${basename%.tex}.aux"
/Library/TeX/texbin/pdflatex --synctex=1 --interaction=batchmode "$basename"
/Library/TeX/texbin/pdflatex --synctex=1 --interaction=batchmode "$basename"
