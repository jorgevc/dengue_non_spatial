include include_make

figure: Figure.tex
	pdflatex -shell-escape -synctex=1 -interaction=nonstopmode Figure.tex
	mv ./Figure.pdf ./$(id)_denisty.pdf
