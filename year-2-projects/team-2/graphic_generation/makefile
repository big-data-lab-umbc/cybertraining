sample.pdf : sample.tex causal_graph_png.png
	pdflatex sample.tex
	
causal_graph_png.png : generate_graph.py graphic_supplemental.py node_locations.csv connections.csv
	generate_graph.py

display : sample.pdf
	(evince sample.pdf) &


clean : 
	rm -f *.pdf
	rm -f *.png
	rm -f *.aux
	rm -f *.log
