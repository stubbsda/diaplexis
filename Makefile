include Makefile.def

install: diaplexis-ntl diaplexis-stl
	mkdir -p $(DIAPLEXIS)/lib
	cp NTL/libdiaplexis-ntl.so $(DIAPLEXIS)/lib/
	cp STL/libdiaplexis-stl.so $(DIAPLEXIS)/lib/
	mkdir -p $(DIAPLEXIS)/include/NTL
	mkdir -p $(DIAPLEXIS)/include/STL
	cp *.h $(DIAPLEXIS)/include/
	cp NTL/*.h $(DIAPLEXIS)/include/NTL/
	cp STL/*.h $(DIAPLEXIS)/include/STL/ 

manual:
	latex diaplexis.tex
	dvips -o diaplexis.ps diaplexis.dvi
	ps2pdf diaplexis.ps userguide.pdf
	mkdir -p $(DIAPLEXIS)/documentation
	mv userguide.pdf $(DIAPLEXIS)/documentation/

diaplexis-ntl: $(OBJECTS) 
	cd NTL; $(MAKE)  

diaplexis-stl: $(OBJECTS)
	cd STL; $(MAKE) 

clean:
	rm -f diaplexis.aux diaplexis.dvi diaplexis.log diaplexis.ps 
	rm -f *~
	rm -f $(DIAPLEXIS)/documentation/*
	rm -f $(DIAPLEXIS)/lib/*
	rm -rf $(DIAPLEXIS)/include/*
	cd NTL; make clean
	cd STL; make clean











