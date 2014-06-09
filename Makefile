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

diaplexis-ntl: $(OBJECTS) 
	cd NTL; $(MAKE)  

diaplexis-stl: $(OBJECTS)
	cd STL; $(MAKE) 

clean:
	rm -f diaplexis.aux diaplexis.dvi diaplexis.log diaplexis.ps 
	rm -f *~
	rm -f $(DIAPLEXIS)/lib/libdiaplexis*
	rm -f $(DIAPLEXIS)/include/diaplexis*
	rm -rf $(DIAPLEXIS)/include/NTL
	rm -rf $(DIAPLEXIS)/include/STL
	cd NTL; make clean
	cd STL; make clean











