build: 
	cd NTL; $(MAKE)
	cd STL; $(MAKE)

NTL: build
	mkdir -p $(BIBLIOTHEK)/lib
	install -p NTL/libdiaplexis.so $(BIBLIOTHEK)/lib/
	mkdir -p $(BIBLIOTHEK)/include/diaplexis
	install -p -m 444 NTL/*.h $(BIBLIOTHEK)/include/diaplexis/

STL: build
	mkdir -p $(BIBLIOTHEK)/lib
	install -p STL/libdiaplexis.so $(BIBLIOTHEK)/lib/
	mkdir -p $(BIBLIOTHEK)/include/diaplexis
	install -p -m 444 STL/*.h $(BIBLIOTHEK)/include/diaplexis/

clean:
	rm -f *~
	rm -f $(BIBLIOTHEK)/lib/libdiaplexis.so
	rm -rf $(BIBLIOTHEK)/include/diaplexis
	cd NTL; make clean
	cd STL; make clean











