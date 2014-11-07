install: diaplexis-ntl diaplexis-stl
	mkdir -p $(DIAPLEXIS)/lib
	install -p NTL/libdiaplexis-ntl.so $(DIAPLEXIS)/lib/
	install -p STL/libdiaplexis-stl.so $(DIAPLEXIS)/lib/
	mkdir -p $(DIAPLEXIS)/include/NTL
	mkdir -p $(DIAPLEXIS)/include/STL
	install -p -m 444 *.h $(DIAPLEXIS)/include/
	install -p -m 444 NTL/*.h $(DIAPLEXIS)/include/NTL/
	install -p -m 444 STL/*.h $(DIAPLEXIS)/include/STL/ 

diaplexis-ntl:
	cd NTL; $(MAKE)  

diaplexis-stl:
	cd STL; $(MAKE) 

clean:
	rm -f *~
	rm -f $(DIAPLEXIS)/lib/libdiaplexis*
	rm -f $(DIAPLEXIS)/include/diaplexis*
	rm -rf $(DIAPLEXIS)/include/NTL
	rm -rf $(DIAPLEXIS)/include/STL
	cd NTL; make clean
	cd STL; make clean











