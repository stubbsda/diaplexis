include Makefile.config

install: polyphyllon monophyllon
	mkdir -p $(INSTALL_DIR)/lib
	install -p src/polyphyllon/libdiaplexis.so $(INSTALL_DIR)/lib/libdiaplexis-polyphyllon.so
	install -p src/monophyllon/libdiaplexis.so $(INSTALL_DIR)/lib/libdiaplexis-monophyllon.so
	mkdir -p $(INSTALL_DIR)/include/diaplexis/polyphyllon
	mkdir -p $(INSTALL_DIR)/include/diaplexis/monophyllon
	install -p -m 444 src/polyphyllon/*.h $(INSTALL_DIR)/include/diaplexis/polyphyllon/
	install -p -m 444 src/monophyllon/*.h $(INSTALL_DIR)/include/diaplexis/monophyllon/

polyphyllon:
	cd src/polyphyllon; $(MAKE)

monophyllon:
	cd src/monophyllon; $(MAKE)

clean:
	cd src/polyphyllon; $(MAKE) clean
	cd src/monophyllon; $(MAKE) clean
	rm -f $(INSTALL_DIR)/lib/libdiaplexis-polyphyllon.so
	rm -f $(INSTALL_DIR)/lib/libdiaplexis-monophyllon.so
	rm -rf $(INSTALL_DIR)/include/diaplexis


