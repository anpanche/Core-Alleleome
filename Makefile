.PHONY: install_blast install_mafft install_all

install_blast:
	@echo "Installing NCBI BLAST..."
	@tar zxvf ./resources/ncbi-blast-2.14.0+-x64-linux.tar.gz -C resources/
	@echo 'export PATH="$$PATH:resources/ncbi-blast-2.14.0+/bin"' >> ~/.bashrc
	@echo "NCBI BLAST installed successfully."

install_mafft:
	@echo "Installing MAFFT..."
	@tar xvf ./resources/mafft-7.505-with-extensions-src.tgz -C resources/
	@cd ./resources/mafft-7.505-with-extensions/core/ && make clean && make && sudo make install
	@echo "MAFFT installed successfully."

install_all: install_blast install_mafft
