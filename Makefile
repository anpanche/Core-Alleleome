.PHONY: install_blast install_mafft install_all

install_blast:
	@echo "Installing NCBI BLAST..."
	@if [ ! -f resources/ncbi-blast-2.14.0+-x64-linux.tar.gz ] || [ `md5sum resources/ncbi-blast-2.14.0+-x64-linux.tar.gz | cut -d ' ' -f 1` != `wget -qO- https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz.md5 | cut -d ' ' -f 1` ]; then wget -O resources/ncbi-blast-2.14.0+-x64-linux.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz; fi
	@tar zxvf ./resources/ncbi-blast-2.14.0+-x64-linux.tar.gz -C resources/
	@echo 'export PATH="$$PATH:resources/ncbi-blast-2.14.0+/bin"' >> ~/.bashrc
	@echo "NCBI BLAST installed successfully."

install_mafft:
	@echo "Installing MAFFT..."
	@wget -O resources/mafft-7.505-with-extensions-src.tgz https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz
	@tar xvf ./resources/mafft-7.505-with-extensions-src.tgz -C resources/
	@cd ./resources/mafft-7.505-with-extensions/core/ && make clean && make && sudo make install
	@echo "MAFFT installed successfully."

install_all: install_blast install_mafft

.PHONY: install_blast install_mafft install_all
