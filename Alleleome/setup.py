from setuptools import setup, find_packages
import subprocess



def install_blast():
    try:
        
        # Extract the downloaded archive
        subprocess.check_call(['tar', 'zxvf', './resources/ncbi-blast-2.14.0+-x64-linux.tar.gz', '-C', 'resources/'])

        # Add BLAST to the PATH
        blast_path = 'resources/ncbi-blast-2.14.0+/bin'
        subprocess.check_call(['echo', 'export PATH="$PATH:' + blast_path + '" >> ~/.bashrc'])
        
    except subprocess.CalledProcessError:
        print("Failed to install NCBI BLAST.")


def install_mafft():
    try:
        
        # Extract the downloaded archive
        subprocess.check_call(['tar', 'xvf', './resources/mafft-7.505-with-extensions-src.tgz', '-C', 'resources/'])
        
        mafft_core_dir='./resources/mafft-7.505-with-extensions/core/'
        # Make MAFFT executable
        subprocess.check_call(['make', 'clean'],cwd=mafft_core_dir)
        subprocess.check_call(['make'],cwd=mafft_core_dir)
        subprocess.check_call(['sudo', 'make', 'install'],cwd=mafft_core_dir)

    except subprocess.CalledProcessError:
        print("Failed to install MAFFT.")


# Call the BLAST installation function before installing your package
install_blast()

# Call the MAFFT installation function before installing your package
install_mafft()


setup(
    name="Alleleome",
    version="0.1",
    packages=find_packages(),
    install_requires=['pandas','numpy','Bio',
    ],  
    package_data={
        'Alleleome':['sample_data/Oenococcus_oeni/*']
    },
    entry_points={
        'console_scripts': ['Alleleome=Alleleome.cli:main']
    }
)
