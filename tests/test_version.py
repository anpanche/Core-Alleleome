import re
import Alleleome

def test_version():
    """Check if version in setup.py and __init__.py match"""
    with open('setup.py') as f:
        setup_py_content = f.read()
    setup_version = re.search(r'version="(\d+\.\d+\.\d+)"', setup_py_content).group(1)
    package_version = Alleleome.__version__
    assert (
        setup_version == package_version
    ), "Version mismatch between setup.py and Alleleome/__init__.py"
