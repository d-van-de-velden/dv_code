import os
from setuptools import setup, find_packages


__version__ = '0.4.99'


def say_welcome(vers_info=None):
    print(
    '       ____\n'+
    '     /     \  __               ___    __\n'+
    '    /       ||  |             /__/   |  |\n'+
    '   /    ___/ |  | __  ______  __     |  | ______\n'+
    '   \   |___  |  |/ / |      ||  | ___|  ||      |\n'+
    '    \___   \ |    /  |  |_| ||  ||      ||  |_| |\n'+
    '       /   \ |    \  |   ___||  || |_|  ||   ___|\n'+
    '    __/    / |  |\ \ |  |___ |  ||      ||  |___\n'+
    '  /_______/  |__| \_\|______||__||______||______|\n'+
    '    ___                   __\n'+
    '   |   |                 |  |\n'+
    '   |   |                 |  |\n'+
    '   |   |                 |  |\n'+
    '   |   |              _  |  |___\n'+
    '   |   |         ____/ \ |      \ \n'+
    '   |   |_____   /      | |  |_| |\n'+
    '   |         \ |  |_|  | |      |\n'+
    '    \________/  \___/|_/ \______/  version: ' + __version__)
    return


def package_tree(pkgroot):
    """Get the submodule list."""
    # adapted from mne-python
    path = os.path.dirname(__file__)
    subdirs = [os.path.relpath(i[0], path).replace(os.path.sep, '.')
    for i in os.walk(os.path.join(path, pkgroot))
        if '__init__.py' in i[2]]
    
    say_welcome(vers_info=__version__)
    
    return sorted(subdirs)

############################################################################
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

    
setup(
    name='dv_code',
    version=__version__,
    author='Daniel van de Velden',
    author_email='daniel.vandevelden@yahoo.de',
    description='Python script repository for neuroimaging analysis',
    url='https://github.com/d-van-de-velden/dv_code/',
    license='BSD (3-clause)',
    packages=find_packages(),
    install_requires=[
        'pyocclient','numpy', 'dcm2bids','dcm2niiX', 'nilearn', 'nibabel'
        ],
    
)


