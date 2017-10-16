from setuptools import setup

setup(
    name='simfiles',
    version='0.1',
    description='Generic interface to HDF5 simulation files.',
    url='',
    author='Kyle Oman',
    author_email='koman@astro.rug.nl',
    license='',
    packages=['simfiles'],
    install_requires=['numpy', 'astropy', 'h5py'],
    include_package_data=True,
    zip_safe=False
)
