from setuptools import setup

setup(
    name='simfiles',
    version='1.0',
    description='Generic interface to HDF5 simulation files.',
    url='',
    author='Kyle Oman',
    author_email='koman@astro.rug.nl',
    license='GNU GPL v3',
    packages=['simfiles'],
    install_requires=['numpy', 'astropy', 'h5py'],
    include_package_data=True,
    zip_safe=False
)
