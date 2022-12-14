from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Complete set of modules for running HypoDD'
LONG_DESCRIPTION = 'I will add later...'

# Setting up
setup(
       # the name must match the folder name !
        name="PyHypoDD", 
        version=VERSION,
        author="Saeed SoltaniMoghadam",
        author_email="<saeed.sltm@gmail.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages
        keywords=['python', 'first package'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: Linux :: LINUX",
        ]
)