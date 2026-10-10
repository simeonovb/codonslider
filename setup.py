from setuptools import setup

setup(
    name='codonslider',
    version='2.0',
    packages=['codonslider'],
    url='https://github.com/simeonovb/codonslider',
    license='MIT',
    author='Boris Simeonov',
    author_email='bsimeonovb@gmail.com',
    description='Bioinformatic tool for codon stretch detection in CCDS sequences',
    install_requires = [
        'pandas',
        'biopython'
    ],
    python_requires='>=3.9',
)

