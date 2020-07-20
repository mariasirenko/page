from setuptools import setup

setup(
    name='page',
    url='https://github.com/mariasirenko4/page',
    author='Maria Sirenko',
    author_email='sirenkom@sloankettering.edu',
    packages=['page'],
    install_requires=['numpy','pandas', 'pysam', 'matplotlib', 'seaborn', 'os', 'chileup'],
    version='0.1',
    license='MIT',
    description='Genotyping of single cell RNA-seq data from 10X libraries and Genotyping of Transcriptomes. ',
)