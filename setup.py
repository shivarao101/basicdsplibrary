from setuptools import setup, find_packages 
with open("README.md", "r", encoding="utf-8") as f:
    readme = f.read()

with open( "CHANGELOG.txt", "r", encoding="utf-8") as f:
    changelog = f.read()

long_description = readme + "\n\n" + changelog

classifiers = ['Development Status :: 3 - Alpha',
               'Environment :: Console',
               'Intended Audience :: Education',
               'License :: OSI Approved :: MIT License',
               'Natural Language :: English',
               'Operating System :: OS Independent',
               'Programming Language :: Python',
               'Topic :: Scientific/Engineering'
]

setup(
    name='basicdsplibrary',
    version='0.0.12',
    description='Basic DSP library without numpy',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/shivarao101/vcet-dsp-bec502',
    author='Shivaprasad',
    author_email='shivarao101@gmail.com',
    license='MIT',
    classifiers=classifiers,
    keywords='DFT, IDFT, DFS, DTFT, DCT, DWT, Radix-2 DIT & DIF, Overlap save and add',
    packages=find_packages(),
)