from setuptools import setup


with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name = 'ipathapy', 
    version='0.0.1', 
    description = 'A package to create metabolic pathmaps with ipath3', 
    long_description=long_description,
    author = 'Lucas Diedrich', 
    author_email = 'lucas.diedrich@embl.de', 
    license = 'MIT'
    zip_safe = False
)