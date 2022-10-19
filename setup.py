from setuptools import find_packages, setup


with open('requirements.txt', 'r') as f:
    requirements = f.read().split('\n')

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name = 'ipathapy', 
    version='0.0.1', 
    description = 'A package to create metabolic pathmaps with ipath3', 
    long_description=long_description,
    author = 'Lucas Diedrich', 
    author_email = 'lucas.diedrich@embl.de', 
    license = 'MIT',
    packages=find_packages(['ipathapy']),
    include_package_data=True,
    install_requires=requirements,


    zip_safe = False
)