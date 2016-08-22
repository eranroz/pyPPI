from distutils.core import setup

setup(
    name='pyPPI',
    version='0.1',
    packages=['pyPPI', 'pyPPI.alignment', 'pyPPI.surfaceComplementarity'],
    scripts=['bin/setupPpiDb.py'],
    url='https://github.com/eranroz/pyPPI',
    license='MIT',
    author='eranroz',
    author_email='eranroz@cs.huji.ac.il',
    description='Protein-Protein interactions calculations', requires=['requests', 'pymysql']
)
