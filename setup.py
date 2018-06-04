import sys
from setuptools import setup, find_packages

if sys.version_info[:2] < (3, 6):
    sys.stdout.write('Python 3.6 or later is required\n')
    sys.exit(1)

setup(
    name='jfrisen1801',
    version='0.1',
    author='',
    author_email='',
    url='',
    description='',
    long_description='',
    license='MIT',
    py_modules=['180425_lineage_loom'],
    install_requires=['pysam', 'numpy', 'loompy'],
    entry_points={'console_scripts': ['180425_lineage_loom = 180425_lineage_loom:main']},
    classifiers=[
        'Development Status :: 1 - Planning',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
