from setuptools import setup, find_packages
import sys

if sys.version_info[:2] < (3, 6):
    sys.stdout.write('Python 3.6 or later is required\n')
    sys.exit(1)

setup(
    name='trex',
    use_scm_version={'write_to': 'src/trex/_version.py'},
    setup_requires=['setuptools_scm'],  # Support pip versions that don't know about pyproject.toml
    author='',
    author_email='',
    url='',
    description='simultaneous lineage TRacking and EXpression profiling of single cells using RNA-seq',
    long_description='',
    license='MIT',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    install_requires=[
        'pysam',
        'numpy',
        'loompy',
        'xopen>=0.5.0',
        'tinyalign',
        'pandas',
    ],
    python_requires='>=3.6',
    entry_points={'console_scripts': ['trex = trex.__main__:main']},
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
