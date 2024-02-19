"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup
import os

if not os.path.exists('scripts'):
    os.makedirs('scripts')

scripts = ['polaris-gen.in']
for script in scripts:
    # Read in the file
    with open(script, 'r') as file:
        filedata = file.read()

    # Replace the target string
    polaris_path = os.path.abspath(os.getcwd().replace('tools', '')) + os.path.sep
    filedata = filedata.replace('@POLARIS_PATH@', polaris_path)

    # Write the file out again
    with open('scripts/' + script.replace('.in', ''), 'w') as file:
        file.write(filedata)


# here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
# with open(os.path.join(here, 'README'), encoding='utf-8') as f:
#     long_description = f.read()

# https://packaging.python.org/en/latest/specifications/core-metadata/
# https://setuptools.pypa.io/en/latest/references/keywords.html
setup(
    name='polaristools',

    version='2.0.0',

    description='PolarisTools is a Python toolkit to create grid files for POLARIS',

    # long_description=long_description,

    url='https://portia.astrophysik.uni-kiel.de/polaris/',

    author='Robert Brauer',

    author_email='robert.brauer@cea.fr',

    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',

        'License :: OSI Approved :: GNU General Public License (GPL)',

        'Programming Language :: Python :: 3 :: Only',
    ],

    keywords='astronomy, astrophysics, monte-carlo, radiative-transfer, simulation, polarization',

    packages=['polaris_tools_modules', 'polaris_tools_custom'],

    python_requires='>=3.6',

    install_requires=['numpy'],

    project_urls={
        'Homepage': 'https://portia.astrophysik.uni-kiel.de/polaris',
        'GitHub': 'https://github.com/polaris-MCRT/POLARIS',
    },

    scripts=['scripts/polaris-gen'],
)
