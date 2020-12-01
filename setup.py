"""
Copy number classifier
"""
from setuptools import find_packages, setup

dependencies = ['click']
# classifycopynumber/metadata/additional_genes.csv
# classifycopynumber/metadata/antigen_presenting_genes.csv
# classifycopynumber/metadata/cancer_gene_census.csv
# classifycopynumber/metadata/Census_ampThu Apr 16 15_35_36 2020.csv
# classifycopynumber/metadata/Census_delsThu Apr 16 15_36_24 2020.csv
# classifycopynumber/metadata/hr_genes.txt
# package_data = package_data={'amp_genes': 'classifycopynumber/metadata/Census_ampThu Apr 16 15_35_36 2020.csv', 
#     'del_genes':'classifycopynumber/metadata/Census_delsThu Apr 16 15_36_24 2020.csv',
#     'additional_genes': 'classifycopynumber/metadata/additional_genes.csv',
#     'antigen_presenting_genes':'classifycopynumber/metadata/antigen_presenting_genes.csv',
#     'hr_genes':'classifycopynumber/metadata/hr_genes.txt'}

setup(
    name='classifycopynumber',
    version='0.1.0',
    url='https://github.com/amcpherson/classifycopynumber',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    description='Copy number classifier',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    # package_data={'amp_genes': ['metadata/Census_ampThu Apr 16 15_35_36 2020.csv']},
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'classifycopynumber = classifycopynumber.cli:main',
        ],
    },
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'Development Status :: 1 - Planning',
        # 'Development Status :: 2 - Pre-Alpha',
        # 'Development Status :: 3 - Alpha',
        'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        # 'Development Status :: 6 - Mature',
        # 'Development Status :: 7 - Inactive',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
