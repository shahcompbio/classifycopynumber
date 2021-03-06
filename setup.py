"""
Copy number classifier
"""
from setuptools import find_packages, setup

install_requires = [
    'click',
    'pandas',
    'PyYAML',
    'seaborn',
    'tables',
    'scgenome@git+https://github.com/shahcompbio/scgenome.git#egg=scgenome',
    'wgs_analysis@git+https://github.com/amcpherson/wgs_analysis.git#egg=wgs_analysis',
]

dependency_links = [
    'https://github.com/shahcompbio/scgenome.git#egg=scgenome',
    'https://github.com/amcpherson/wgs_analysis.git#egg=wgs_analysis',
]

setup(
    name='classifycopynumber',
    version='0.1.0',
    url='https://github.com/amcpherson/classifycopynumber',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    description='Copy number classifier',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=install_requires,
    dependency_links=dependency_links,
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
