#!/usr/bin/env python
"""
PySparse: A Fast Sparse Matrix Library for Python

Pysparse is a fast sparse matrix library for Python. It provides several sparse
matrix storage formats and conversion methods. It also implements a number of
iterative solvers, preconditioners, and interfaces to efficient factorization
packages. Both low-level and high-level interfaces are available, each with
different strengths. PySparse is distributed under the FreeBSD license.

R. Geus    <hamsel@users.sf.net>
D. Orban   <d-orban@users.sf.net>
D. Wheeler <wd15@users.sf.net>
"""

DOCLINES = __doc__.split("\n")

import os
import sys
from pkg_resources import require

try:
    import setuptools  # If you want to enable 'python setup.py develop'
    #pass
except:
    print 'setuptools module not found.'
    print "Install setuptools if you want to enable 'python setup.py develop'."

require('numpy>=1.2')

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Programming Language :: C
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS :: MacOS X
Natural Language :: English
"""

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_include_dirs(numpy.get_include())
    config.add_include_dirs(os.path.join(config.top_path,'pysparse','include'))
    config.add_subpackage('pysparse')

    # Set config.version
    config.get_version(os.path.join('pysparse','version.py'))

    return config

##Hacked from numpy
def svn_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        import subprocess
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['svn', 'info'])
    except OSError:
        print(" --- Could not run svn info --- ")
        return ""

    import re
    r = re.compile('Revision: ([0-9]+)')
    svnver = ""

    out = out.decode()

    for line in out.split('\n'):
        m = r.match(line.strip())
        if m:
            svnver = m.group(1)

    if not svnver:
        print("Error while parsing svn version")

    return svnver

def get_version():
    version = '1.2'
    release = False
    if not release:
        version += '-dev' + svn_version()
    return version

def write_version_py(filename=os.path.join('pysparse', 'version.py')):
    cnt = """
## THIS FILE IS GENERATED FROM PYSPARSE SETUP.PY
version = '%(version)s'
"""

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': get_version()})
    finally:
        a.close()

def setup_package():

    from numpy.distutils.core import setup, Extension
    from numpy.distutils.misc_util import Configuration

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0,local_path)
    sys.path.insert(0,os.path.join(local_path,'pysparse')) # to retrieve version

    try:
        setup(
            name = 'pysparse',
            author = "Roman Geus, Dominique Orban, Daniel Wheeler",
            author_email = "{hamsel,d-orban,wd15}@sf.net,",
            maintainer = "PySparse Developers",
            maintainer_email = "{hamsel,d-orban,wd15}@sf.net,",
            summary = "Fast sparse matrix library for Python",
            description = DOCLINES[0],
            long_description = "\n".join(DOCLINES[2:]),
            url = "pysparse.sf.net",
            download_url = "sf.net/projects/pysparse",
            license = 'BSD-style',
            classifiers=filter(None, CLASSIFIERS.split('\n')),
            platforms = ["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
            install_requires=['numpy>=1.2'],
            configuration=configuration,
            )
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    write_version_py()
    setup_package()
