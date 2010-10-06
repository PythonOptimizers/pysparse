#!/usr/bin/env python

def getoption(config, section, option):
    try:
        val = config.get(section,option)
    except:
        val = None
    return val

def configuration(parent_package='',top_path=None):
    import numpy
    import os
    import ConfigParser
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, NotFoundError

    # Read relevant PySparse-specific configuration options.
    pysparse_config = ConfigParser.SafeConfigParser()
    pysparse_config.read(os.path.join(top_path, 'site.cfg'))
    dflt_lib_dirs = getoption( pysparse_config, 'DEFAULT', 'library_dirs')
    if dflt_lib_dirs is None:
        dflt_lib_dirs = []
    dflt_libs = getoption(pysparse_config, 'DEFAULT', 'libraries')
    if dflt_libs is None:
        dflt_libs = []

    print 'Using dflt_lib_dirs = ', dflt_lib_dirs
    print 'Using dflt_libs = ', dflt_libs

    config = Configuration('eigen', parent_package, top_path)

    # Get BLAS info from site.cfg
    blas_info = get_info('blas_opt',0)
    if not blas_info:
        blas_info = get_info('blas',0)
        if not blas_info:
            print 'No blas info found'
    print 'Eigen:: Using BLAS info:' ; print blas_info

    # Get LAPACK info from site.cfg
    lapack_info = get_info('lapack_opt',0)
    if not lapack_info:
        lapack_info = get_info('lapack',0)
        if not lapack_info:
            print 'No lapack info found'
    print 'Eigen:: Using LAPACK info:' ; print lapack_info

    jdsym_src = ['jdsymmodule.c']
    config.add_extension(
        name='jdsym',
        sources=[os.path.join('src',name) for name in jdsym_src],
        libraries=dflt_libs,
        library_dirs=dflt_lib_dirs,
        include_dirs=['src'],
        extra_info=[blas_info, lapack_info],
        )

    config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
