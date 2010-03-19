#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    import numpy
    import os
    import ConfigParser
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, NotFoundError

    # Read relevant PySparse-specific configuration options.
    pysparse_config = ConfigParser.SafeConfigParser()
    pysparse_config.read(os.path.join(top_path, 'site.cfg'))
    umfpack_libdir = pysparse_config.get('UMFPACK', 'umfpack_libdir')
    umfpack_include = pysparse_config.get('UMFPACK', 'umfpack_include')
    superlu_libdir = pysparse_config.get('SuperLU', 'superlu_libdir')
    superlu_include = pysparse_config.get('SuperLU', 'superlu_include')
    amd_libdir = pysparse_config.get('AMD', 'amd_libdir')
    amd_include = pysparse_config.get('AMD', 'amd_include')

    config = Configuration('direct', parent_package, top_path)

    # Get BLAS info from site.cfg
    blas_info = get_info('blas_opt',0)
    if not blas_info:
        print 'No blas info found'

    umfpack_src= ['umfpackmodule.c']
    config.add_extension(
        name='umfpack',
        sources=[os.path.join('src',name) for name in umfpack_src],
        libraries=['amd', 'umfpack'],
        library_dirs=[amd_libdir] + [umfpack_libdir],
        include_dirs=['src'] + [umfpack_include] + [amd_include],
        extra_info=blas_info,
        )

    superlu_src= ['superlu3module.c']
    config.add_extension(
        name='superlu',
        sources=[os.path.join('src',name) for name in superlu_src],
        libraries=['superlu'],
        library_dirs=[superlu_libdir],
        include_dirs=['src'] + [superlu_include],
        extra_info=blas_info,
        )

    config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
