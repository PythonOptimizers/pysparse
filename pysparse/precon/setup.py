#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    import numpy
    import os
    import ConfigParser
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info, NotFoundError

    # Read relevant PySparse-specific configuration options.
    #pysparse_config = ConfigParser.SafeConfigParser()
    #pysparse_config.read(os.path.join(top_path, 'site.cfg'))
    #hsl_dir = pysparse_config.get('HSL', 'hsl_dir')

    config = Configuration('precon', parent_package, top_path)

    # Get BLAS info from site.cfg
    blas_info = get_info('blas_opt',0)
    if not blas_info:
        print 'No blas info found'

    precon_src = ['preconmodule.c']
    config.add_extension(
        name='precon',
        sources=[os.path.join('src',name) for name in precon_src],
        libraries=[],
        include_dirs=['src'],
        extra_info=blas_info,
        )

    config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
