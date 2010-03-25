#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    import numpy
    import fnmatch
    import os
    from numpy.distutils.misc_util import Configuration

    config = Configuration('pysparse', parent_package, top_path)

    config.add_subpackage('sparse')
    config.add_subpackage('eigen')
    config.add_subpackage('direct')
    config.add_subpackage('itsolvers')
    config.add_subpackage('precon')
    config.add_subpackage('tools')
    config.add_subpackage('misc')
    #config.add_data_dir('tests')

    # Pysparse header files
    config.add_data_dir(os.path.join(config.top_path,'pysparse','include'))

    config.make_config_py()
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
