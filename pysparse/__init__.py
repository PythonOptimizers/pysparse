"PySparse: A Fast Sparse Matrix Library for Python"

__docformat__ = 'restructuredtext'

# Imports
from numpy._import_tools import PackageLoader
try:
    from version import version as __version__
except ImportError:
    Warning("Run setup.py to generate version number.")
    __version__ = 'undefined'
    
from sparse import spmatrix
#from sparse import *
from misc import get_include

pkgload = PackageLoader()
pkgload(verbose=False,postpone=True)

if __doc__:
    __doc__ += """

Available subpackages
---------------------
"""
if __doc__:
    __doc__ += pkgload.get_pkgdocs()

__all__ = filter(lambda s: not s.startswith('_'), dir())
__all__ += '__version__'

__doc__ += """

Miscellaneous
-------------

    __version__  :  pysparse version string
"""

from pysparse.misc import Deprecated

class _superlu:
    @Deprecated('Use pysparse.direct.superlu.factorize instead.')
    def factorize(self, *args, **kwargs):
        import pysparse.direct.superlu
        self.factorizeFnc = pysparse.direct.superlu.factorize
        return self.factorizeFnc(*args, **kwargs)
    
superlu = _superlu()
