"""
Preconditioners.
"""

__all__ = filter(lambda s:not s.startswith('_'), dir())

from precon import ssor
from precon import jacobi
