import os
execfile(os.path.join(__path__[0], '__version__.py'))

try:
    import spmatrix
    import itsolvers
    import jdsym
    import precon
    import superlu
    import umfpack
    import sparray
except:
    pass
