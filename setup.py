#!/usr/bin/env python

from distutils.core import setup, Extension
import glob
import os, socket

# default settings
library_dirs_list = []
superlu_defs = [('USE_VENDOR_BLAS',1)]
f77_defs = []

hostname = socket.gethostname()
if hostname in ['bree', 'brokoli', 'givens2', 'givens4', 'nedelec', 'gondor']:
    # SuSE Linux 8.x and 9.0
    libraries_list = ['lapack', 'blas', 'g2c']
elif hostname == 'maxwell':
    # AMD Opteron 'x86-64' architecture using ACML
    #
    # use --install-purelib=/usr/lib64/python2.2/site-packages 
    # and --install-platlib=/usr/lib64/python2.2/site-packages options
    # when installing
    libraries_list = ['acml', 'g2c']
    library_dirs_list = ['/opt/acml/gnu64']
elif hostname == 'sim0':
    # Linux sim0 2.2.19-7.0.16enterprise #1 SMP Wed Mar 13 13:23:22 EST 2002 i686 unknown
    library_dirs_list = ['/home/geus/linux/lib']
    libraries_list = ['lapack', 'blas', 'g2c']
elif hostname == 'stardust':
    # HP-UX stardust B.11.11 U 9000/800 3761215035 unlimited-user license
    library_dirs_list = ['/home/infk/geus/lib/pa20_64']
    libraries_list = ['lapack']
elif hostname == 'zuse':
    # SunOS zuse 5.6 Generic_105181-21 sun4u sparc SUNW,Ultra-Enterprise
    library_dirs_list = ['/software/SunOS/5.X/opt/SUNWspro/WS6/lib']
    libraries_list = ['F77', 'sunperf', 'fui', 'fsu', 'sunmath']
elif hostname == 'Rivendell':
    # WinXP, using VC++ and Intel MKL
    #
    # Uses BLAS and LAPACK from the Intel Math Kernel Library 5.2
    # DLL directory must be in the PATH
    library_dirs_list = [r'C:\Program Files\Intel\MKL\ia32\lib']
    libraries_list = ['mkl_c_dll']
    superlu_defs = [("NO_TIMER",1), ("NoChange",1), ('USE_VENDOR_BLAS',1)]
    f77_defs = [("NOF77UNDERSCORE",1)]
elif hostname == 'nedelec-vmware':
    # Win32 using MinGW
    libraries_list = ['lapack', 'blas', 'g2c']
    superlu_defs = [("NO_TIMER",1), ('USE_VENDOR_BLAS',1)]
elif hostname == 'psw283.psi.ch':
    # OSF1 psw283.psi.ch V4.0 1530 alpha
    library_dirs_list = ['/data/geus/lib']
    libraries_list = ['dxml']
    

ext_modules = [Extension('spmatrix', ['Src/spmatrixmodule.c']),
               Extension('itsolvers', ['Src/itsolversmodule.c',
                                       'Src/pcg.c',
                                       'Src/minres.c',
                                       'Src/qmrs.c',
                                       'Src/bicgstab.c',
                                       'Src/cgs.c'],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         define_macros=f77_defs),
               Extension('precon',  [os.path.join('Src', 'preconmodule.c')],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         define_macros=f77_defs),
               Extension('superlu', [os.path.join('Src', 'superlumodule.c'),
                                     "superlu/dcolumn_bmod.c",
                                     "superlu/dcolumn_dfs.c",
                                     "superlu/dcomplex.c",
                                     "superlu/scomplex.c",
                                     "superlu/dcopy_to_ucol.c",
                                     "superlu/dgscon.c",
                                     "superlu/dgsequ.c",
                                     "superlu/dgsrfs.c",
                                     "superlu/dgssv.c",
                                     "superlu/dgssvx.c",
                                     "superlu/dgstrf.c",
                                     "superlu/dgstrs.c",
                                     "superlu/dlacon.c",
                                     "superlu/dlamch.c",
                                     "superlu/dlangs.c",
                                     "superlu/dlaqgs.c",
                                     "superlu/dmemory.c",
                                     "superlu/colamd.c",
                                     "superlu/dpanel_bmod.c",
                                     "superlu/dpanel_dfs.c",
                                     "superlu/dpivotL.c",
                                     "superlu/dpivotgrowth.c",
                                     "superlu/dpruneL.c",
                                     "superlu/dreadhb.c",
                                     "superlu/dsnode_bmod.c",
                                     "superlu/dsnode_dfs.c",
                                     "superlu/dsp_blas2.c",
                                     "superlu/dsp_blas3.c",
                                     "superlu/superlu_timer.c",
                                     "superlu/dutil.c",
                                     "superlu/dzsum1.c",
                                     "superlu/get_perm_c.c",
                                     "superlu/icmax1.c",
                                     "superlu/izmax1.c",
                                     "superlu/lsame.c",
                                     "superlu/memory.c",
                                     "superlu/mmd.c",
                                     "superlu/relax_snode.c",
                                     "superlu/sp_coletree.c",
                                     "superlu/sp_ienv.c",
                                     "superlu/sp_preorder.c",
                                     "superlu/util.c",
                                     "superlu/xerbla.c"],
                         define_macros=superlu_defs,
                         include_dirs=["superlu"],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list),
               Extension('jdsym', [os.path.join('Src', 'jdsymmodule.c')],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         define_macros=f77_defs),
               ]

execfile(os.path.join('Lib', 'pysparse_version.py'))

setup(name = 'pysparse',
      version = version,
      description = 'Python Sparse Matrix Package',
      author = 'Roman Geus',
      author_email = 'roman@geus.ch',
      license = 'BSD style',
      url = 'http://www.geus.ch',
      packages = [''],
      package_dir = {'': 'Lib'},
      include_dirs = ['Include'],
      headers = glob.glob(os.path.join ("Include","pysparse","*.h")),
      ext_modules = ext_modules
      )
