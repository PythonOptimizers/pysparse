#!/usr/bin/env python

from distutils import sysconfig, version
from distutils.core import setup, Extension
import glob
import os, socket
import sys

# default settings
library_dirs_list= []
libraries_list = ['lapack', 'blas', 'g2c']
superlu_defs = [('USE_VENDOR_BLAS',1)]
if sys.platform == 'win32':
    superlu_defs += [('NO_TIMER', 1)]

f77_defs = []
linky=[]
compily=[]

# Specify whether to link against user's SuperLU library
use_users_superlu = False #True
umfpack_defs = [('DINT', 1), ('NBLAS', 1)] # No BLAS. Ok if using your system's BLAS.
#umfpack_defs = [('DINT', 1), ('CBLAS', 1)] # with atlas c-blas (http://math-atlas.sf.net)

if use_users_superlu:
    # Specify location of include files
    superlu_include_dirs = ['/Users/dpo/local/linalg/SuperLU/SuperLU_3.1/SRC']
    # AND specify ONE of the following:
    # 1) Location of source files (overrides linking with existing library)
    superlu_src_dir = '/Users/dpo/local/linalg/SuperLU/SuperLU_3.1/SRC'
    # 2) Location of precompiled library
    superlu_lib_dir = ['']
    superlu_libraries = ['']

# Specify whether to link against user's UMFPACK library
use_users_umfpack = False #True
if use_users_umfpack:
    # Specify location of include files
    umfpack_include_dirs = ['/Users/dpo/local/linalg/UMFPACK/UMFPACK/Include']
    # AND specify ONE of the following:
    # 1) Location of source files (overrides linking with existing library)
    umfpack_src_dir = '/Users/dpo/local/linalg/UMFPACK/UMFPACK/Source'
    # 2) Location of precompiled library
    umfpack_lib_dir = ['']
    umfpack_libraries = ['']

    # Do the same for AMD
    amd_include_dirs = ['/Users/dpo/local/linalg/UMFPACK/AMD/Include']
    # AND specify ONE of the following:
    # 1) Location of source files (overrides linking with existing library)
    amd_src_dir = '/Users/dpo/local/linalg/UMFPACK/AMD/Source'
    # 2) Location of precompiled library
    amd_lib_dir = ['']
    amd_libraries = ['']
    
else:
    umfpack_libraries = []
    umfpack_include_dirs = ['amd', 'umfpack']
    umfpack_library_dirs = []

    #umfpack_libraries = ['atlas', 'cblas', 'm']
    #umfpack_include_dirs = ['amd', 'umfpack'] # you may need to set this to find the atlas
    #umfpack_library_dirs = []

# Obtain Python version
pyver = sysconfig.get_python_version()  # Returns something like '2.5'
current = version.StrictVersion(pyver)  # Object that can be compored
target = version.StrictVersion('2.6')
newpy_defs = []
if current >= target:
    newpy_defs = [('LENFUNC_OK', 1)]    # For PyMappingMethods types

## numpy / Numeric settings
try:
    import numpy
    use_numpy = True
except:
    try:
        import Numeric
        use_numpy = False
    except:
        raise ImportError, "Failed to import either numpy or Numeric"

package_name = 'pysparse'
if use_numpy:
    numerix_macro = [('NUMPY', '1')]
    numerix_include_dirs = [numpy.get_include()]
else:
    numerix_macro = []
    numerix_include_dirs = []

if sys.platform == 'darwin':
    superlu_defs = [('USE_VENDOR_BLAS',1)]
    library_dirs_list = ['/System/Library/Frameworks']
    libraries_list = []
    f77_defs = []
    
    # the following 'linky' arguments must not be concatenated together into a single
    # string, c.f. <http://mail.python.org/pipermail/distutils-sig/2003-December/003532.html>
    
    if sys.exec_prefix == '/sw':
        # fink python
        linky=["-faltivec","-framework","vecLib","-bundle_loader","/sw/bin/python"]
    else:
        # Apple python
        #linky=["-faltivec","-framework","Accelerate"]
        linky=["-framework","Accelerate"]
        # The python Framework build is compiled with (and propagates to all
        # library builds) the '-fno-common' flag. Nobody seems to know why.
        # (c.f. <https://sourceforge.net/tracker/?func=detail&atid=105470&aid=768306&group_id=5470>)
        # This flag wreaks havoc with the nightmarishly circular declarations
        # in the itsolvers module. We reset it by appending this flag:
        compily=["-fcommon"]

elif sys.platform == 'linux2':
    ## fix for Fedora core 4, 'g2c' doesn't exist and isn't required
    if 'redhat-release' in os.listdir('/etc'):
        f = open('/etc/redhat-release', 'r')
        if 'release 4' in f.read():
            libraries_list = ['lapack', 'blas']
        f.close()

from distutils.core import Command
    
ext_modules = [Extension(package_name + '.spmatrix', ['Src/spmatrixmodule.c'],
                         define_macros=numerix_macro + newpy_defs),
               Extension(package_name + '.itsolvers', ['Src/itsolversmodule.c',
                                       'Src/pcg.c',
                                       'Src/gmres.c',
                                       'Src/minres.c',
                                       'Src/qmrs.c',
                                       'Src/bicgstab.c',
                                       'Src/cgs.c'],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         define_macros=f77_defs + numerix_macro,
                         extra_compile_args=compily,
                         extra_link_args=linky),
               Extension(package_name + '.precon',  [os.path.join('Src', 'preconmodule.c')],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         define_macros=f77_defs + numerix_macro,
                         extra_compile_args=compily,
                         extra_link_args=linky),
               Extension(package_name + '.jdsym', [os.path.join('Src', 'jdsymmodule.c')],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         define_macros=f77_defs + numerix_macro,
                         extra_compile_args=compily,
                         extra_link_args=linky)
               ]

# Schedule UMFPACK
if not use_users_umfpack:
    ext_modules += [Extension(package_name + '.umfpack',
                              sources=[os.path.join('Src', 'umfpackmodule.c'),
                                     'amd/amd_aat.c',
                                     'amd/amd_1.c',
                                     'amd/amd_2.c',
                                     'amd/amd_dump.c',
                                     'amd/amd_postorder.c',
                                     'amd/amd_post_tree.c',
                                     'amd/amd_defaults.c',
                                     'amd/amd_order.c',
                                     'amd/amd_control.c',
                                     'amd/amd_info.c',
                                     'amd/amd_valid.c',
                                     'umfpack/umf_analyze.c',
                                     'umfpack/umf_apply_order.c',
                                     'umfpack/umf_colamd.c',
                                     'umfpack/umf_free.c',
                                     'umfpack/umf_fsize.c',
                                     'umfpack/umf_is_permutation.c',
                                     'umfpack/umf_malloc.c',
                                     'umfpack/umf_realloc.c',
                                     'umfpack/umf_report_perm.c',
                                     'umfpack/umf_singletons.c',
                                     'umfpack/umfpack_timer.c',
                                     'umfpack/umfpack_tictoc.c',
                                     'umfpack/umf_lhsolve.c',
                                     'umfpack/umf_uhsolve.c',
                                     'umfpack/umf_triplet_map_nox.c',
                                     'umfpack/umf_triplet_nomap_x.c',
                                     'umfpack/umf_triplet_nomap_nox.c',
                                     'umfpack/umf_triplet_map_x.c',
                                     'umfpack/umf_assemble_fixq.c',
                                     'umfpack/umf_assemble.c',
                                     'umfpack/umf_blas3_update.c',
                                     'umfpack/umf_build_tuples.c',
                                     'umfpack/umf_create_element.c',
                                     'umfpack/umf_dump.c',
                                     'umfpack/umf_extend_front.c',
                                     'umfpack/umf_garbage_collection.c',
                                     'umfpack/umf_get_memory.c',
                                     'umfpack/umf_init_front.c',
                                     'umfpack/umf_kernel.c',
                                     'umfpack/umf_kernel_init.c',
                                     'umfpack/umf_kernel_wrapup.c',
                                     'umfpack/umf_local_search.c',
                                     'umfpack/umf_lsolve.c',
                                     'umfpack/umf_ltsolve.c',
                                     'umfpack/umf_mem_alloc_element.c',
                                     'umfpack/umf_mem_alloc_head_block.c',
                                     'umfpack/umf_mem_alloc_tail_block.c',
                                     'umfpack/umf_mem_free_tail_block.c',
                                     'umfpack/umf_mem_init_memoryspace.c',
                                     'umfpack/umf_report_vector.c',
                                     'umfpack/umf_row_search.c',
                                     'umfpack/umf_scale_column.c',
                                     'umfpack/umf_set_stats.c',
                                     'umfpack/umf_solve.c',
                                     'umfpack/umf_symbolic_usage.c',
                                     'umfpack/umf_transpose.c',
                                     'umfpack/umf_tuple_lengths.c',
                                     'umfpack/umf_usolve.c',
                                     'umfpack/umf_utsolve.c',
                                     'umfpack/umf_valid_numeric.c',
                                     'umfpack/umf_valid_symbolic.c',
                                     'umfpack/umf_grow_front.c',
                                     'umfpack/umf_start_front.c',
                                     'umfpack/umf_2by2.c',
                                     'umfpack/umf_store_lu.c',
                                     'umfpack/umf_scale.c',
                                     'umfpack/umfpack_wsolve.c',
                                     'umfpack/umfpack_col_to_triplet.c',
                                     'umfpack/umfpack_defaults.c',
                                     'umfpack/umfpack_free_numeric.c',
                                     'umfpack/umfpack_free_symbolic.c',
                                     'umfpack/umfpack_get_numeric.c',
                                     'umfpack/umfpack_get_lunz.c',
                                     'umfpack/umfpack_get_symbolic.c',
                                     'umfpack/umfpack_numeric.c',
                                     'umfpack/umfpack_qsymbolic.c',
                                     'umfpack/umfpack_report_control.c',
                                     'umfpack/umfpack_report_info.c',
                                     'umfpack/umfpack_report_matrix.c',
                                     'umfpack/umfpack_report_numeric.c',
                                     'umfpack/umfpack_report_perm.c',
                                     'umfpack/umfpack_report_status.c',
                                     'umfpack/umfpack_report_symbolic.c',
                                     'umfpack/umfpack_report_triplet.c',
                                     'umfpack/umfpack_report_vector.c',
                                     'umfpack/umfpack_solve.c',
                                     'umfpack/umfpack_symbolic.c',
                                     'umfpack/umfpack_transpose.c',
                                     'umfpack/umfpack_triplet_to_col.c',
                                     'umfpack/umfpack_scale.c',
                                     'umfpack/umfpack_load_numeric.c',
                                     'umfpack/umfpack_save_numeric.c',
                                     'umfpack/umfpack_load_symbolic.c',
                                     'umfpack/umfpack_save_symbolic.c'],
                         define_macros=umfpack_defs + numerix_macro,
                         include_dirs=umfpack_include_dirs,
                         libraries=umfpack_libraries,
                         extra_compile_args=compily,
                         library_dirs=umfpack_library_dirs)
               ]
else:
    if umfpack_src_dir.strip() != '':
        umf_src_files = [ 'umf_2by2.c', 'umf_analyze.c', 'umf_apply_order.c', 'umf_assemble.c', 'umf_blas3_update.c', 'umf_build_tuples.c', 'umf_colamd.c', 'umf_create_element.c', 'umf_dump.c', 'umf_extend_front.c', 'umf_free.c', 'umf_fsize.c', 'umf_garbage_collection.c', 'umf_get_memory.c', 'umf_grow_front.c', 'umf_init_front.c', 'umf_is_permutation.c', 'umf_kernel.c', 'umf_kernel_init.c', 'umf_kernel_wrapup.c', 'umf_local_search.c', 'umf_lsolve.c', 'umf_ltsolve.c', 'umf_malloc.c', 'umf_mem_alloc_element.c', 'umf_mem_alloc_head_block.c', 'umf_mem_alloc_tail_block.c', 'umf_mem_free_tail_block.c', 'umf_mem_init_memoryspace.c', 'umf_multicompile.c', 'umf_realloc.c', 'umf_report_perm.c', 'umf_report_vector.c', 'umf_row_search.c', 'umf_scale.c', 'umf_scale_column.c', 'umf_set_stats.c', 'umf_singletons.c', 'umf_solve.c', 'umf_start_front.c', 'umf_store_lu.c', 'umf_symbolic_usage.c', 'umf_transpose.c', 'umf_triplet.c', 'umf_tuple_lengths.c', 'umf_usolve.c', 'umf_utsolve.c', 'umf_valid_numeric.c', 'umf_valid_symbolic.c', 'umfpack_col_to_triplet.c', 'umfpack_defaults.c', 'umfpack_free_numeric.c', 'umfpack_free_symbolic.c', 'umfpack_get_determinant.c', 'umfpack_get_lunz.c', 'umfpack_get_numeric.c', 'umfpack_get_symbolic.c', 'umfpack_global.c', 'umfpack_load_numeric.c', 'umfpack_load_symbolic.c', 'umfpack_numeric.c', 'umfpack_qsymbolic.c', 'umfpack_report_control.c', 'umfpack_report_info.c', 'umfpack_report_matrix.c', 'umfpack_report_numeric.c', 'umfpack_report_perm.c', 'umfpack_report_status.c', 'umfpack_report_symbolic.c', 'umfpack_report_triplet.c', 'umfpack_report_vector.c', 'umfpack_save_numeric.c', 'umfpack_save_symbolic.c', 'umfpack_scale.c', 'umfpack_solve.c', 'umfpack_symbolic.c', 'umfpack_tictoc.c', 'umfpack_timer.c', 'umfpack_transpose.c', 'umfpack_triplet_to_col.c' ]
        umf_src_files = map(lambda s: os.path.join(umfpack_src_dir, s), umf_src_files)
        src_files = umf_src_files + [os.path.join('Src', 'umfpackmodule.c')]
    else:
        src_files = [os.path.join('Src', 'umfpackmodule.c')]
        library_dirs_list += umfpack_lib_dir
        libraries_list += umfpack_libraries

    if amd_src_dir.strip() != '':
        amd_src_files = [ 'amd_1.c', 'amd_2.c', 'amd_aat.c', 'amd_control.c', 'amd_defaults.c', 'amd_dump.c', 'amd_global.c', 'amd_info.c', 'amd_order.c', 'amd_post_tree.c', 'amd_postorder.c', 'amd_preprocess.c', 'amd_valid.c' ]
        amd_src_files = map(lambda s: os.path.join(amd_src_dir, s), amd_src_files)
        src_files += amd_src_files

    else:
        library_dirs_list += amd_lib_dir
        libraries_list += amd_libraries

    ext_modules += [ Extension(package_name + '.umfpack',
                         src_files,
                         define_macros=umfpack_defs + numerix_macro,
                         include_dirs=umfpack_include_dirs + amd_include_dirs,
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         extra_compile_args=compily,
                         extra_link_args=linky) ]

# Schedule SuperLU
if not use_users_superlu:
    ext_modules += [ Extension(package_name + '.superlu',
                             [os.path.join('Src', 'superlumodule.c'),
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
                         define_macros=superlu_defs + numerix_macro,
                         include_dirs=["superlu"],
                         library_dirs=library_dirs_list,
                         extra_compile_args=compily,
                         libraries=libraries_list,extra_link_args=linky) ]
else:
    if superlu_src_dir.strip() != '':
        src_files = ['superlu_timer.c', 'util.c', 'memory.c', 'get_perm_c.c', 'mmd.c', 'sp_coletree.c', 'sp_preorder.c', 'sp_ienv.c', 'relax_snode.c', 'heap_relax_snode.c', 'colamd.c', 'lsame.c', 'xerbla.c', 'dlacon.c', 'dlamch.c', 'dgssv.c', 'dgssvx.c', 'dsp_blas2.c', 'dsp_blas3.c', 'dgscon.c', 'dlangs.c', 'dgsequ.c', 'dlaqgs.c', 'dpivotgrowth.c', 'dgsrfs.c', 'dgstrf.c', 'dgstrs.c', 'dcopy_to_ucol.c', 'dsnode_dfs.c', 'dsnode_bmod.c', 'dpanel_dfs.c', 'dpanel_bmod.c', 'dreadhb.c', 'dcolumn_dfs.c', 'dcolumn_bmod.c', 'dpivotL.c', 'dpruneL.c', 'dmemory.c', 'dutil.c']
        src_files = map(lambda s: os.path.join(superlu_src_dir, s), src_files)
        src_files += [os.path.join('Src', 'superlu3module.c')]
    else:
        src_files = [os.path.join('Src', 'superlu3module.c')]
        library_dirs_list += superlu_lib_dir
        libraries_list += superlu_libraries
    ext_modules += [ Extension(package_name + '.superlu',
                         src_files,
                         define_macros=superlu_defs + numerix_macro,
                         include_dirs=superlu_include_dirs,
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         extra_compile_args=compily,
                         extra_link_args=linky) ]

execfile(os.path.join('Lib', '__version__.py'))

setup(name = 'pysparse',
      version = __version__,
      description = 'Python Sparse Matrix Package',
      author = 'Roman Geus',
      author_email = 'roman@geus.ch',
      license = 'BSD style',
      url = 'http://www.geus.ch',
      packages = [package_name],
      package_dir = {package_name: 'Lib'},
      package_data = {package_name: package_data},
      include_dirs = ['Include'] + numerix_include_dirs,
      headers = glob.glob(os.path.join ("Include","pysparse","*.h")),
      ext_modules = ext_modules
      )
