o1 = genmodule.object()
o1.name = 'SAINVObject'
o1.abbrev = 'SAINVObject'
o1.methodlist = ['precon']
o1.funclist = ['new', 'tp_dealloc', 'tp_getattr']
o1.typelist = []

m = genmodule.module()
m.name = 'sparslab'
m.abbrev = 'sparslab'
m.methodlist = ['sainv']
m.objects = [o1]

