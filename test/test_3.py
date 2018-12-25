# Simple inheritance case
class super1():
    def __new__(typ, *args, **kwargs):
        obj = object.__new__(typ)
        obj.attr1 = []
        print('super1_new')
        return obj


class derived1(super1):
    def __new__(typ, *args, **kwargs):
        obj = object.__new__(typ)
        print(typ)
        print(obj)
        obj.attr1 = ['dammy']
        print('derived1_new')
        return obj

    def __init__(self, arg4, **kwargs):
        self.attr4 = arg4
        self.attr5 = kwargs['arg5']
        print('derived1_init')


if '__main__' == __name__:
    print('before derived1() called')
    d1 = derived1(222, arg5=333)
    d1.attr1.append(111)
    print(d1.attr1, d1.attr4, d1.attr5)
    print(isinstance(d1, super1))
    d2 = derived1(22222, arg5=33333)
