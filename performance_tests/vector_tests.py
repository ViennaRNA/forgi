from __future__ import print_function
from orthonormal_test import create_orthonormal_basis, create_orthonormal_basis2
import numpy as np
from forgi.threedee.utilities.vector import magnitude
import forgi.threedee.utilities.vector as ftuv
import forgi
import cpp_vect, brokenml1
import sys
from cppvect import transposed_inverted
def create_orthonormal_basis_numpy(vec1, vec2):
    vec3 = np.cross(vec1, vec2)
    assert abs(magnitude(vec3)-1)<10**-10, magnitude(vec3)
    return np.array([vec1,vec2,vec3])




if __name__ == "__main__":
    print("=====  test inv ========")
    a=np.random.rand(3);
    a=a/ftuv.magnitude(a)
    b=ftuv.get_orthogonal_unit_vector(a)
    basis = ftuv.create_orthonormal_basis(a,b)
    print(np.linalg.inv(basis.T))
    print(transposed_inverted(basis))
    assert np.all(np.linalg.inv(basis.T)-transposed_inverted(basis)  < 10**-8)
    from timeit import timeit
    setup = ("import numpy as np; import forgi.threedee.utilities.vector as ftuv;"
             "a=np.random.rand(3); a=a/ftuv.magnitude(a);"
             "b=ftuv.get_orthogonal_unit_vector(a);"
             "basis = ftuv.create_orthonormal_basis(a,b);"
             "from cppvect import transposed_inverted;")
    t_cyt =  timeit("transposed_inverted(basis)", setup=setup)
    t_np  =  timeit("np.linalg.inv(basis.T)", setup=setup)
    print("cyt, np", t_cyt, t_np)
    print("====== testdot =========")
    print("Basis dot product")
    a=np.array([[1.,2,3],[4.,5,6],[7.,8,9]])
    b=np.array([[1.,2,1],[2.,1,2],[1.,2,3]])
    print(np.dot(a,b))
    cpp_vect.testdot1()
    print("Vect-basis-dot")
    print(np.dot(a.T,np.array([1.,2,5])))
    print("====== testdot done =========")

    print("broken ML deviation")
    rna, = forgi.load_rna("../test/forgi/threedee/data/1GID_A.cg")
    class virtual_stat:
         r1 = 1.
         u1 = 1.
         v1 = 1.
         u  = 1.
         v  = 1.
         t  = 1.

    a = brokenml1.get_broken_ml_deviation(rna, 'm3', 's8', virtual_stat)
    b = cpp_vect.get_broken_ml_deviation(rna, 'm3', 's8', virtual_stat)
    print(a,b)
    assert a==b
    print("TIMINGS")

    from timeit import timeit
    setup = ("import numpy as np; import forgi.threedee.utilities.vector as ftuv;"
             " a=np.random.rand(3); a=a/ftuv.magnitude(a);"
             "b=ftuv.get_orthogonal_unit_vector(a);"
             "from __main__ import create_orthonormal_basis_numpy;"
             "from orthonormal_test import create_orthonormal_basis, create_orthonormal_basis2, create_orthonormal_basis3")
    t_np = 0#timeit("create_orthonormal_basis_numpy(a,b)", setup=setup)
    t_cy1 = 0#timeit("create_orthonormal_basis(a,b)", setup=setup)
    t_cy2 = 0#timeit("create_orthonormal_basis2(a,b)", setup=setup)
    t_cy3 = 0#timeit("create_orthonormal_basis3(a,b)", setup=setup)
    print(t_np, t_cy1, t_cy2, t_cy3)

    setup = ("from brokenml2 import test_vec, test_np; import numpy as np; a=np.random.rand(3); b=np.random.rand(3)")
    t_vec =  timeit("test_vec(a,b)", setup=setup)
    t_np  =  timeit("test_vec(a,b)", setup=setup)

    print("Vect, np", t_vec, t_np)
    setup = ("import brokenml1; import brokenml2; import cpp_vect; from __main__ import rna, virtual_stat")
    t_broken1 = timeit("brokenml1.get_broken_ml_deviation(rna, 'm3', 's8', virtual_stat)", setup=setup, number=100000)
    t_broken3 = timeit("cpp_vect.get_broken_ml_deviation(rna, 'm3', 's8', virtual_stat)", setup=setup, number=100000)
    print(t_broken1, t_broken3)
