from __future__ import print_function
from orthonormal_test import create_orthonormal_basis, create_orthonormal_basis2
import numpy as np
from forgi.threedee.utilities.vector import magnitude
import forgi


def create_orthonormal_basis_numpy(vec1, vec2):
    vec3 = np.cross(vec1, vec2)
    assert abs(magnitude(vec3)-1)<10**-10, magnitude(vec3)
    return np.array([vec1,vec2,vec3])


if __name__ == "__main__":
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
    rna, = forgi.load_rna("../test/forgi/threedee/data/1GID_A.cg")
    class virtual_stat:
         r1 = 5
         u1 = 1.2
         v1 = 0.4
         u  = 0.7
         v  = 0.25
         t  = 0.92
    setup = ("import brokenml1; import brokenml2; import cppvect; from __main__ import rna, virtual_stat")
    t_broken1 = timeit("brokenml1.get_broken_ml_deviation(rna, 'm3', 's8', virtual_stat)", setup=setup, number=10000)
    t_broken2 = timeit("brokenml2.get_broken_ml_deviation(rna, 'm3', 's8', virtual_stat)", setup=setup, number=10000)
    t_broken2 = timeit("cppvect.get_broken_ml_deviation(rna, 'm3', 's8', virtual_stat)", setup=setup, number=10000)
    print(t_broken1, t_broken2)
