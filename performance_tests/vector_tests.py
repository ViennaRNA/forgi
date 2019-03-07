from orthonormal_test import create_orthonormal_basis, create_orthonormal_basis2
import numpy as np
from forgi.threedee.utilities.vector import magnitude

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
             "from orthonormal_test import create_orthonormal_basis, create_orthonormal_basis2")
    t_np = timeit("create_orthonormal_basis_numpy(a,b)", setup=setup)
    t_cy1 = timeit("create_orthonormal_basis(a,b)", setup=setup)
    t_cy2 = timeit("create_orthonormal_basis2(a,b)", setup=setup)
    print(t_np, t_cy1, t_cy2)
