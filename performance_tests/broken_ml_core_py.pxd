
cdef extern from "broken_ml_core.h":
    cdef cppclass Vect:
        Vect (double,double,double) except +
        Vect() except +
        double x,y,z

    Vect _get_broken_ml_dev_core(
                    Vect fixed_s_pos,
                    Vect fixed_s_vec,
                    Vect orig_coords0,
                    Vect orig_coords1,
                    Vect orig_twist,
                    Vect s_twist,
                    double r1, double u1, double v1,
                    double u, double v, double t)

    void testdot()
