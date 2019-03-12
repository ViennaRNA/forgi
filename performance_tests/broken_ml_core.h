#include<array>
#include<vector>

const double PI = 3.141592653589793;

class Vect{
    public: double x,y,z;
    Vect(double x, double y, double z);
    Vect();
    ~Vect();
    Vect cross(const Vect &b);
    void normalize();
    double distance_to(const Vect &b);
    void print();

};

class Basis{
    public: std::array<std::array<double,3>,3> coords;
    Basis(Vect a, Vect b, Vect c);
    Basis(double ax, double bx, double cx, double ay, double by, double cy, double az, double bz, double cz);
    Vect apply_to(const Vect &a);
    void print();
};

Vect transform_coords(Basis b, double r, double u, double v);
Basis rotation_matrix_y(const double theta);
Basis rotation_matrix_z(const double theta);
Vect twist2_orient_from_stem1_1(Basis stem1_basis,
                                double u, double v, double t);
Vect spherical_polar_to_cartesian(double r, double u, double v);
Basis create_orthonormal_basis(Vect a, Vect b);
Vect operator+(const Vect &a, const Vect &b);
Vect operator-(const Vect &a, const Vect &b);
double operator*(const Vect &a, const Vect &b); // dot product
Vect operator*(const Vect &a, const double &t);
Basis operator*(const Basis &a, const Basis &t);
double vec_angle(Vect a, Vect b);
Vect _get_broken_ml_dev_core(
                            Vect fixed_s_pos,
                            Vect fixed_s_vec,
                            Vect orig_coords0,
                            Vect orig_coords1,
                            Vect orig_twist,
                            Vect s_twist,
                            double r1, double u1, double v1,
                            double u, double v, double t);

void testdot();
