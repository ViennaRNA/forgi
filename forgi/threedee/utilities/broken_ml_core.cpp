#include <math.h>       /* sqrt, pow */
#include <vector>
#include "broken_ml_core.h"
#include <iostream>
#include <assert.h>     /* assert */
Vect::Vect ( double x, double y, double z){
    this->x=x;
    this->y=y;
    this->z=z;
}
Vect::Vect(){}

void Vect::normalize(){
    double mag = sqrt(pow(this->x,2)+pow(this->y, 2)+pow(this->z,2));
    this->x/=mag;
    this->y/=mag;
    this->z/=mag;
}

Vect::~Vect(){}

Vect Vect::cross(const Vect &b){
    return Vect(
        this->y*b.z-this->z*b.y,
        this->z*b.x-this->x*b.z,
        this->x*b.y-this->y*b.x
    );
}

double Vect::distance_to(const Vect &b){
    return sqrt( pow(this->x-b.x,2)+ pow(this->y-b.y,2) + pow(this->z-b.z,2) );

}

Basis::Basis(Vect a, Vect b, Vect c){
    this->coords[0][0]=a.x;
    this->coords[0][1]=a.y;
    this->coords[0][2]=a.z;
    this->coords[1][0]=b.x;
    this->coords[1][1]=b.y;
    this->coords[1][2]=b.z;
    this->coords[2][0]=c.x;
    this->coords[2][1]=c.y;
    this->coords[2][2]=c.z;

}

Basis::Basis(double ax, double bx, double cx, double ay, double by, double cy, double az, double bz, double cz){
    this->coords[0][0]=ax;
    this->coords[1][0]=ay;
    this->coords[2][0]=az;
    this->coords[0][1]=bx;
    this->coords[1][1]=by;
    this->coords[2][1]=bz;
    this->coords[0][2]=cx;
    this->coords[1][2]=cy;
    this->coords[2][2]=cz;

}

Vect Basis::apply_to(const Vect& a){
    return Vect (
        this->coords[0][0]*a.x+this->coords[1][0]*a.y+this->coords[2][0]*a.z,
        this->coords[0][1]*a.x+this->coords[1][1]*a.y+this->coords[2][1]*a.z,
        this->coords[0][2]*a.x+this->coords[1][2]*a.y+this->coords[2][2]*a.z
    );
}

Vect spherical_polar_to_cartesian(double r, double u, double v){
    return Vect(
        r*sin(u)*cos(v),
        r*sin(u)*sin(v),
        r*cos(u)
    );
}

Vect transform_coords(Basis b, double r, double u, double v){
    Vect xyz = spherical_polar_to_cartesian(r,u,v);
    return b.apply_to(xyz);

}

Vect twist2_orient_from_stem1_1(Basis stem1_basis,
                                double u, double v, double t){
    Vect twist2_new = Vect(0., cos(t), sin(t));
    Basis rot_mat1 = rotation_matrix_z(v);
    Basis rot_mat2 = rotation_matrix_y(u - PI / 2.);

    Basis rot_mat = rot_mat2 * rot_mat1;
    twist2_new = rot_mat.apply_to(twist2_new);

    Vect twist2_new_basis = stem1_basis.apply_to(twist2_new);

    return twist2_new_basis;
}

Basis rotation_matrix_y(const double theta){
    double c = cos(theta);
    double s = sin(theta);
    return Basis(
        c, 0.,-s,
        0.,1.,0.,
        s, 0., c
    );
}
Basis rotation_matrix_z(const double theta){
    double c = cos(theta);
    double s = sin(theta);
    return Basis(
        c ,s ,0.,
        -s,c ,0.,
        0.,0.,1.
    );
}
Basis create_orthonormal_basis(Vect a, Vect b){
    a.normalize();
    b.normalize();
    // std::cout<<"normalixed a ";
    // a.print();
    Vect c=a.cross(b);
    // std::cout<<"c: ";
    // c.print();
    return Basis(a,b,c);
}

double vec_angle(Vect a, Vect b){
    a.normalize();
    b.normalize();
    double d=a*b;
    if (d<-1.){d=-1.;}
    if (d>1.){d=1.;}
    return acos(d);
}

Vect operator+(const Vect &a, const Vect&b){
    return Vect(a.x+b.x, a.y+b.y, a.z+b.z);
}

Vect operator-(const Vect &a, const Vect&b){
    return Vect(a.x-b.x, a.y-b.y, a.z-b.z);
}

double operator*(const Vect &a, const Vect&b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}
Vect operator*(const Vect &a, const double &t){
    return Vect(a.x*t, a.y*t, a.z*t);
}
Basis operator*(const Basis &a, const Basis &b){
    return Basis(
            a.coords[0][0]*b.coords[0][0]+a.coords[0][1]*b.coords[1][0]+a.coords[0][2]*b.coords[2][0],
            a.coords[0][0]*b.coords[0][1]+a.coords[0][1]*b.coords[1][1]+a.coords[0][2]*b.coords[2][1],
            a.coords[0][0]*b.coords[0][2]+a.coords[0][1]*b.coords[1][2]+a.coords[0][2]*b.coords[2][2],
            a.coords[1][0]*b.coords[0][0]+a.coords[1][1]*b.coords[1][0]+a.coords[1][2]*b.coords[2][0],
            a.coords[1][0]*b.coords[0][1]+a.coords[1][1]*b.coords[1][1]+a.coords[1][2]*b.coords[2][1],
            a.coords[1][0]*b.coords[0][2]+a.coords[1][1]*b.coords[1][2]+a.coords[1][2]*b.coords[2][2],
            a.coords[2][0]*b.coords[0][0]+a.coords[2][1]*b.coords[1][0]+a.coords[2][2]*b.coords[2][0],
            a.coords[2][0]*b.coords[0][1]+a.coords[2][1]*b.coords[1][1]+a.coords[2][2]*b.coords[2][1],
            a.coords[2][0]*b.coords[0][2]+a.coords[2][1]*b.coords[1][2]+a.coords[2][2]*b.coords[2][2]
    );
}

void Vect::print(){
        std::cout<<this->x<<", "<<this->y<<", "<<this->z<<"\n";
}

Vect _get_broken_ml_dev_core(
                            Vect fixed_s_pos,
                            Vect fixed_s_vec,
                            Vect orig_coords0,
                            Vect orig_coords1,
                            Vect orig_twist,
                            Vect s_twist,
                            double r1, double u1, double v1,
                            double u, double v, double t){
    // std::cout<<"Creating Fixed stem basis from\n";
    // fixed_s_pos.print();
    // fixed_s_vec.print();
    Basis fixed_stem_basis = create_orthonormal_basis(fixed_s_vec, s_twist);
    // std::cout<<"Fixed stem basis ";
    //fixed_stem_basis.print();
    Vect vbulge_vec = transform_coords(fixed_stem_basis, r1, u1, v1);
    Vect vstem_vec = transform_coords(fixed_stem_basis, 1., u, v);
    Vect vstem_twist = twist2_orient_from_stem1_1(fixed_stem_basis,
                                                   u,v,t);

    // std::cout<<"vbulge_vec";
    // vbulge_vec.print();
    // std::cout<<"vstem_vec";
    // vstem_vec.print();
    // std::cout<<"vstem_twist";
    // vstem_twist.print();

    vstem_vec = vstem_vec*5;
    Vect vstem_coords0 = fixed_s_pos + vbulge_vec;
    Vect vstem_coords1 = vstem_coords0 + vstem_vec;

    Vect orig_stem_vec = orig_coords1 - orig_coords0;
    Vect true_bulge_vec = orig_coords0 - fixed_s_pos;

    double pos_dev = orig_coords0.distance_to(vstem_coords0);
    double ang_dev = vec_angle(vstem_vec, orig_stem_vec);
    double twist_dev = vec_angle(orig_twist, vstem_twist);
    return Vect(pos_dev, ang_dev, twist_dev);
}

void testdot(){
    Basis b1 = Basis(1.,2.,3.,4.,5.,6.,7.,8.,9);
    std::cout<<"B1:\n";
    b1.print();
    Basis b2 = Basis(1.,2.,1.,2.,1.,2.,1.,2.,3.);
    std::cout<<"B2:\n";
    b2.print();

    Basis b3 = b1*b2;
    std::cout<<"B3:\n";
    b3.print();
    assert (b3.coords[0][0]==8);
    assert (b3.coords[1][0]==20);
    assert (b3.coords[1][2]==32);

    Vect a=Vect(1.,2.,5.);
    Vect b=b1.apply_to(a);
    std::cout<<"Vector b";
    b.print();

}

void Basis::print(){
    std::cout<<this->coords[0][0]<<" "<<this->coords[0][1]<<" "<<this->coords[0][2]<<"\n";
    std::cout<<this->coords[1][0]<<" "<<this->coords[1][1]<<" "<<this->coords[1][2]<<"\n";
    std::cout<<this->coords[2][0]<<" "<<this->coords[2][1]<<" "<<this->coords[2][2]<<"\n";
}
