#ifndef _IBMarkerMD_H_
#define _IBMarkerMD_H_


#include <AMReX.H>
#include <AMReX_RealVect.H>


using namespace amrex;



struct Vertex {
    RealVect r;
    RealVect f;
    int  name;
    bool fixed;

    Vertex(const RealVect & r_in);

    void Clear();
};



class Edge {

public:

    Edge(Vertex & v1, Vertex & v2, Real length_0);
    Edge(Vertex & v1, Vertex & v2);

    void Update();

    // Do not let vertices be mutable because changes to Vertex::r can affect
    // Edge's internals
    const Vertex & start() const {return m_start;};
    const Vertex & end() const {return m_end;};
    const RealVect & link() const {return m_link;};
    const RealVect & normal() const {return m_normal;};
    Real length() const {return m_length;};
    Real length_0() const {return m_length_0;};

    // Vertex forces are mutable because they do not affect Edge's internals
    RealVect & start_f() const {return m_start.f;};
    RealVect & end_f() const {return m_end.f;};

private:

    Vertex & m_start;
    Vertex & m_end;

    RealVect m_normal;
    RealVect m_link;
    Real m_length;
    Real m_length_0;
};



// struct Curve {
//  Edge * edges;
//  int length;
//  int capacity;
// };



void add_bending_forces(Edge & e_ref, Edge & e, Real k, Real cos_theta_0);

void bending_f(      RealVect & f,       RealVect & f_p,       RealVect & f_m,
               const RealVect & r, const RealVect & r_p, const RealVect & r_m,
               Real k, Real cos_theta_0);





#endif
