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
};



class Edge {

public:

    Edge(Vertex & v1, Vertex & v2, Real length_0);
    Edge(Vertex & v1, Vertex & v2);

    void Update();

    Vertex & start() const {return m_start;};
    Vertex & end() const {return m_end;};
    const RealVect & link() const {return m_link;};
    const RealVect & normal() const {return m_normal;};
    Real length() const {return m_length;};
    Real length_0() const {return m_length_0;};

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



void assign_bending_forces(Edge & e_ref, Edge & e, Real k,
                           Real cos_theta_0, Real cos_theta);



void assign_bending_forces(Edge & e_ref, Edge & e, Real k);


#endif
