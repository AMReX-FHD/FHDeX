#ifndef _IBMarkerMD_H_
#define _IBMarkerMD_H_


#include <AMReX.H>
#include <AMReX_RealVect.H>


using namespace amrex;


/********************************************************************************
 *                                                                              *
 * NAMESPACE immbdy_geom: a namspace containing classes to characterize the     *
 *                        geometry of immersed boundaries                       *
 *                                                                              *
 *******************************************************************************/


namespace immbdy_geom{

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

};



/********************************************************************************
 *                                                                              *
 * NAMESPACE immbdy_md: a namespace containign function to time-step (MD) the   *
 *                      immersed boundary markers                               *
 *                                                                              *
 *******************************************************************************/



namespace immbdy_md {

typedef enum {
    ARG_r_m,
    ARG_r,
    ARG_r_p //,
} NDERIV__coordinate;


using namespace immbdy_geom;

// Bending rigidity model
void add_bending_forces(Edge & e_ref, Edge & e, Real k, Real cos_theta_0);

void bending_f(      RealVect & f,       RealVect & f_p,       RealVect & f_m,
               const RealVect & r, const RealVect & r_p, const RealVect & r_m,
               Real k, Real cos_theta_0);


// Active driving model
Real UW(const RealVect & r_m, const RealVect & r, const RealVect & r_p,
        const RealVect & u, Real theta);

Real ndrUW(const RealVect & r_m, const RealVect & r, const RealVect & r_p,
           const RealVect & u, Real theta,
           NDERIV__coordinate arg, const RealVect & dx, Real delta);

void driving_f(      RealVect & f,       RealVect & f_p,       RealVect & f_m,
               const RealVect & r, const RealVect & r_p, const RealVect & r_m,
               const RealVect & u, Real theta, Real k);
}


#endif
