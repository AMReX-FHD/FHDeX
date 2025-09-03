#include <cfloat>

#include <AMReX.H>
#include <AMReX_Print.H>

#include <IBMarkerMD.H>



namespace immbdy_geom {

Vertex::Vertex(const RealVect & r_in) : r(r_in) {

    Clear();
}



void Vertex::Clear() {

    f = RealVect{AMREX_D_DECL(0., 0., 0.)};
}



Edge::Edge(Vertex & v1, Vertex & v2, Real length_0)
    : m_start(v1), m_end(v2), m_length_0(length_0) {

    Update();
}



Edge::Edge(Vertex & v1, Vertex & v2) : m_start(v1), m_end (v2) {

    RealVect link_init = m_end.r - m_start.r;
    m_length_0 = link_init.vectorLength();

    Update();
}



void Edge::Update() {

    // Construct link vector
    m_link = m_end.r - m_start.r;

    // Edge length
    m_length = m_link.vectorLength();

    // Construct normal vector
    m_normal = m_link;
    if (m_length > 0)
        m_normal *= 1./m_length;
}

};



namespace immbdy_md {

using namespace immbdy_geom;


void add_bending_forces(Edge & e_ref, Edge & e, Real k, Real cos_theta_0) {


    //___________________________________________________________________________
    // edge lengths (and squares)

    Real l_ref = e_ref.length();
    Real l_e   = e.length();

    Real l2M = l_ref*l_ref;
    Real l2P = l_e*l_e;


    //___________________________________________________________________________
    // r, r_p, r_m <= vectors representing the central, previous (m) and next (p)
    //                vertex positions

    const RealVect & r = e.start().r;
    const RealVect & r_p = e.end().r;
    const RealVect & r_m = e_ref.start().r;


    //___________________________________________________________________________
    // angle between two edges

    Real cos_theta = e_ref.normal().dotProduct(e.normal());


    //___________________________________________________________________________
    // cartesian coordinates of vertices (m: previous, p: next)

    Real x = r[0], xM = r_m[0], xP = r_p[0];
#if (AMREX_SPACEDIM > 1)
    Real y = r[1], yM = r_m[1], yP = r_p[1];
#endif
#if (AMREX_SPACEDIM > 2)
    Real z = r[2], zM = r_m[2], zP = r_p[2];
#endif


    //___________________________________________________________________________
    // compute bending forces analytically

    Real fx2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(x - xM))/l2M
            + (cos_theta*(x - xP))/l2P - (-2*x + xM + xP)/(l_ref*l_e));

    Real fPx2 = -2*(-cos_theta + cos_theta_0)*((-x + xM)/(l_ref*l_e)
            + (cos_theta*(-x + xP))/l2P);

    Real fMx2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(-x + xM))/l2M
            + (-x + xP)/(l_ref*l_e));


#if (AMREX_SPACEDIM > 1)
    Real fy2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(y - yM))/l2M
            + (cos_theta*(y - yP))/l2P - (-2*y + yM + yP)/(l_ref*l_e));

    Real fPy2 = -2*(-cos_theta + cos_theta_0)*((-y + yM)/(l_ref*l_e)
            + (cos_theta*(-y + yP))/l2P);

    Real fMy2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(-y + yM))/l2M
            + (-y + yP)/(l_ref*l_e));
#endif

#if (AMREX_SPACEDIM > 2)
    Real fz2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(z - zM))/l2M
            + (cos_theta*(z - zP))/l2P - (-2*z + zM + zP)/(l_ref*l_e));

    Real fPz2 = -2*(-cos_theta + cos_theta_0)*((-z + zM)/(l_ref*l_e)
            + (cos_theta*(-z + zP))/l2P);

    Real fMz2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(-z + zM))/l2M
            + (-z + zP)/(l_ref*l_e));
#endif


    //___________________________________________________________________________
    // Add bending forces to the vertex data

    RealVect f;
    RealVect f_p;
    RealVect f_m;

    f[0]   = k*fx2/2;
    f_p[0] = k*fPx2/2;
    f_m[0] = k*fMx2/2;
#if (AMREX_SPACEDIM > 1)
    f[1]   = k*fy2/2;
    f_p[1] = k*fPy2/2;
    f_m[1] = k*fMy2/2;
#endif
#if (AMREX_SPACEDIM > 2)
    f[2]   = k*fz2/2;
    f_p[2] = k*fPz2/2;
    f_m[2] = k*fMz2/2;
#endif


    e.start_f()     += f;
    e.end_f()       += f_p;
    e_ref.start_f() += f_m;
}



void bending_f(      RealVect & f,       RealVect & f_p,       RealVect & f_m,
               const RealVect & r, const RealVect & r_p, const RealVect & r_m,
               Real k, Real cos_theta_0) {

    // Construct vertices (note that Vertex::f == 0 due to constructor)
    Vertex v(r), v_p(r_p), v_m(r_m);

    // Set up topology
    Edge e_ref(v_m, v), e(v, v_p);

    // Add bending forces to vertices
    add_bending_forces(e_ref, e, k, cos_theta_0);


    // Add result to bending forces to output variables
    f   += e.start().f;
    f_p += e.end().f;
    f_m += e_ref.start().f;
}



Real UW(const RealVect & r_m, const RealVect & r, const RealVect & r_p,
        const RealVect & u, Real theta) {

    Real x = r[0], xM = r_m[0], xP = r_p[0], ux = u[0];
#if (AMREX_SPACEDIM > 1)
    Real y = r[1], yM = r_m[1], yP = r_p[1], uy = u[1];
#endif
#if (AMREX_SPACEDIM > 2)
    Real z = r[2], zM = r_m[2], zP = r_p[2], uz = u[2];
#endif

    Real cosTh = cos(theta);
#if (AMREX_SPACEDIM > 2)
    Real sinTh = sin(theta);
#endif

#if   (AMREX_SPACEDIM == 1)
#error incompatible with DIM == 1
#elif (AMREX_SPACEDIM == 2)

    // TODO: shouldn't this just be the 2D rotation matrix (around "z")?

    Real A1 = cosTh + (1 - cosTh)*ux*ux;
    Real A2 = (1 - cosTh)*ux*uy;

    Real B1 = (1 - cosTh)*ux*uy;
    Real B2 = cosTh + (1 - cosTh)*uy*uy;

    Real l_p = 1.;
    Real l_m = 1.;

    // NOTE: this might not be necessary
    // Real l_p = std::sqrt( (xP-x)*(xP-x) + (yP-y)*(yP-y) );
    // Real l_m = std::sqrt( (x-xM)*(x-xM) + (y-yM)*(y-yM) );
    // if (l_p == 0 ) l_p = 100*DBL_EPSILON;
    // if (l_m == 0 ) l_m = 100*DBL_EPSILON;

    Real Y1 = (xP-x)/l_p - ( A1*(x-xM) + B1*(y-yM) )/l_m;
    Real Y2 = (yP-y)/l_p - ( A2*(x-xM) + B2*(y-yM) )/l_m;


#elif (AMREX_SPACEDIM == 3)

    Real A1 = cosTh + (1 - cosTh)*ux*ux;
    Real A2 = (1 - cosTh)*ux*uy + sinTh*uz;
    Real A3 = -(sinTh*uy) + (1 - cosTh)*ux*uz;

    Real B1 = (1 - cosTh)*ux*uy - sinTh*uz;
    Real B2 = cosTh + (1 - cosTh)*uy*uy;
    Real B3 = sinTh*ux + (1 - cosTh)*uy*uz;

    Real C1 = sinTh*uy + (1 - cosTh)*ux*uz;
    Real C2 = -(sinTh*ux) + (1 - cosTh)*uy*uz;
    Real C3 = cosTh + (1 - cosTh)*uz*uz;

    Real l_p = 1.;
    Real l_m = 1.;

    // NOTE: this might not be necessary
    // Real l_p = std::sqrt( (xP-x)*(xP-x) + (yP-y)*(yP-y) + (zP-z)*(zP-z) );
    // Real l_m = std::sqrt( (x-xM)*(x-xM) + (y-yM)*(y-yM) + (z-zM)*(z-zM) );
    // if (l_p == 0 ) l_p = 100*DBL_EPSILON;
    // if (l_m == 0 ) l_m = 100*DBL_EPSILON;

    Real Y1 = (xP-x)/l_p - ( A1*(x-xM) + B1*(y-yM) + C1*(z-zM) )/l_m;
    Real Y2 = (yP-y)/l_p - ( A2*(x-xM) + B2*(y-yM) + C2*(z-zM) )/l_m;
    Real Y3 = (zP-z)/l_p - ( A3*(x-xM) + B3*(y-yM) + C3*(z-zM) )/l_m;

#else
#error incompatible with DIM > 3
#endif


#if   (AMREX_SPACEDIM == 2)
    return (Y1*Y1 + Y2*Y2);
#elif (AMREX_SPACEDIM == 3)
    return (Y1*Y1 + Y2*Y2 + Y3*Y3);
#endif
}


Real ndrUW(const RealVect & r_m, const RealVect & r, const RealVect & r_p,
           const RealVect & u, Real theta,
           NDERIV__coordinate arg, const RealVect & dx, Real delta) {

    RealVect r_dx;

    if      (arg == ARG_r_m) r_dx = r_m + dx;
    else if (arg == ARG_r)   r_dx = r   + dx;
    else if (arg == ARG_r_p) r_dx = r_p + dx;

    else {
        Abort();
    }

    Real uw_r_dx_p=0;

    if      (arg == ARG_r_m) uw_r_dx_p = UW(r_dx, r, r_p, u, theta);
    else if (arg == ARG_r)   uw_r_dx_p = UW(r_m, r_dx, r_p, u, theta);
    else if (arg == ARG_r_p) uw_r_dx_p = UW(r_m, r, r_dx, u, theta);

    else {
        Abort();
    }

    if      (arg == ARG_r_m) r_dx = r_m - dx;
    else if (arg == ARG_r)   r_dx = r   - dx;
    else if (arg == ARG_r_p) r_dx = r_p - dx;

    else {
        Abort();
    }

    Real uw_r_dx_m=0;

    if      (arg == ARG_r_m) uw_r_dx_m = UW(r_dx, r, r_p, u, theta);
    else if (arg == ARG_r)   uw_r_dx_m = UW(r_m, r_dx, r_p, u, theta);
    else if (arg == ARG_r_p) uw_r_dx_m = UW(r_m, r, r_dx, u, theta);

    else {
        Abort();
    }

    return (uw_r_dx_p - uw_r_dx_m)/(2*delta);
}


void driving_f(      RealVect & f,       RealVect & f_p,       RealVect & f_m,
               const RealVect & r, const RealVect & r_p, const RealVect & r_m,
               const RealVect & u, Real theta, Real Kw) {

    BL_PROFILE_VAR("immbdy_md::driving_f", BendingForce);

    Real delta = 100*DBL_EPSILON;

    RealVect dx = {AMREX_D_DECL(delta, 0.0, 0.0)};
    RealVect dy = {AMREX_D_DECL(0.0, delta, 0.0)};
#if (AMREX_SPACEDIM > 2)
    RealVect dz = {AMREX_D_DECL(0.0, 0.0, delta)};
#endif


    Real fx2  = ndrUW(r_m, r, r_p, u, theta, ARG_r,   dx, delta);
    Real fPx2 = ndrUW(r_m, r, r_p, u, theta, ARG_r_p, dx, delta);
    Real fMx2 = ndrUW(r_m, r, r_p, u, theta, ARG_r_m, dx, delta);

#if (AMREX_SPACEDIM > 1)
    Real fy2  = ndrUW(r_m, r, r_p, u, theta, ARG_r,   dy, delta);
    Real fPy2 = ndrUW(r_m, r, r_p, u, theta, ARG_r_p, dy, delta);
    Real fMy2 = ndrUW(r_m, r, r_p, u, theta, ARG_r_m, dy, delta);
#endif
#if (AMREX_SPACEDIM > 2)
    Real fz2  = ndrUW(r_m, r, r_p, u, theta, ARG_r,   dz, delta);
    Real fPz2 = ndrUW(r_m, r, r_p, u, theta, ARG_r_p, dz, delta);
    Real fMz2 = ndrUW(r_m, r, r_p, u, theta, ARG_r_m, dz, delta);
#endif

    RealVect f_loc;
    RealVect f_p_loc;
    RealVect f_m_loc;

    f_loc[0]   = -Kw*fx2/2;
    f_p_loc[0] = -Kw*fPx2/2;
    f_m_loc[0] = -Kw*fMx2/2;
#if (AMREX_SPACEDIM > 1)
    f_loc[1]   = -Kw*fy2/2;
    f_p_loc[1] = -Kw*fPy2/2;
    f_m_loc[1] = -Kw*fMy2/2;
#endif
#if (AMREX_SPACEDIM > 2)
    f_loc[2]   = -Kw*fz2/2;
    f_p_loc[2] = -Kw*fPz2/2;
    f_m_loc[2] = -Kw*fMz2/2;
#endif


    f   += f_loc;
    f_p += f_p_loc;
    f_m += f_m_loc;

    BL_PROFILE_VAR_STOP(BendingForce);
}



void rotate_z(Real & rx, Real & ry,
              Real   tx, Real   ty,
              Real theta) {

    Real ct = cos(theta);
    Real st = sin(theta);

    rx = ct*tx - st*ty;
    ry = st*tx + ct*ty;
}



void next_node_z(Real & nx, Real & ny,
                 Real   px, Real   py,
                 Real   tx, Real   ty,
                 Real ds) {

    nx = px + ds*tx;
    ny = py + ds*ty;
}

};
