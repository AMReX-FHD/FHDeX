#include <IBMarkerMD.H>




Edge::Edge(Vertex & v1, Vertex & v2, Real length_0)
    : m_start(v1), m_end(v2), m_length_0(length_0) {

    Update();
}



Edge::Edge(Vertex & v1, Vertex & v2) : m_start(v1), m_end (v2) {

    RealVect link_init = m_end.r - m_start.r;
    Real length_init = link_init.vectorLength();

    Edge(v1, v2, length_init);
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




void assign_bending_forces(Edge & e_ref, Edge & e,
                           Real k, Real cos_theta_0, Real cos_theta) {


    //___________________________________________________________________________
    // edge lengths (and squares)
    Real l_ref = e_ref.length();
    Real l_e   = e.length();

    Real l2M = l_ref*l_ref;
    Real l2P = l_e*l_e;


    //___________________________________________________________________________
    // r, r_p, r_m <= vectors representing the central, previous (m) and next (p)
    //                vertex positions
    RealVect r, r_p, r_m;

    // vcpy(& r, & e->start->r);
    r = e.start().r;

    // vcpy(& r_p, & r);
    // vmuladd_ip(&r_p, e->length, &(e->normal));
    r_p = e.end().r;

    // vcpy(& r_m, & r);
    // vmuladd_ip(&r_m, -e_ref->length, &(e_ref->normal));
    r_m = e_ref.start().r;


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
    // compute bending forces

    // Real fx2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(x - xM))/l2M
    //         + (cos_theta*(x - xP))/l2P - (-2*x + xM + xP)/(sqrt(l2M)*sqrt(l2P)));

    // Real fPx2 = -2*(-cos_theta + cos_theta_0)*((-x + xM)/(sqrt(l2M)*sqrt(l2P))
    //         + (cos_theta*(-x + xP))/l2P);

    // Real fMx2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(-x + xM))/l2M
    //         + (-x + xP)/(sqrt(l2M)*sqrt(l2P)));


    Real fx2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(x - xM))/l2M
            + (cos_theta*(x - xP))/l2P - (-2*x + xM + xP)/(l_ref*l_e));

    Real fPx2 = -2*(-cos_theta + cos_theta_0)*((-x + xM)/(l_ref*l_e)
            + (cos_theta*(-x + xP))/l2P);

    Real fMx2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(-x + xM))/l2M
            + (-x + xP)/(l_ref*l_e));


#if (AMREX_SPACEDIM > 1)
    // Real fy2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(y - yM))/l2M
    //         + (cos_theta*(y - yP))/l2P - (-2*y + yM + yP)/(sqrt(l2M)*sqrt(l2P)));

    // Real fPy2 = -2*(-cos_theta + cos_theta_0)*((-y + yM)/(sqrt(l2M)*sqrt(l2P))
    //         + (cos_theta*(-y + yP))/l2P);

    // Real fMy2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(-y + yM))/l2M
    //         + (-y + yP)/(sqrt(l2M)*sqrt(l2P)));

    Real fy2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(y - yM))/l2M
            + (cos_theta*(y - yP))/l2P - (-2*y + yM + yP)/(l_ref*l_e));

    Real fPy2 = -2*(-cos_theta + cos_theta_0)*((-y + yM)/(l_ref*l_e)
            + (cos_theta*(-y + yP))/l2P);

    Real fMy2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(-y + yM))/l2M
            + (-y + yP)/(l_ref*l_e));
#endif

#if (AMREX_SPACEDIM > 2)
    // Real fz2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(z - zM))/l2M
    //         + (cos_theta*(z - zP))/l2P - (-2*z + zM + zP)/(sqrt(l2M)*sqrt(l2P)));

    // Real fPz2 = -2*(-cos_theta + cos_theta_0)*((-z + zM)/(sqrt(l2M)*sqrt(l2P))
    //         + (cos_theta*(-z + zP))/l2P);

    // Real fMz2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(-z + zM))/l2M
    //         + (-z + zP)/(sqrt(l2M)*sqrt(l2P)));

    Real fz2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(z - zM))/l2M
            + (cos_theta*(z - zP))/l2P - (-2*z + zM + zP)/(l_ref*l_e));

    Real fPz2 = -2*(-cos_theta + cos_theta_0)*((-z + zM)/(l_ref*l_e)
            + (cos_theta*(-z + zP))/l2P);

    Real fMz2 = -2*(-cos_theta + cos_theta_0)*((cos_theta*(-z + zM))/l2M
            + (-z + zP)/(l_ref*l_e));
#endif


    // apply bending stiffness

    RealVect f   = RealVect(AMREX_D_DECL(k*fx2/2,  k*fy2/2,  k*fz2/2));
    RealVect f_p = RealVect(AMREX_D_DECL(k*fPx2/2, k*fPy2/2, k*fPz2/2));
    RealVect f_m = RealVect(AMREX_D_DECL(k*fMx2/2, k*fMy2/2, k*fMz2/2));

    // RealVect f, f_p, f_m;

    // f.x = k*fx2/2;
    // f_p.x = k*fPx2/2;
    // f_m.x = k*fMx2/2;

    // f.y = k*fy2/2;
    // f_p.y = k*fPy2/2;
    // f_m.y = k*fMy2/2;

    // f.z = k*fz2/2;
    // f_p.z = k*fPz2/2;
    // f_m.z = k*fMz2/2;


    //___________________________________________________________________________
    // Add bending forces to the vertex data

    // vadd_ip(&(e->start->f), &f);
    // vadd_ip(&(e->end->f), &f_p);
    // vadd_ip(&(e_ref->start->f), &f_m);

    e.start().f     += f;
    e.end().f       += f_p;
    e_ref.start().f += f_m;

}
