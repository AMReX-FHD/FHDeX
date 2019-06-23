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

    // Construct normal vector
    m_length = m_link.vectorLength();
    m_normal = m_link;
    if (m_length > 0)
        m_normal*= 1./m_length;

}
