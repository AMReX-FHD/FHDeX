#include "AMReX_PlotFileUtil.H"
#include "compressible_functions.H"
#include "compressible_functions_stag.H"
#include "common_functions.H"

#if defined(TURB)
void GetTurbQty(std::array< MultiFab, AMREX_SPACEDIM >& vel,
                std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                MultiFab& prim,
                MultiFab& eta,
                MultiFab& zeta,
                const amrex::Geometry& geom,
                Real& turbKE, Real& c_speed,
                Real& u_rms,
                Real& taylor_len, Real& taylor_Re, Real& taylor_Ma,
                Real& skew, Real& kurt,
                Real& eps_s, Real& eps_d, Real& eps_ratio,
                Real& kolm_s, Real& kolm_d, Real& kolm_t)

{
    BL_PROFILE_VAR("GetTurbQty()",GetTurbQty);

    Real dProb = (AMREX_SPACEDIM==2) ? n_cells[0]*n_cells[1] :
                                       n_cells[0]*n_cells[1]*n_cells[2];
    dProb = 1./dProb;

    // Setup temp MultiFabs
    std::array< MultiFab, AMREX_SPACEDIM > macTemp;
    MultiFab gradU;
    MultiFab eta_bulk_diss;
    MultiFab sound_speed;
    MultiFab ccTemp;
    MultiFab ccTempA;
    MultiFab ccTempDiv;
    std::array< MultiFab, NUM_EDGE > curlU;
    std::array< MultiFab, NUM_EDGE > eta_edge;
    std::array< MultiFab, NUM_EDGE > curlUtemp;
    AMREX_D_TERM(macTemp[0].define(convert(prim.boxArray(),nodal_flag_x), prim.DistributionMap(), 1, 1);,
                 macTemp[1].define(convert(prim.boxArray(),nodal_flag_y), prim.DistributionMap(), 1, 1);,
                 macTemp[2].define(convert(prim.boxArray(),nodal_flag_z), prim.DistributionMap(), 1, 1););
    gradU.define(prim.boxArray(),prim.DistributionMap(),AMREX_SPACEDIM,0);
    sound_speed.define(prim.boxArray(),prim.DistributionMap(),1,0);
    ccTemp.define(prim.boxArray(),prim.DistributionMap(),1,0);
    ccTempA.define(prim.boxArray(),prim.DistributionMap(),1,0);
    ccTempDiv.define(prim.boxArray(),prim.DistributionMap(),1,0);
    if (visc_type == 3) eta_bulk_diss.define(prim.boxArray(),prim.DistributionMap(),1,0);
#if (AMREX_SPACEDIM == 3)
    curlU[0].define(convert(prim.boxArray(),nodal_flag_xy), prim.DistributionMap(), 1, 0);
    curlU[1].define(convert(prim.boxArray(),nodal_flag_xz), prim.DistributionMap(), 1, 0);
    curlU[2].define(convert(prim.boxArray(),nodal_flag_yz), prim.DistributionMap(), 1, 0);
    eta_edge[0].define(convert(prim.boxArray(),nodal_flag_xy), prim.DistributionMap(), 1, 0);
    eta_edge[1].define(convert(prim.boxArray(),nodal_flag_xz), prim.DistributionMap(), 1, 0);
    eta_edge[2].define(convert(prim.boxArray(),nodal_flag_yz), prim.DistributionMap(), 1, 0);
#elif (AMREX_SPACEDIM == 2)
    curlU[0].define(convert(prim.boxArray(),nodal_flag_xy), prim.DistributionMap(), 1, 0);
    eta_edge[0].define(convert(prim.boxArray(),nodal_flag_xy), prim.DistributionMap(), 1, 0);
#endif

#if (AMREX_SPACEDIM == 3)
    curlUtemp[0].define(convert(prim.boxArray(),nodal_flag_xy), prim.DistributionMap(), 1, 0);
    curlUtemp[1].define(convert(prim.boxArray(),nodal_flag_xz), prim.DistributionMap(), 1, 0);
    curlUtemp[2].define(convert(prim.boxArray(),nodal_flag_yz), prim.DistributionMap(), 1, 0);
#elif (AMREX_SPACEDIM == 2)
    curlUtemp[0].define(convert(prim.boxArray(),nodal_flag_xy), prim.DistributionMap(), 1, 0);
#endif

    // Setup temp variables
    Real temp;
    Vector<Real> tempvec(3);
    Vector<Real> rhouu(3);
    Vector<Real> uu(3);
    Vector<Real> gradU2(3);
    Vector<Real> gradU3(3);
    Vector<Real> gradU4(3);
    Vector<Real> eps_s_vec(3); // solenoidal dissipation

    // turbulent kinetic energy
//   StagInnerProd(cumom,0,vel,0,macTemp,rhouu);
    {
        auto mask = cumom[0].OwnerMask(geom.periodicity());
        rhouu[0] = MultiFab::Dot(cumom[0],0,vel[0],0,1,0);
    }
    {
        auto mask = cumom[1].OwnerMask(geom.periodicity());
        rhouu[1] = MultiFab::Dot(cumom[1],0,vel[1],0,1,0);
    }
    {
        auto mask = cumom[2].OwnerMask(geom.periodicity());
        rhouu[2] = MultiFab::Dot(cumom[2],0,vel[2],0,1,0);
    }
    rhouu[0] /= (n_cells[0]+1)*n_cells[1]*n_cells[2];
    rhouu[1] /= (n_cells[1]+1)*n_cells[2]*n_cells[0];
    rhouu[2] /= (n_cells[2]+1)*n_cells[0]*n_cells[1];
    turbKE = 0.5*( rhouu[0] + rhouu[1] + rhouu[2] );

    // RMS velocity
//    StagInnerProd(vel,0,vel,0,macTemp,uu);
    {
        auto mask = vel[0].OwnerMask(geom.periodicity());
        uu[0] = MultiFab::Dot(vel[0],0,vel[0],0,1,0);
    }
    {
        auto mask = vel[1].OwnerMask(geom.periodicity());
        uu[1] = MultiFab::Dot(vel[1],0,vel[1],0,1,0);
    }
    {
        auto mask = vel[2].OwnerMask(geom.periodicity());
        uu[2] = MultiFab::Dot(vel[2],0,vel[2],0,1,0);
    }
    uu[0] /= (n_cells[0]+1)*n_cells[1]*n_cells[2];
    uu[1] /= (n_cells[1]+1)*n_cells[2]*n_cells[0];
    uu[2] /= (n_cells[2]+1)*n_cells[0]*n_cells[1];
    u_rms = sqrt((uu[0] + uu[1] + uu[2])/3.0);

    // Compute sound speed
    ComputeSoundSpeed(sound_speed,prim);
    c_speed = ComputeSpatialMean(sound_speed, 0);

    // compute gradU = [du/dx dv/dy dw/dz] at cell-centers
    ComputeCentredGradFC(vel,gradU,geom);
    ccTemp.setVal(0.0);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Add(ccTemp,gradU,d,0,1,0);
    }
    CCInnerProd(ccTemp,0,ccTemp,0,ccTempDiv,temp); // store (\sum_i du_i/dx_i)^2 MFab

    // Compute Velocity gradient moment sum
    // 2nd moment
    ccTemp.setVal(0.0);
    ccTempA.setVal(0.0);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        CCMoments(gradU,d,ccTempA,2,gradU2[d]);
        MultiFab::Add(ccTemp,ccTempA,0,0,1,0);
        gradU2[d] *= dProb; // <(du_i/dx_i)^2> each component
    }
    Real avg_mom2 = ComputeSpatialMean(ccTemp, 0); // <\sum_i (du_i/dx_i)^2>

    // 3rd moment
    ccTemp.setVal(0.0);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        CCMoments(gradU,d,ccTempA,3,gradU3[d]);
        MultiFab::Add(ccTemp,ccTempA,0,0,1,0);
        gradU3[d] *= dProb; // <(du_i/dx_i)^3> each component
    }
    Real avg_mom3 = ComputeSpatialMean(ccTemp, 0); //  <\sum_i (du_i/dx_i)^3>

    // 4th moment
    ccTemp.setVal(0.0);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        CCMoments(gradU,d,ccTempA,4,gradU4[d]);
        MultiFab::Add(ccTemp,ccTempA,0,0,1,0);
        gradU4[d] *= dProb; // <(du_i/dx_i)^4> each component
    }
    Real avg_mom4 = ComputeSpatialMean(ccTemp, 0); //  <\sum_i (du_i/dx_i)^4>

    // Taylor Microscale
    taylor_len = sqrt(3.0)*u_rms/sqrt(avg_mom2); // from Wang et al., JFM, 2012

    // Taylor Reynolds Number & Turbulent Mach number
    Real rho_avg = ComputeSpatialMean(prim, 0);
    Real eta_avg = ComputeSpatialMean(eta, 0);
    taylor_Re = rho_avg*taylor_len*u_rms/eta_avg; // from from John, Donzis, Sreenivasan, PRL 2019
    taylor_Ma = sqrt(3.0)*u_rms/c_speed; // from John, Donzis, Sreenivasan, PRL 2019

    // Skewness
    //Real skew1 = gradU3[0]/pow(gradU2[0],1.5); // <(du_1/dx_1)^3>/<(du_1/dx_1)^2>^1.5
    //Real skew2 = gradU3[1]/pow(gradU2[1],1.5); // <(du_2/dx_2)^3>/<(du_2/dx_2)^2>^1.5
    //Real skew3 = gradU3[2]/pow(gradU2[2],1.5); // <(du_3/dx_3)^3>/<(du_3/dx_3)^2>^1.5
    // <\sum_i (du_i/dx_i)^3> / (\sum_i <(du_i/dx_i)^2>^1.5)
    skew = avg_mom3/(pow(gradU2[0],1.5) + pow(gradU2[1],1.5) + pow(gradU2[2],1.5));

    // Kurtosis
    //Real kurt1 = gradU4[0]/pow(gradU2[0],2); // <(du_1/dx_1)^4>/<(du_1/dx_1)^2>^2
    //Real kurt2 = gradU4[1]/pow(gradU2[1],2); // <(du_2/dx_2)^4>/<(du_2/dx_2)^2>^2
    //Real kurt3 = gradU4[2]/pow(gradU2[2],2); // <(du_3/dx_3)^4>/<(du_3/dx_3)^2>^2
    // <\sum_i (du_i/dx_i)^4> / (\sum_i <(du_i/dx_i)^2>^2)
    kurt =  avg_mom4/(pow(gradU2[0],2) + pow(gradU2[1],2) + pow(gradU2[2],2));

    // Compute \omega (curl)
    ComputeCurlFaceToEdge(vel,curlU,geom);

    // Solenoidal dissipation: <eta \omega_i \omega_i>
    AverageCCToEdge(eta,eta_edge,0,1,SPEC_BC_COMP,geom);
    EdgeInnerProd(curlU,0,curlU,0,curlUtemp,tempvec);
//    EdgeInnerProd(curlUtemp,0,eta_edge,0,curlU,eps_s_vec);
    {
        auto mask = curlUtemp[0].OwnerMask(geom.periodicity());
        eps_s_vec[0] = MultiFab::Dot(curlUtemp[0],0,eta_edge[0],0,1,0);
    }
    {
        auto mask = curlUtemp[1].OwnerMask(geom.periodicity());
        eps_s_vec[1] = MultiFab::Dot(curlUtemp[1],0,eta_edge[1],0,1,0);
    }
    {
        auto mask = curlUtemp[2].OwnerMask(geom.periodicity());
        eps_s_vec[2] = MultiFab::Dot(curlUtemp[2],0,eta_edge[2],0,1,0);
    }
    eps_s_vec[0] /= (n_cells[0]+1)*(n_cells[1]+1)*n_cells[2];
    eps_s_vec[1] /= (n_cells[0]+1)*(n_cells[2]+1)*n_cells[1];
    eps_s_vec[2] /= (n_cells[1]+1)*(n_cells[2]+1)*n_cells[0];
    eps_s = (eps_s_vec[0] + eps_s_vec[1] + eps_s_vec[2]);

    // Dilational dissipation (4/3)*<eta (\sum_i du_i/dx_i)^2>
//    CCInnerProd(ccTempDiv,0,eta,0,ccTemp,eps_d);
    if (visc_type == 3) {
        // get eta_bulk_diss = kappa + 4/3 eta
        MultiFab::LinComb(eta_bulk_diss, 1.0, zeta, 0,
                          1.3333333333, eta, 0,
                          0, 1, 0);
        eps_d = MultiFab::Dot(eta_bulk_diss, 0, ccTempDiv, 0, 1, 0);
        eps_d *= dProb;
    }
    else {
        eps_d = MultiFab::Dot(eta, 0, ccTempDiv, 0, 1, 0);
        eps_d *= dProb*(4.0/3.0);
    }

    // Ratio of Dilational to Solenoidal dissipation
    eps_ratio = eps_d/eps_s;
    Real eps_t = eps_s + eps_d;

    // Kolmogorov scales
    kolm_s = pow((eta_avg*eta_avg*eta_avg/(rho_avg*rho_avg*eps_s)),0.25);
    kolm_d = pow((eta_avg*eta_avg*eta_avg/(rho_avg*rho_avg*eps_d)),0.25);
    kolm_t = pow((eta_avg*eta_avg*eta_avg/(rho_avg*rho_avg*eps_t)),0.25);
//    kolm_s = pow((eta_avg*eta_avg*eta_avg/eps_s),0.25);
//    kolm_d = pow((eta_avg*eta_avg*eta_avg/eps_d),0.25);
//    kolm_t = pow((eta_avg*eta_avg*eta_avg/eps_t),0.25);

}
#endif

#if defined(TURB)
void GetTurbQtyDecomp(const MultiFab& vel_decomp_in, // contains 6 components for solenoidal and dilataional velocities
                      const MultiFab& prim,
                      const amrex::Geometry& geom,
                      Real& turbKE_s, Real& turbKE_d, Real& delta_turbKE,
                      Real& u_rms_s, Real& u_rms_d, Real& delta_u_rms,
                      Real& taylor_Ma_d,
                      Real& skew_s, Real& kurt_s,
                      Real& skew_d, Real& kurt_d)

{
    BL_PROFILE_VAR("GetTurbQtyDecomp()",GetTurbQtyDecomp);

    MultiFab vel_decomp;
    vel_decomp.define(prim.boxArray(),prim.DistributionMap(),6,1); // need a ghost cell for gradients
    vel_decomp.ParallelCopy(vel_decomp_in,0,0,6);
    vel_decomp.FillBoundary(geom.periodicity());

    Vector<Real> dProb(3);
    dProb[0] = 1.0/((n_cells[0]+1)*n_cells[1]*n_cells[2]);
    dProb[1] = 1.0/((n_cells[1]+1)*n_cells[2]*n_cells[0]);
    dProb[2] = 1.0/((n_cells[2]+1)*n_cells[0]*n_cells[1]);

    // Setup temp MultiFabs
    std::array< MultiFab, AMREX_SPACEDIM > gradU;
    std::array< MultiFab, AMREX_SPACEDIM > faceTemp;
    MultiFab sound_speed;
    MultiFab ccTemp;
    AMREX_D_TERM(gradU[0].define(convert(prim.boxArray(),nodal_flag_x), prim.DistributionMap(), 6, 0);,
                 gradU[1].define(convert(prim.boxArray(),nodal_flag_y), prim.DistributionMap(), 6, 0);,
                 gradU[2].define(convert(prim.boxArray(),nodal_flag_z), prim.DistributionMap(), 6, 0););
    AMREX_D_TERM(faceTemp[0].define(convert(prim.boxArray(),nodal_flag_x), prim.DistributionMap(), 1, 0);,
                 faceTemp[1].define(convert(prim.boxArray(),nodal_flag_y), prim.DistributionMap(), 1, 0);,
                 faceTemp[2].define(convert(prim.boxArray(),nodal_flag_z), prim.DistributionMap(), 1, 0););
    sound_speed.define(prim.boxArray(),prim.DistributionMap(),1,0);
    ccTemp.define(prim.boxArray(),prim.DistributionMap(),1,0);

    // Setup temp variables
    Vector<Real> gradU2_temp(3);
    Vector<Real> gradU2_s(3);
    Vector<Real> gradU3_s(3);
    Vector<Real> gradU4_s(3);
    Vector<Real> gradU2_d(3);
    Vector<Real> gradU3_d(3);
    Vector<Real> gradU4_d(3);

    Vector<int> comps_s{0,1,2};
    Vector<int> comps_d{3,4,5};

    // turbulent kinetic energy (solenoidal)
    ccTemp.setVal(0.0);
    MultiFab::AddProduct(ccTemp,vel_decomp,0,vel_decomp,0,0,1,0); //uu
    MultiFab::AddProduct(ccTemp,vel_decomp,1,vel_decomp,1,0,1,0); //vv
    MultiFab::AddProduct(ccTemp,vel_decomp,2,vel_decomp,2,0,1,0); //ww
    u_rms_s = ComputeSpatialMean(ccTemp, 0);
    u_rms_s = sqrt(u_rms_s/3.0);
    MultiFab::Multiply(ccTemp,prim,0,0,1,0); // rho*(uu+vv+ww)
    turbKE_s = ComputeSpatialMean(ccTemp,0);
    turbKE_s = 0.5*turbKE_s;

    // turbulent kinetic energy (dilatational)
    ccTemp.setVal(0.0);
    MultiFab::AddProduct(ccTemp,vel_decomp,3,vel_decomp,3,0,1,0); //uu
    MultiFab::AddProduct(ccTemp,vel_decomp,4,vel_decomp,4,0,1,0); //vv
    MultiFab::AddProduct(ccTemp,vel_decomp,5,vel_decomp,5,0,1,0); //ww
    u_rms_d = ComputeSpatialMean(ccTemp, 0);
    u_rms_d = sqrt(u_rms_d/3.0);
    MultiFab::Multiply(ccTemp,prim,0,0,1,0); // rho*(uu+vv+ww)
    turbKE_d = ComputeSpatialMean(ccTemp,0);
    turbKE_d = 0.5*turbKE_d;

    // ratio of turbulent kinetic energies
    delta_u_rms  = u_rms_d/u_rms_s;
    delta_turbKE = turbKE_d/turbKE_s;

    // Taylor Mach (dilatational)
    ComputeSoundSpeed(sound_speed,prim);
    Real c_speed = ComputeSpatialMean(sound_speed, 0);
    taylor_Ma_d = sqrt(3.0)*u_rms_d/c_speed;

    // compute gradU = [du/dx dv/dy dw/dz] at cell-centers
    ComputeGrad(vel_decomp,gradU,0,0,6,-1,geom,0);

    // Compute Velocity gradient moment sum
    // 2nd moment (solenoidal)
    FCMoments(gradU,comps_s,faceTemp,2,gradU2_temp);
    gradU2_s[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
    gradU2_s[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
    gradU2_s[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
    // 2nd moment (dilatational)
    FCMoments(gradU,comps_d,faceTemp,2,gradU2_temp);
    gradU2_d[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
    gradU2_d[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
    gradU2_d[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));

    // Compute Velocity gradient moment sum
    // 3rd moment (solenoidal)
    FCMoments(gradU,comps_s,faceTemp,3,gradU2_temp);
    gradU3_s[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
    gradU3_s[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
    gradU3_s[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
    // 3rd moment (dilatational)
    FCMoments(gradU,comps_d,faceTemp,3,gradU2_temp);
    gradU3_d[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
    gradU3_d[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
    gradU3_d[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));

    // Compute Velocity gradient moment sum
    // 4th moment (solenoidal)
    FCMoments(gradU,comps_s,faceTemp,4,gradU2_temp);
    gradU4_s[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
    gradU4_s[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
    gradU4_s[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
    // 4th moment (dilatational)
    FCMoments(gradU,comps_d,faceTemp,4,gradU2_temp);
    gradU4_d[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
    gradU4_d[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
    gradU4_d[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));

    // Skewness
    // <\sum_i (du_i/dx_i)^3> / (\sum_i <(du_i/dx_i)^2>^1.5)
    skew_s = (gradU3_s[0] + gradU3_s[1] + gradU3_s[2])/
             (pow(gradU2_s[0],1.5) + pow(gradU2_s[1],1.5) + pow(gradU2_s[2],1.5));
    skew_d = (gradU3_d[0] + gradU3_d[1] + gradU3_d[2])/
             (pow(gradU2_d[0],1.5) + pow(gradU2_d[1],1.5) + pow(gradU2_d[2],1.5));

    // Kurtosis
    // <\sum_i (du_i/dx_i)^4> / (\sum_i <(du_i/dx_i)^2>^2)
    kurt_s = (gradU4_s[0] + gradU4_s[1] + gradU4_s[2])/
             (pow(gradU2_s[0],2.0) + pow(gradU2_s[1],2.0) + pow(gradU2_s[2],2.0));
    kurt_d = (gradU4_d[0] + gradU4_d[1] + gradU4_d[2])/
             (pow(gradU2_d[0],2.0) + pow(gradU2_d[1],2.0) + pow(gradU2_d[2],2.0));

}
#endif


void EvaluateWritePlotFileVelGrad(int step,
                                 const amrex::Real time,
                                 const amrex::Geometry& geom,
                                 const std::array<MultiFab, AMREX_SPACEDIM>& vel)
{
    BL_PROFILE_VAR("EvaluateWritePlotFileVelGrad()",EvaluateWritePlotFileVelGrad);

    // Evaluate velocity gradient components and divergence and vorticity
    MultiFab vel_grad;

    // Cell-Centered Velocity Gradient Stats (1,2,3 are directions)
    // 0: u_1,1
    // 1: u_2,2
    // 2: u_3,3
    // 3: u_1,2
    // 4: u_1,3
    // 5: u_2,3
    // 6: divergence = u_1,1 + u_2,2 + u_3,3
    // 7: vorticity = sqrt(wx + wy + wz)
    vel_grad.define(convert(vel[0].boxArray(),IntVect(AMREX_D_DECL(0,0,0))), vel[0].DistributionMap(), 8, 0);
    vel_grad.setVal(0.0);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(vel_grad,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & vgrad = vel_grad.array(mfi);

        AMREX_D_TERM(Array4<Real const> const& velx = vel[0].array(mfi);,
                     Array4<Real const> const& vely = vel[1].array(mfi);,
                     Array4<Real const> const& velz = vel[2].array(mfi););

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // u_1,1
            vgrad(i,j,k,0) = (velx(i+1,j,k) - velx(i,j,k))/dx[0];

            // u_2,2
            vgrad(i,j,k,1) = (vely(i,j+1,k) - vely(i,j,k))/dx[1];

            // u_3,3
            vgrad(i,j,k,2) = (velz(i,j,k+1) - velz(i,j,k))/dx[2];

            // divergence
            vgrad(i,j,k,6) = vgrad(i,j,k,0) + vgrad(i,j,k,1) + vgrad(i,j,k,2);

            // on edges: u_1,2 and u_2,1 and curl w1 = u_2,1 - u_1,2
            Real u12_mm = (velx(i,j,k) - velx(i,j-1,k))/dx[1];
            Real u21_mm = (vely(i,j,k) - vely(i-1,j,k))/dx[0];
            Real w1_mm  = u21_mm - u12_mm;
            Real u12_mp = (velx(i,j+1,k) - velx(i,j,k))/dx[1];
            Real u21_mp = (vely(i,j+1,k) - vely(i-1,j+1,k))/dx[0];
            Real w1_mp  = u21_mp - u12_mp;
            Real u12_pm = (velx(i+1,j,k) - velx(i+1,j-1,k))/dx[1];
            Real u21_pm = (vely(i+1,j,k) - vely(i,j,k))/dx[0];
            Real w1_pm  = u21_pm - u12_pm;
            Real u12_pp = (velx(i+1,j+1,k) - velx(i+1,j,k))/dx[1];
            Real u21_pp = (vely(i+1,j+1,k) - vely(i,j+1,k))/dx[0];
            Real w1_pp  = u21_pp - u12_pp;

            // u_1,2
            vgrad(i,j,k,3) = 0.25*(u12_mm + u12_mp + u12_pm + u12_pp);

            // on edges: u_1,3 and u_3,1 and curl w2 = u_1,3 - u_3,1
            Real u13_mm = (velx(i,j,k) - velx(i,j,k-1))/dx[2];
            Real u31_mm = (velz(i,j,k) - velz(i-1,j,k))/dx[0];
            Real w2_mm  = u13_mm - u31_mm;
            Real u13_mp = (velx(i,j,k+1) - velx(i,j,k))/dx[2];
            Real u31_mp = (velz(i,j,k+1) - velz(i-1,j,k+1))/dx[0];
            Real w2_mp  = u13_mp - u31_mp;
            Real u13_pm = (velx(i+1,j,k) - velx(i+1,j,k-1))/dx[2];
            Real u31_pm = (velz(i+1,j,k) - velz(i,j,k))/dx[0];
            Real w2_pm  = u13_pm - u31_pm;
            Real u13_pp = (velx(i+1,j,k+1) - velx(i+1,j,k))/dx[2];
            Real u31_pp = (velz(i+1,j,k+1) - velz(i,j,k+1))/dx[0];
            Real w2_pp  = u13_pp - u31_pp;

            // u_1,3
            vgrad(i,j,k,4) = 0.25*(u13_mm + u13_mp + u13_pm + u13_pp);

            // on edges: u_2,3 and u_3,2 and curl w2 = u_3,2 - u_2,3
            Real u23_mm = (vely(i,j,k) - vely(i,j,k-1))/dx[2];
            Real u32_mm = (velz(i,j,k) - velz(i,j-1,k))/dx[1];
            Real w3_mm  = u32_mm - u23_mm;
            Real u23_mp = (vely(i,j,k+1) - vely(i,j,k))/dx[2];
            Real u32_mp = (velz(i,j,k+1) - velz(i,j-1,k+1))/dx[1];
            Real w3_mp  = u32_mp - u23_mp;
            Real u23_pm = (vely(i,j+1,k) - vely(i,j+1,k-1))/dx[2];
            Real u32_pm = (velz(i,j+1,k) - velz(i,j,k))/dx[1];
            Real w3_pm  = u32_pm - u23_pm;
            Real u23_pp = (vely(i,j+1,k+1) - vely(i,j+1,k))/dx[2];
            Real u32_pp = (velz(i,j+1,k+1) - velz(i,j,k+1))/dx[1];
            Real w3_pp  = u32_pp - u23_pp;

            // u_2,3
            vgrad(i,j,k,5) = 0.25*(u23_mm + u23_mp + u23_pm + u23_pp);

            // vorticity magnitude: sqrt(w1*w1 + w2*w2 + w3*w3)
            vgrad(i,j,k,7) = sqrt(0.25*(w1_mm*w1_mm + w1_mp*w1_mp + w1_pm*w1_pm + w1_pp*w1_pp +
                                        w2_mm*w2_mm + w2_mp*w2_mp + w2_pm*w2_pm + w2_pp*w2_pp +
                                        w3_mm*w3_mm + w3_mp*w3_mp + w3_pm*w3_pm + w3_pp*w3_pp));
        });
    }

    // Write on a plotfile
    std::string plotfilename = amrex::Concatenate("vel_grad",step,9);
    amrex::Vector<std::string> varNames(8);
    varNames[0] = "ux_x";
    varNames[1] = "uy_y";
    varNames[2] = "uz_z";
    varNames[3] = "ux_y";
    varNames[4] = "ux_z";
    varNames[5] = "uy_z";
    varNames[6] = "div";
    varNames[7] = "vort";
    WriteSingleLevelPlotfile(plotfilename,vel_grad,varNames,geom,time,step);
}

#if defined(TURB)
void EvaluateWritePlotFileVelGrad(int step,
                                 const amrex::Real time,
                                 const amrex::Geometry& geom,
                                 const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                                 const amrex::MultiFab& vel_decomp_in)
{
    BL_PROFILE_VAR("EvaluateWritePlotFileVelGrad()",EvaluateWritePlotFileVelGrad);

    MultiFab output;

    // Cell-Centered Velocity Gradient Stats (1,2,3 are directions)
    // 0: ux_s
    // 1: uy_s
    // 2: uz_s
    // 3: ux_d
    // 4: uy_d
    // 5: uz_d
    // 6: umag_s
    // 7: umag_d
    // 8: divergence = u_1,1 + u_2,2 + u_3,3
    // 9: vorticity = sqrt(wx + wy + wz)
    output.define(convert(vel[0].boxArray(),IntVect(AMREX_D_DECL(0,0,0))), vel[0].DistributionMap(), 10, 0);
    output.setVal(0.0);

    MultiFab vel_decomp;
    vel_decomp.define(convert(vel[0].boxArray(),IntVect(AMREX_D_DECL(0,0,0))), vel[0].DistributionMap(), 6, 0);
    vel_decomp.ParallelCopy(vel_decomp_in,0,0,6);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(output,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<      Real>&             out   = output.array(mfi);

        const Array4<const Real>&  v_decomp         = vel_decomp.array(mfi);

        AMREX_D_TERM(Array4<Real const> const& velx = vel[0].array(mfi);,
                     Array4<Real const> const& vely = vel[1].array(mfi);,
                     Array4<Real const> const& velz = vel[2].array(mfi););

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            out(i,j,k,0) = v_decomp(i,j,k,0);
            out(i,j,k,1) = v_decomp(i,j,k,1);
            out(i,j,k,2) = v_decomp(i,j,k,2);
            out(i,j,k,3) = v_decomp(i,j,k,3);
            out(i,j,k,4) = v_decomp(i,j,k,4);
            out(i,j,k,5) = v_decomp(i,j,k,5);

            out(i,j,k,6) = std::sqrt(out(i,j,k,0)*out(i,j,k,0) + out(i,j,k,1)*out(i,j,k,1) + out(i,j,k,2)*out(i,j,k,2)); // mag solenoidal
            out(i,j,k,7) = std::sqrt(out(i,j,k,3)*out(i,j,k,3) + out(i,j,k,4)*out(i,j,k,4) + out(i,j,k,5)*out(i,j,k,5)); // mag dilatational

            // divergence
            out(i,j,k,8) = (velx(i+1,j,k) - velx(i,j,k))/dx[0] +
                           (vely(i,j+1,k) - vely(i,j,k))/dx[1] +
                           (velz(i,j,k+1) - velz(i,j,k))/dx[2] ;

            // on edges: u_1,2 and u_2,1 and curl w1 = u_2,1 - u_1,2
            Real u12_mm = (velx(i,j,k) - velx(i,j-1,k))/dx[1];
            Real u21_mm = (vely(i,j,k) - vely(i-1,j,k))/dx[0];
            Real w1_mm  = u21_mm - u12_mm;
            Real u12_mp = (velx(i,j+1,k) - velx(i,j,k))/dx[1];
            Real u21_mp = (vely(i,j+1,k) - vely(i-1,j+1,k))/dx[0];
            Real w1_mp  = u21_mp - u12_mp;
            Real u12_pm = (velx(i+1,j,k) - velx(i+1,j-1,k))/dx[1];
            Real u21_pm = (vely(i+1,j,k) - vely(i,j,k))/dx[0];
            Real w1_pm  = u21_pm - u12_pm;
            Real u12_pp = (velx(i+1,j+1,k) - velx(i+1,j,k))/dx[1];
            Real u21_pp = (vely(i+1,j+1,k) - vely(i,j+1,k))/dx[0];
            Real w1_pp  = u21_pp - u12_pp;

            // on edges: u_1,3 and u_3,1 and curl w2 = u_1,3 - u_3,1
            Real u13_mm = (velx(i,j,k) - velx(i,j,k-1))/dx[2];
            Real u31_mm = (velz(i,j,k) - velz(i-1,j,k))/dx[0];
            Real w2_mm  = u13_mm - u31_mm;
            Real u13_mp = (velx(i,j,k+1) - velx(i,j,k))/dx[2];
            Real u31_mp = (velz(i,j,k+1) - velz(i-1,j,k+1))/dx[0];
            Real w2_mp  = u13_mp - u31_mp;
            Real u13_pm = (velx(i+1,j,k) - velx(i+1,j,k-1))/dx[2];
            Real u31_pm = (velz(i+1,j,k) - velz(i,j,k))/dx[0];
            Real w2_pm  = u13_pm - u31_pm;
            Real u13_pp = (velx(i+1,j,k+1) - velx(i+1,j,k))/dx[2];
            Real u31_pp = (velz(i+1,j,k+1) - velz(i,j,k+1))/dx[0];
            Real w2_pp  = u13_pp - u31_pp;

            // on edges: u_2,3 and u_3,2 and curl w2 = u_3,2 - u_2,3
            Real u23_mm = (vely(i,j,k) - vely(i,j,k-1))/dx[2];
            Real u32_mm = (velz(i,j,k) - velz(i,j-1,k))/dx[1];
            Real w3_mm  = u32_mm - u23_mm;
            Real u23_mp = (vely(i,j,k+1) - vely(i,j,k))/dx[2];
            Real u32_mp = (velz(i,j,k+1) - velz(i,j-1,k+1))/dx[1];
            Real w3_mp  = u32_mp - u23_mp;
            Real u23_pm = (vely(i,j+1,k) - vely(i,j+1,k-1))/dx[2];
            Real u32_pm = (velz(i,j+1,k) - velz(i,j,k))/dx[1];
            Real w3_pm  = u32_pm - u23_pm;
            Real u23_pp = (vely(i,j+1,k+1) - vely(i,j+1,k))/dx[2];
            Real u32_pp = (velz(i,j+1,k+1) - velz(i,j,k+1))/dx[1];
            Real w3_pp  = u32_pp - u23_pp;

            // vorticity magnitude: sqrt(w1*w1 + w2*w2 + w3*w3)
            out(i,j,k,9) = std::sqrt(0.25*(w1_mm*w1_mm + w1_mp*w1_mp + w1_pm*w1_pm + w1_pp*w1_pp +
                                      w2_mm*w2_mm + w2_mp*w2_mp + w2_pm*w2_pm + w2_pp*w2_pp +
                                      w3_mm*w3_mm + w3_mp*w3_mp + w3_pm*w3_pm + w3_pp*w3_pp));
        });
    }

    // Write on a plotfile
    std::string plotfilename = amrex::Concatenate("vel_grad_decomp",step,9);
    amrex::Vector<std::string> varNames(10);
    varNames[0] = "ux_s";
    varNames[1] = "uy_s";
    varNames[2] = "uz_s";
    varNames[3] = "ux_d";
    varNames[4] = "uy_d";
    varNames[5] = "uz_d";
    varNames[6] = "umag_s";
    varNames[7] = "umag_d";
    varNames[8] = "div";
    varNames[9] = "vort";
    WriteSingleLevelPlotfile(plotfilename,output,varNames,geom,time,step);
}
#endif

#if defined(TURB)
void EvaluateWritePlotFileVelGradTiny(int step,
                                 const amrex::Real time,
                                 const amrex::Geometry& geom,
                                 const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                                 const amrex::MultiFab& vel_decomp_in)
{
    BL_PROFILE_VAR("EvaluateWritePlotFileVelGradTiny()",EvaluateWritePlotFileVelGradTiny);

    MultiFab output;

    // 0: vorticity wx_sifted
    // 1: vorticity wy_shifted
    // 2: vorticity wz_shifted
    // 3: vorticity wx_avg
    // 4: vorticity wy_avg
    // 5: vorticity wz_avg
    // 6: vorticity_mag_shft_then_sq = sqrt(wx + wy + wz)
    // 7: vorticity_mag_avg_then_sq = sqrt(wx + wy + wz)
    // 8: vorticity_mag_sq_then_avg = sqrt(wx + wy + wz)
    // 9: divergence =  u_1,1 + u_2,2 + u_3,3
    output.define(convert(vel[0].boxArray(),IntVect(AMREX_D_DECL(0,0,0))), vel[0].DistributionMap(), 10, 0);
    output.setVal(0.0);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(output,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<      Real>&             out   = output.array(mfi);

        AMREX_D_TERM(Array4<Real const> const& velx = vel[0].array(mfi);,
                     Array4<Real const> const& vely = vel[1].array(mfi);,
                     Array4<Real const> const& velz = vel[2].array(mfi););

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            // divergence
            out(i,j,k,9) = (velx(i+1,j,k) - velx(i,j,k))/dx[0] +
                           (vely(i,j+1,k) - vely(i,j,k))/dx[1] +
                           (velz(i,j,k+1) - velz(i,j,k))/dx[2] ;

            // on edges: u_1,2 and u_2,1 and curl w1 = u_2,1 - u_1,2
            Real u12_mm = (velx(i,j,k) - velx(i,j-1,k))/dx[1];
            Real u21_mm = (vely(i,j,k) - vely(i-1,j,k))/dx[0];
            Real w1_mm  = u21_mm - u12_mm;
            Real u12_mp = (velx(i,j+1,k) - velx(i,j,k))/dx[1];
            Real u21_mp = (vely(i,j+1,k) - vely(i-1,j+1,k))/dx[0];
            Real w1_mp  = u21_mp - u12_mp;
            Real u12_pm = (velx(i+1,j,k) - velx(i+1,j-1,k))/dx[1];
            Real u21_pm = (vely(i+1,j,k) - vely(i,j,k))/dx[0];
            Real w1_pm  = u21_pm - u12_pm;
            Real u12_pp = (velx(i+1,j+1,k) - velx(i+1,j,k))/dx[1];
            Real u21_pp = (vely(i+1,j+1,k) - vely(i,j+1,k))/dx[0];
            Real w1_pp  = u21_pp - u12_pp;
            out(i,j,k,0) = w1_mm;
            out(i,j,k,3) = 0.5*(w1_mm+w1_mp+w1_pm+w1_pp);

            // on edges: u_1,3 and u_3,1 and curl w2 = u_1,3 - u_3,1
            Real u13_mm = (velx(i,j,k) - velx(i,j,k-1))/dx[2];
            Real u31_mm = (velz(i,j,k) - velz(i-1,j,k))/dx[0];
            Real w2_mm  = u13_mm - u31_mm;
            Real u13_mp = (velx(i,j,k+1) - velx(i,j,k))/dx[2];
            Real u31_mp = (velz(i,j,k+1) - velz(i-1,j,k+1))/dx[0];
            Real w2_mp  = u13_mp - u31_mp;
            Real u13_pm = (velx(i+1,j,k) - velx(i+1,j,k-1))/dx[2];
            Real u31_pm = (velz(i+1,j,k) - velz(i,j,k))/dx[0];
            Real w2_pm  = u13_pm - u31_pm;
            Real u13_pp = (velx(i+1,j,k+1) - velx(i+1,j,k))/dx[2];
            Real u31_pp = (velz(i+1,j,k+1) - velz(i,j,k+1))/dx[0];
            Real w2_pp  = u13_pp - u31_pp;
            out(i,j,k,1) = w2_mm;
            out(i,j,k,4) = 0.5*(w2_mm+w2_mp+w2_pm+w2_pp);

            // on edges: u_2,3 and u_3,2 and curl w2 = u_3,2 - u_2,3
            Real u23_mm = (vely(i,j,k) - vely(i,j,k-1))/dx[2];
            Real u32_mm = (velz(i,j,k) - velz(i,j-1,k))/dx[1];
            Real w3_mm  = u32_mm - u23_mm;
            Real u23_mp = (vely(i,j,k+1) - vely(i,j,k))/dx[2];
            Real u32_mp = (velz(i,j,k+1) - velz(i,j-1,k+1))/dx[1];
            Real w3_mp  = u32_mp - u23_mp;
            Real u23_pm = (vely(i,j+1,k) - vely(i,j+1,k-1))/dx[2];
            Real u32_pm = (velz(i,j+1,k) - velz(i,j,k))/dx[1];
            Real w3_pm  = u32_pm - u23_pm;
            Real u23_pp = (vely(i,j+1,k+1) - vely(i,j+1,k))/dx[2];
            Real u32_pp = (velz(i,j+1,k+1) - velz(i,j,k+1))/dx[1];
            Real w3_pp  = u32_pp - u23_pp;
            out(i,j,k,2) = w3_mm;
            out(i,j,k,5) = 0.5*(w3_mm+w3_mp+w3_pm+w3_pp);

            // vorticity magnitude: sqrt(w1*w1 + w2*w2 + w3*w3)
            out(i,j,k,6) = sqrt(w1_mm*w1_mm + w2_mm*w2_mm + w3_mm*w3_mm);
            out(i,j,k,7) = sqrt(out(i,j,k,4)*out(i,j,k,4) + out(i,j,k,5)*out(i,j,k,5)
                                + out(i,j,k,6)*out(i,j,k,6));
            out(i,j,k,8) = std::sqrt(0.25*(w1_mm*w1_mm + w1_mp*w1_mp + w1_pm*w1_pm + w1_pp*w1_pp +
                                      w2_mm*w2_mm + w2_mp*w2_mp + w2_pm*w2_pm + w2_pp*w2_pp +
                                      w3_mm*w3_mm + w3_mp*w3_mp + w3_pm*w3_pm + w3_pp*w3_pp));
        });
    }

    // Write on a plotfile
    std::string plotfilename = amrex::Concatenate("vort_div",step,9);
    amrex::Vector<std::string> varNames(10);
    varNames[0] = "w1_shift";
    varNames[1] = "w2_shift";
    varNames[2] = "w3_shift";
    varNames[3] = "w1_avg";
    varNames[4] = "w2_avg";
    varNames[5] = "w3_avg";
    varNames[6] = "vort_mag_shft";
    varNames[7] = "vort_mag_shft_avg";
    varNames[8] = "vort_mag_avg";
    varNames[9] = "div";
    WriteSingleLevelPlotfile(plotfilename,output,varNames,geom,time,step);
}
#endif
