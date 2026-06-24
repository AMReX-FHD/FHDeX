#include "surfchem_mui_functions.H"
#include "AMReX_ParmParse.H"

#define BETA 0.5

AMREX_GPU_MANAGED int surfchem_mui::nspec_mui;

void InitializeSurfChemMUINamespace()
{
    // extract inputs parameters
    ParmParse pp;

    // number of species involved in mui (via adsorption/desorption)
    pp.get("nspec_mui",nspec_mui);

    return;
}

void mui_print_MF_bottom(MultiFab& mf,int n,double factor,int step, const char *str1,const char *str2)
{
    Print() << "** mui_print_MF_bottom: " << str1 << " at step " << step << " **\n";

    for (MFIter mfi(mf,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & mf_arr = mf.array(mfi);

        // unless bx contains cells at the interface, skip
        int k = 0;
        if (k<lo.z || k>hi.z) continue;

        for (int i = lo.x; i<=hi.x; ++i)
        {
            for (int j = lo.y; j<= hi.y; ++j)
            {
                AllPrint() << str2 << " " << i << " " << j << " " << factor*mf_arr(i,j,k,n) << "\n";
            }
        }
    }
}

void mui_push(MultiFab& cu, MultiFab& prim, const amrex::Real* dx, mui::uniface2d &uniface, const int step)
// this routine pushes the following information to MUI
// - species number densities and temperature of FHD cells contacting the interface
// it assumes that the interface is perpendicular to the z-axis
// and includes cells with the smallest value of z (i.e. k=0)
{
    for (MFIter mfi(cu,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & cu_arr = cu.array(mfi);
        const Array4<Real> & prim_arr = prim.array(mfi);

        // unless bx contains cells at the interface, skip
        int k = 0;
        if (k<lo.z || k>hi.z) continue;

        for (int i = lo.x; i<=hi.x; ++i)
        {
            for (int j = lo.y; j<= hi.y; ++j)
            {
                double x = prob_lo[0]+(i+0.5)*dx[0];
                double y = prob_lo[1]+(j+0.5)*dx[1];

                std::string channel;

                for (int n = 0; n < nspec_mui; ++n)
                {
                    channel = "CH_density";
                    channel += '0'+(n+1);   // assuming nspec_mui<10

                    double dens = cu_arr(i,j,k,5+n);    // mass density
                    dens *= AVONUM/molmass[n];          // number density

                    uniface.push(channel,{x,y},dens);
                }

                channel = "CH_temp";
                double temp = prim_arr(i,j,k,4);
                uniface.push(channel,{x,y},temp);
            }
        }
    }

    return;
}

void mui_commit(mui::uniface2d &uniface, const int step)
{
    uniface.commit(step);

    return;
}

void mui_fetch(MultiFab& cu, MultiFab& prim, const amrex::Real* dx, mui::uniface2d &uniface, const int step)
// this routine fetches the following information from MUI:
// - adsoprtion and desoprtion counts of each species between time points
// it assumes that the interface is perpendicular to the z-axis
// and includes cells with the smallest value of z (i.e. k=0)
{
    mui::sampler_kmc_fhd2d<int> s({dx[0],dx[1]});
    mui::chrono_sampler_exact2d t;

    for (MFIter mfi(cu,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & cu_arr = cu.array(mfi);
        const Array4<Real> & prim_arr = prim.array(mfi);

        // unless bx contains cells at the interface, skip
        // ad-hoc fix to avoid memory leakage
        int k = 0;
        if (k<lo.z || k>hi.z)
        {
            double x = prob_lo[0]+(lo.x+0.5)*dx[0];
            double y = prob_lo[1]+(lo.y+0.5)*dx[1];

            uniface.fetch("CH_ac1",{x,y},step,s,t);

            continue;
        }

        for (int i = lo.x; i<=hi.x; ++i)
        {
            for (int j = lo.y; j<= hi.y; ++j)
            {
                double x = prob_lo[0]+(i+0.5)*dx[0];
                double y = prob_lo[1]+(j+0.5)*dx[1];

                for (int n = 0; n < nspec_mui; ++n)
                {
                    std::string channel;
                    int ac,dc;

                    channel = "CH_ac";
                    channel += '0'+(n+1);   // assuming nspec_mui<10
                    ac = uniface.fetch(channel,{x,y},step,s,t);

                    channel = "CH_dc";
                    channel += '0'+(n+1);   // assuming nspec_mui<10
                    dc = uniface.fetch(channel,{x,y},step,s,t);

                    // update

                    amrex::Real dN = ac-dc;
                    amrex::Real T_inst = prim_arr(i,j,k,4);
                    amrex::Real factor1 = molmass[n]/AVONUM/(dx[0]*dx[1]*dx[2]);
                    amrex::Real factor2 = (BETA*k_B*T_inst+(e0[n]+hcv[n]*T_inst)*molmass[n]/AVONUM)/(dx[0]*dx[1]*dx[2]);

                    cu_arr(i,j,k,0) -= factor1*dN;
                    cu_arr(i,j,k,5+n) -= factor1*dN;
                    cu_arr(i,j,k,4) -= factor2*dN;
                }
            }
        }
    }

    return;
}

void mui_fetch_Ntot(MultiFab& Ntot, const amrex::Real* dx, mui::uniface2d &uniface, const int step)
{
    mui::sampler_kmc_fhd2d<int> s({dx[0],dx[1]});
    mui::chrono_sampler_exact2d t;

    for (MFIter mfi(Ntot,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & Ntot_arr = Ntot.array(mfi);

        // unless bx contains the bottom layer, there is nothing to update
        // however, to avoid memory leakage, we do an ad-hoc fix
        int k = 0;
        if (k<lo.z || k>hi.z)
        {
            double x = prob_lo[0]+(lo.x+0.5)*dx[0];
            double y = prob_lo[1]+(lo.y+0.5)*dx[1];

            uniface.fetch("CH_one",{x,y},step,s,t);

            continue;
        }

        // if bx contains the bottom layer, need to update Ntot
        for (int i = lo.x; i<=hi.x; ++i)
        {
            for (int j = lo.y; j<= hi.y; ++j)
            {
                double x = prob_lo[0]+(i+0.5)*dx[0];
                double y = prob_lo[1]+(j+0.5)*dx[1];

                Ntot_arr(i,j,k,0) = uniface.fetch("CH_one",{x,y},step,s,t);

                AllPrint() << "Ntot(" << i << "," << j << ")= " << Ntot_arr(i,j,k,0) << "\n";
            }
        }
    }

    return;
}

void mui_fetch_surfcov(MultiFab& Ntot, MultiFab& surfcov, const amrex::Real* dx, mui::uniface2d &uniface, const int step)
{
    // reset surfcov to zero
    surfcov.setVal(0.);

    mui::sampler_kmc_fhd2d<int> s({dx[0],dx[1]});
    mui::chrono_sampler_exact2d t;

    for (MFIter mfi(surfcov,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);
        const Array4<Real> & surfcov_arr = surfcov.array(mfi);
        const Array4<Real> & Ntot_arr = Ntot.array(mfi);

        // unless bx contains the bottom layer, there is nothing to update
        // however, to avoid memory leakage, we do an ad-hoc fix
        int k = 0;
        if (k<lo.z || k>hi.z)
        {
            double x = prob_lo[0]+(lo.x+0.5)*dx[0];
            double y = prob_lo[1]+(lo.y+0.5)*dx[1];

            uniface.fetch("CH_occ1",{x,y},step,s,t);

            continue;
        }

        // if bx contains the bottom layer, need to update surfcov
        for (int i = lo.x; i<=hi.x; ++i)
        {
            for (int j = lo.y; j<= hi.y; ++j)
            {
                double x = prob_lo[0]+(i+0.5)*dx[0];
                double y = prob_lo[1]+(j+0.5)*dx[1];

                // to compute surfcov for each species
                // get the number of sites occupied by the species
                for (int n = 0; n < nspec_mui; ++n)
                {
                    std::string channel;

                    channel = "CH_occ";
                    channel += '0'+(n+1);   // assuming nspec_mui<10
                    int Nocc = uniface.fetch(channel,{x,y},step,s,t);
                    surfcov_arr(i,j,k,n) = Nocc/Ntot_arr(i,j,k,0);
                }
            }
        }
    }

    return;
}

void mui_forget(mui::uniface2d &uniface, const int step)
{
    uniface.forget(step);

    return;
}


void mui_announce_send_recv_span(mui::uniface2d &uniface,MultiFab& mf,const Real* dx)
{
    // find the lo and hi points of a square that covers all boxes assigned to the MPI process

    int lox,loy,loz,hix,hiy,hiz;

    bool isfirst = true;

    for (MFIter mfi(mf,false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Dim3 lo = lbound(bx);
        Dim3 hi = ubound(bx);

        if (isfirst)
        {
            lox = lo.x;
            loy = lo.y;
            loz = lo.z;
            hix = hi.x;
            hiy = hi.y;
            hiz = hi.z;

            isfirst = false;
        }
        else
        {
            lox = (lox<lo.x) ? lox : lo.x;
            loy = (loy<lo.y) ? loy : lo.y;
            loz = (loz<lo.z) ? loz : lo.z;
            hix = (hix>hi.x) ? hix : hi.x;
            hiy = (hiy>hi.y) ? hiy : hi.y;
            hiz = (hiz>hi.z) ? hiz : hi.z;
        }
    }

    // announce span to MUI
    // we assume that the FHD layer contacting the KMC is k = 0
    int k = 0;

    if (k>=loz && k<=hiz)
    {
        double tmp[2];

        tmp[0] = prob_lo[0] + lox*dx[0];
        tmp[1] = prob_lo[1] + loy*dx[1];
        point<double,2> span_lo(tmp);

        tmp[0] = prob_lo[0] + (hix+1)*dx[0];
        tmp[1] = prob_lo[1] + (hiy+1)*dx[1];
        point<double,2> span_hi(tmp);

        mui::geometry::box<config_2d> span(span_lo,span_hi);

        uniface.announce_send_span(0.,(double)max_step,span);
        uniface.announce_recv_span(0.,(double)max_step,span);
    }
    else
    {
        double tmp[2];

        tmp[0] = -1.;
        tmp[1] = -1.;
        point<double,2> span_lo(tmp);

        tmp[0] = -0.9;
        tmp[1] = -0.9;
        point<double,2> span_hi(tmp);

        mui::geometry::box<config_2d> span(span_lo,span_hi);

        uniface.announce_send_span(0.,(double)max_step,span);
        uniface.announce_recv_span(0.,(double)max_step,span);
    }

    return;
}
