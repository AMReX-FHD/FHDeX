#include "INS_functions.H"
#include "common_functions.H"
#include "FhdParticleContainer.H"

//#include <AMReX.H>
//#include <AMReX_AmrParGDB.H>
//#include <AMReX_ParmParse.H>
//#include <AMReX_Particles.H>
//#include <AMReX_NeighborParticles.H>

//#include <common_functions.H>

//#include <ib_functions_F.H>

//#include <iostream>


void FhdParticleContainer::InitParticles(species* particleInfo, const Real* dxp)
{
    
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi();

    int qcount = 0;

    double cosTheta, sinTheta, cosPhi, sinPhi,sep, th;    

    int pcount = 0;
        
    int ll =0;
    //double initTemp = 0;
    //double pc = 0;

    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();       


        if(ParallelDescriptor::MyProc() == 0 && mfi.LocalTileIndex() == 0)
        {

            for(int i_spec=0; i_spec < nspecies; i_spec++)
            {

            for (int i_part=0; i_part<particleInfo[i_spec].total;i_part++)
            {
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.idata(FHD_intData::sorted) = 0;
                
                p.pos(0) = prob_lo[0] + get_uniform_func()*(prob_hi[0]-prob_lo[0]);
                p.pos(1) = prob_lo[1] + get_uniform_func()*(prob_hi[1]-prob_lo[1]);
#if (BL_SPACEDIM == 3)
                p.pos(2) = prob_lo[2] + get_uniform_func()*(prob_hi[2]-prob_lo[2]);
#endif

// SC p3m testing:  two particles near each other and both near the wall
                //p.pos(0) = (geom.ProbHi(0) - geom.ProbLo(0))*0.5 + i_spec*1.82*dxp[0]; 
                //p.pos(1) = geom.ProbLo(1) + 0.21*dxp[0];

                //p.pos(1) = geom.ProbHi(1) - 0.21*dxp[0] - i_spec*1.34*dxp[0];
                //p.pos(1) = (geom.ProbHi(1) - geom.ProbLo(1))*0.501 - i_spec*1.72*dxp[0];
                //p.pos(1) = (geom.ProbHi(1) - geom.ProbLo(1))*0.5;
//#if (BL_SPACEDIM == 3)
//                p.pos(2) = (geom.ProbHi(2) - geom.ProbLo(2))*0.5; 
//#endif



// SC P3M Test #(34A):  one particle near wall, other is near the first but not the wall
//                      make sure prob_lo(2) = 0  for this one!
//                      Make sure particle count = 1 1 
//                      Make sure bc_es_lo/hi(2) = 1 aka dirichlet 
//                int ref = 2;
//                p.pos(0) = (geom.ProbHi(0))*0.5; // + i_spec*1.82*dxp[0]; 
//                p.pos(1) = ref*2.5*dxp[0] + ref*i_spec*1.82*dxp[0]; 
//#if (BL_SPACEDIM == 3)
//                p.pos(2) = (geom.ProbHi(2) - geom.ProbLo(2))*0.5; 
//#endif

//// SC P3M Test #(34B): same configuration as above but with no walls; also place charges 
////                     reflected across y=0 line. Make sure prob_lo(2) = -1*prob_hi(2) 
////                     Make sure particle count = 2 2 
////                     Make sure bc_es_lo/hi(2) = -1  aka periodic 
//                int ref = 4;
//                p.pos(0) = (geom.ProbHi(0) - geom.ProbLo(0))*0.5; // + i_spec*1.82*dxp[0]; 
//                std::cout << "i_part = "<< i_part << " i_spec = " << i_spec << std::endl; 
//                if (i_part == 0) 
//		{
//                	p.pos(1) = 2.5*ref*dxp[0] + ref*i_spec*1.82*dxp[0]; 
//		}
//		else 
//		{
//                	p.pos(1) = (-1.)*(ref*2.5*dxp[0] + ref*i_spec*1.82*dxp[0]); 
//		}

//#if (BL_SPACEDIM == 3)
//                p.pos(2) = (geom.ProbHi(2) - geom.ProbLo(2))*0.5; 
//#endif

// SC P3M: two particles in free space -- turn walls off
//                p.pos(0) = (geom.ProbHi(0) - geom.ProbLo(0))*0.5; // + i_spec*1.82*dxp[0]; 
//                //p.pos(1) = (geom.ProbHi(1) - geom.ProbLo(1))*0.10 + i_spec*1.82*dxp[0]; 
//                p.pos(1) = (geom.ProbLo(1))+5*dxp[0] + i_spec*10*dxp[0]; 
//#if (BL_SPACEDIM == 3)
//                p.pos(2) = (geom.ProbHi(2) - geom.ProbLo(2))*0.5; //  + i_spec*2.82*dxp[0];
//#endif


//                sep = 4.05;
//                th = 3.14159/6.0;
//                //th = 0;
//                cosTheta = cos(th);
//                sinTheta = sin(th);

//                p.pos(0) = phi[0]/2.0+0.08*dx[0] - (sep/2)*cosTheta*dx[0] + ll*(sep)*cosTheta*dx[0];
//                p.pos(1) = phi[1]/2.0+0.31*dx[1] - (sep/2)*sinTheta*dx[1] + ll*(sep)*sinTheta*dx[1];
//#if (BL_SPACEDIM == 3)
//                p.pos(2) = phi[2]/2.0+0.77*dx[2];
//#endif

//                p.pos(0) = phi[0]/2.0-10*dx[0] + ll*10*dx[0];
//                p.pos(1) = phi[1]/2.0;
//#if (BL_SPACEDIM == 3)
//                p.pos(2) = phi[2]/2.0;
//#endif

                ll++;
                
                p.rdata(FHD_realData::q) = particleInfo[i_spec].q;

                //Print() << "Pos: " << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2) << ", " << p.rdata(FHD_realData::q) << "\n" ;

                //original position stored for MSD calculations
                p.rdata(FHD_realData::ox) = p.pos(0);
                p.rdata(FHD_realData::oy) = p.pos(1);
#if (BL_SPACEDIM == 3)
                p.rdata(FHD_realData::oz) = p.pos(2);
#endif

                //p.rdata(FHD_realData::vx) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                //p.rdata(FHD_realData::vy) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                //p.rdata(FHD_realData::vz) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();

                p.rdata(FHD_realData::vx) = 0;
                p.rdata(FHD_realData::vy) = 0;
                p.rdata(FHD_realData::vz) = 0;

                p.rdata(FHD_realData::ux) = 0;
                p.rdata(FHD_realData::uy) = 0;
                p.rdata(FHD_realData::uz) = 0;

                p.rdata(FHD_realData::ax) = 0;
                p.rdata(FHD_realData::ay) = 0;
                p.rdata(FHD_realData::az) = 0;

                p.rdata(FHD_realData::fx) = 0;
                p.rdata(FHD_realData::fy) = 0;
                p.rdata(FHD_realData::fz) = 0;

                p.rdata(FHD_realData::travelTime) = 0;
                p.rdata(FHD_realData::diffAv) = 0;
                p.rdata(FHD_realData::stepCount) = 0;

                p.rdata(FHD_realData::mass) = particleInfo[i_spec].m; //mass
                p.rdata(FHD_realData::R) = particleInfo[i_spec].R; //R
                p.rdata(FHD_realData::radius) = particleInfo[i_spec].d/2.0; //radius
                p.rdata(FHD_realData::accelFactor) = -6*3.14159265359*p.rdata(FHD_realData::radius)/p.rdata(FHD_realData::mass); //acceleration factor (replace with amrex c++ constant for pi...)
                p.rdata(FHD_realData::dragFactor) = 6*3.14159265359*p.rdata(FHD_realData::radius); //drag factor
                //p.rdata(FHD_realData::dragFactor) = 0; //drag factor
                //p.rdata(FHD_realData::dragFactor) = 6*3.14159265359*dx[0]*1.322; //drag factor

                p.rdata(FHD_realData::wetDiff) = particleInfo[i_spec].wetDiff;
                p.rdata(FHD_realData::dryDiff) = particleInfo[i_spec].dryDiff;
                p.rdata(FHD_realData::totalDiff) = particleInfo[i_spec].totalDiff;

                p.rdata(FHD_realData::sigma) = particleInfo[i_spec].sigma;
                p.rdata(FHD_realData::eepsilon) = particleInfo[i_spec].eepsilon;

                p.idata(FHD_intData::species) = i_spec +1;
		p.rdata(FHD_realData::potential) = 0;                 

		// set distance for which we do direct, short range coulomb force calculation
		// in p3m to be 6.5*dx_poisson_grid
                p.rdata(FHD_realData::p3m_radius) = (pkernel_es + 0.5)*dxp[0];  

                particle_tile.push_back(p);

                pcount++;
            }
           }
        }
    }

    Redistribute();
    UpdateCellVectors();
    ReBin();
    clearNeighbors();
    fillNeighbors();

}
