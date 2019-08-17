#include "INS_functions.H"
#include "common_functions.H"
#include "FhdParticleContainer.H"
#include <iostream>
#include <fstream>

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

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        //Assuming tile=box for now, i.e. no tiling.
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();       

        if(ParallelDescriptor::MyProc() == 0)
        {
                int i_spec = 0;

                  std::ifstream membraneFile("membrane.dat");



//                int memsize = (int)sqrt(particle_count[0]);
//                double div = (geom.ProbHi(0) - geom.ProbLo(0))/(double)memsize;


//                for(int ii=0; ii<memsize; ii++)
//                {
//                for(int jj=0; jj<memsize; jj++)
//                {
                  for(int i=0; i<particle_count[0]; i++)
                  {

                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.idata(IntData::sorted) = 0;

                    membraneFile >> p.pos(0);                       
                    membraneFile >> p.pos(1);
                    membraneFile >> p.pos(2);

//                    p.pos(0) = div/2 + ii*div;
//                    p.pos(1) = div/2 + jj*div;
//    #if (BL_SPACEDIM == 3)
//                    p.pos(2) = (geom.ProbHi(2) - geom.ProbLo(2))*0.5;
//    #endif
                    
                    //Print() << "ID: " << p.id() << "\n";

                    p.rdata(RealData::q) = particleInfo[i_spec].q;

                    Print() << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2) << "\n" ;

                    //original position stored for MSD calculations
                    p.rdata(RealData::ox) = p.pos(0);
                    p.rdata(RealData::oy) = p.pos(1);
    #if (BL_SPACEDIM == 3)
                    p.rdata(RealData::oz) = p.pos(2);
    #endif

                    //p.rdata(RealData::vx) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                    //p.rdata(RealData::vy) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                    //p.rdata(RealData::vz) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();

                    p.rdata(RealData::vx) = 0;
                    p.rdata(RealData::vy) = 0;
                    p.rdata(RealData::vz) = 0;

                    p.rdata(RealData::ux) = 0;
                    p.rdata(RealData::uy) = 0;
                    p.rdata(RealData::uz) = 0;

                    p.rdata(RealData::ax) = 0;
                    p.rdata(RealData::ay) = 0;
                    p.rdata(RealData::az) = 0;

                    p.rdata(RealData::fx) = 0;
                    p.rdata(RealData::fy) = 0;
                    p.rdata(RealData::fz) = 0;

                    p.rdata(RealData::travelTime) = 0;
                    p.rdata(RealData::diffAv) = 0;
                    p.rdata(RealData::stepCount) = 0;

                    p.rdata(RealData::mass) = particleInfo[i_spec].m; //mass
                    p.rdata(RealData::R) = particleInfo[i_spec].R; //R
                    p.rdata(RealData::radius) = particleInfo[i_spec].d/2.0; //radius
                    p.rdata(RealData::accelFactor) = -6*3.14159265359*p.rdata(RealData::radius)/p.rdata(RealData::mass); //acceleration factor (replace with amrex c++ constant for pi...)
                    p.rdata(RealData::dragFactor) = 6*3.14159265359*p.rdata(RealData::radius); //drag factor
                    //p.rdata(RealData::dragFactor) = 0; //drag factor
                    //p.rdata(RealData::dragFactor) = 6*3.14159265359*dx[0]*1.322; //drag factor

                    p.rdata(RealData::wetDiff) = particleInfo[i_spec].wetDiff;
                    p.rdata(RealData::dryDiff) = particleInfo[i_spec].dryDiff;
                    p.rdata(RealData::totalDiff) = particleInfo[i_spec].totalDiff;

                    p.rdata(RealData::sigma) = particleInfo[i_spec].sigma;
                    p.rdata(RealData::eepsilon) = particleInfo[i_spec].eepsilon;

                    p.idata(IntData::species) = i_spec +1;
		            p.rdata(RealData::potential) = 0;                 

		            p.rdata(RealData::p3m_radius) = 6.5*dxp[0];   

                    particle_tile.push_back(p);

                    pcount++;
                  }


                  membraneFile.close();
//                }
//                }
        }

//        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
 //       {

            for(int i_spec=1; i_spec < nspecies; i_spec++)
            {
            for (int i_part=0; i_part<particleInfo[i_spec].ppb;i_part++)
            {
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.idata(IntData::sorted) = 0;
                
                p.pos(0) = smallEnd[0]*dx[0] + get_uniform_func()*dx[0]*(bigEnd[0]-smallEnd[0]+1);
                p.pos(1) = smallEnd[1]*dx[1] + get_uniform_func()*dx[1]*(bigEnd[1]-smallEnd[1]+1);
#if (BL_SPACEDIM == 3)
                p.pos(2) = smallEnd[2]*dx[2] + get_uniform_func()*dx[2]*(bigEnd[2]-smallEnd[2]+1);
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
//               
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
                
                p.rdata(RealData::q) = particleInfo[i_spec].q;

                //Print() << "Pos: " << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2) << ", " << p.rdata(RealData::q) << "\n" ;

                //original position stored for MSD calculations
                p.rdata(RealData::ox) = p.pos(0);
                p.rdata(RealData::oy) = p.pos(1);
#if (BL_SPACEDIM == 3)
                p.rdata(RealData::oz) = p.pos(2);
#endif

                //p.rdata(RealData::vx) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                //p.rdata(RealData::vy) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();
                //p.rdata(RealData::vz) = sqrt(particleInfo.R*particleInfo.T)*get_particle_normal_func();

                p.rdata(RealData::vx) = 0;
                p.rdata(RealData::vy) = 0;
                p.rdata(RealData::vz) = 0;

                p.rdata(RealData::ux) = 0;
                p.rdata(RealData::uy) = 0;
                p.rdata(RealData::uz) = 0;

                p.rdata(RealData::ax) = 0;
                p.rdata(RealData::ay) = 0;
                p.rdata(RealData::az) = 0;

                p.rdata(RealData::fx) = 0;
                p.rdata(RealData::fy) = 0;
                p.rdata(RealData::fz) = 0;

                p.rdata(RealData::travelTime) = 0;
                p.rdata(RealData::diffAv) = 0;
                p.rdata(RealData::stepCount) = 0;

                p.rdata(RealData::mass) = particleInfo[i_spec].m; //mass
                p.rdata(RealData::R) = particleInfo[i_spec].R; //R
                p.rdata(RealData::radius) = particleInfo[i_spec].d/2.0; //radius
                p.rdata(RealData::accelFactor) = -6*3.14159265359*p.rdata(RealData::radius)/p.rdata(RealData::mass); //acceleration factor (replace with amrex c++ constant for pi...)
                p.rdata(RealData::dragFactor) = 6*3.14159265359*p.rdata(RealData::radius); //drag factor
                //p.rdata(RealData::dragFactor) = 0; //drag factor
                //p.rdata(RealData::dragFactor) = 6*3.14159265359*dx[0]*1.322; //drag factor

                p.rdata(RealData::wetDiff) = particleInfo[i_spec].wetDiff;
                p.rdata(RealData::dryDiff) = particleInfo[i_spec].dryDiff;
                p.rdata(RealData::totalDiff) = particleInfo[i_spec].totalDiff;

                p.rdata(RealData::sigma) = particleInfo[i_spec].sigma;
                p.rdata(RealData::eepsilon) = particleInfo[i_spec].eepsilon;

                p.idata(IntData::species) = i_spec +1;
		p.rdata(RealData::potential) = 0;                 

		// set distance for which we do direct, short range coulomb force calculation
		// in p3m to be 6.5*dx_poisson_grid
		p.rdata(RealData::p3m_radius) = 6.5*dxp[0];   

                particle_tile.push_back(p);

                pcount++;
            }
 //           }
        }

        
    }

//    std::cout << "pcount: " << pcount << "\n";
    UpdateCellVectors();
    Redistribute();
    ReBin();

}
