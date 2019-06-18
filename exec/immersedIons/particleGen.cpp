#include "INS_functions.H"
#include "common_functions.H"
#include "FhdParticleContainer.H"

void FhdParticleContainer::InitParticles(species* particleInfo)
{
    
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi();

    int qcount = 0;

    double cosTheta, sinTheta, cosPhi, sinPhi;    

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


//        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
 //       {

            for(int i_spec=0; i_spec < nspecies; i_spec++)
            {
            for (int i_part=0; i_part<particleInfo[i_spec].ppb;i_part++)
            {
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                p.idata(IntData::sorted) = 0;
                
//                p.pos(0) = smallEnd[0]*dx[0] + get_uniform_func()*dx[0]*(bigEnd[0]-smallEnd[0]+1);
//                p.pos(1) = smallEnd[1]*dx[1] + get_uniform_func()*dx[1]*(bigEnd[1]-smallEnd[1]+1);
#if (BL_SPACEDIM == 3)
//                p.pos(2) = smallEnd[2]*dx[2] + get_uniform_func()*dx[2]*(bigEnd[2]-smallEnd[2]+1);
#endif

                p.pos(0) = phi[0]/2.0 - 0.75*dx[0] + ll*(1.5*dx[0]);
                p.pos(1) = phi[1]/2.0;
#if (BL_SPACEDIM == 3)
                p.pos(2) = phi[2]/2.0;
#endif
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
                //p.rdata(RealData::dragFactor) = 6*3.14159265359*dx[0]*1.322; //drag factor

                p.rdata(RealData::wetDiff) = particleInfo[i_spec].wetDiff;
                p.rdata(RealData::dryDiff) = particleInfo[i_spec].dryDiff;
                p.rdata(RealData::totalDiff) = particleInfo[i_spec].totalDiff;

                p.rdata(RealData::sigma) = particleInfo[i_spec].sigma;
                p.rdata(RealData::eepsilon) = particleInfo[i_spec].eepsilon;

                p.idata(IntData::species) = i_spec +1;
		p.rdata(RealData::potential) = 0;                 
		// SPC: temporary hack--set distance for which we do direct coulomb force calculation
		//      to be same as that of the SR leonard jones
		p.rdata(RealData::coulombRadiusFactor) = particleInfo[i_spec].sigma;                 
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
