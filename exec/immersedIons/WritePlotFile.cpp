#include "INS_functions.H"
#include "common_functions.H"

#include "AMReX_PlotFileUtil.H"

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const amrex::Geometry cgeom,
                   const amrex::Geometry egeom,
                   const MultiFab& particleInstant,
                   const MultiFab& particleMeans,
                   FhdParticleContainer& particles,
                   const MultiFab& charge,
                   const MultiFab& chargeM,
                   const MultiFab& potential,
                   const MultiFab& potentialM,
                   const std::array< MultiFab, AMREX_SPACEDIM >& efield) 
{
    BL_PROFILE_VAR("WritePlotFile()",WritePlotFile);

    std::string cplotfilename = Concatenate("cplt",step,9);
    std::string eplotfilename = Concatenate("eplt",step,9);
    std::string pplotfilename = Concatenate("parplt",step,9);

    BoxArray cba = particleInstant.boxArray();
    DistributionMapping cdmap = particleInstant.DistributionMap();

    BoxArray eba = charge.boxArray();
    DistributionMapping edmap = charge.DistributionMap();

     int cnPlot = 16+2*nspecies;

    // charge, chargeM
    // pot, potM
    // Ex, Ey, Ez
    int enPlot = 6+AMREX_SPACEDIM;

    int fnPlot = nspecies*AMREX_SPACEDIM;

    MultiFab cplotfile(cba, cdmap, cnPlot, 0);

    MultiFab eplotfile(eba, edmap, enPlot, 0);

    Vector<std::string> cvarNames(cnPlot);
    Vector<std::string> evarNames(enPlot);

    Vector<std::string> fvarNames(fnPlot);

    amrex::MultiFab::Copy(eplotfile,charge    ,0,0,1,0);
    amrex::MultiFab::Copy(eplotfile,chargeM   ,0,1,1,0);
    amrex::MultiFab::Copy(eplotfile,potential ,0,2,1,0);
    amrex::MultiFab::Copy(eplotfile,potentialM,0,3,1,0);
    amrex::MultiFab::Copy(eplotfile,efield[0] ,0,4,1,0);
    amrex::MultiFab::Copy(eplotfile,efield[1] ,0,5,1,0);
    amrex::MultiFab::Copy(eplotfile,efield[2] ,0,6,1,0);

    amrex::MultiFab::Copy(cplotfile,particleInstant,0,0 ,8+nspecies,0);
    amrex::MultiFab::Copy(cplotfile,particleMeans  ,0,8+nspecies,8+nspecies,0);


    for (int l=0; l<nspecies; ++l)
    {
        for (int d=0; d<AMREX_SPACEDIM; ++d)
        {
            fvarNames[l*AMREX_SPACEDIM + d] = Concatenate("spec ", l, 0) + Concatenate(" dim ", d, 0);

          //  Print() << fvarNames[l*AMREX_SPACEDIM + d] << "\n";
        }
    }

    evarNames[0] = "chargeInstant";
    evarNames[1] = "chargeMean";
    evarNames[2] = "potentialInstant";
    evarNames[3] = "potentialMean";
    evarNames[4] = "ExInstant";
    evarNames[5] = "EyInstant";

#if (AMREX_SPACEDIM==3)
    evarNames[6] = "EzInstant";
#endif

    cvarNames[0] = "membersInstant";
    cvarNames[1] = "densityInstant";
    cvarNames[2] = "velxInstant";
    cvarNames[3] = "velyInstant";
    cvarNames[4] = "velzInstant";
    cvarNames[5] = "ixInstant";
    cvarNames[6] = "iyInstant";
    cvarNames[7] = "izInstant";

    int ccount = 8;

    for(int i=0;i<nspecies;i++)
    {
        std::string specname = Concatenate("chargeDensityInstantSpecies",i);
        cvarNames[ccount+i]= specname;
    }

    cvarNames[8+nspecies] = "membersMean";
    cvarNames[9+nspecies] = "densityMean";
    cvarNames[10+nspecies] = "velxMean";
    cvarNames[11+nspecies] = "velyMean";
    cvarNames[12+nspecies] = "velzMean";
    cvarNames[13+nspecies] = "ixMean";
    cvarNames[14+nspecies] = "iyMean";
    cvarNames[15+nspecies] = "izMean";

    ccount = 16+nspecies;

    for(int i=0;i<nspecies;i++)
    {
        std::string specname = Concatenate("chargeDensityMeanSpecies",i);
        cvarNames[ccount+i]= specname;
    }

    // timer
    Real t1 = ParallelDescriptor::second();
    
    WriteSingleLevelPlotfile(cplotfilename,cplotfile,cvarNames,cgeom,time,step);
    
    Real t2 = ParallelDescriptor::second() - t1;
    ParallelDescriptor::ReduceRealMax(t2);
    amrex::Print() << "Time spent writing plotfile cplt " << t2 << std::endl;
    
    t1 = ParallelDescriptor::second();
    
    WriteSingleLevelPlotfile(eplotfilename,eplotfile,evarNames,egeom,time,step);
    
    t2 = ParallelDescriptor::second() - t1;
    ParallelDescriptor::ReduceRealMax(t2);
    amrex::Print() << "Time spent writing plotfile eplt " << t2 << std::endl;

    //particles.Checkpoint(pplotfilename, "particle0");

    // particle ASCII
    if(plot_ascii == 1)
    {
//        MultiFab explotout(eba, edmap, 1, 0);
//        MultiFab eyplotout(eba, edmap, 1, 0);
//        MultiFab ezplotout(eba, edmap, 1, 0);

//        amrex::MultiFab::Copy(explotout,efield[0],0,0,1,0);
//        amrex::MultiFab::Copy(eyplotout,efield[1],0,0,1,0);
//        amrex::MultiFab::Copy(ezplotout,efield[2],0,0,1,0);

//        std::string asciiName1 = Concatenate("asciiCharge",step,9);
//        std::string asciiName2 = Concatenate("asciiPotential",step,9);
//        std::string asciiName3 = Concatenate("asciiEx",step,9);
//        std::string asciiName4 = Concatenate("asciiEy",step,9);
//        std::string asciiName5 = Concatenate("asciiEz",step,9);
        //std::string asciiName6 = Concatenate("asciiIx",step,9);
        //std::string asciiName7 = Concatenate("asciiIy",step,9);
        //std::string asciiName8 = Concatenate("asciiIz",step,9);

        

        //outputMFAscii(charge, asciiName1);
        //outputMFAscii(potential, asciiName2);
        //outputMFAscii(efield[0], asciiName3);
        //outputMFAscii(efield[1], asciiName4);
        //outputMFAscii(efield[2], asciiName5);

        //MultiFab ix(cba, cdmap, 1, 0);
       // MultiFab iy(cba, cdmap, 1, 0);
       // MultiFab iz(cba, cdmap, 1, 0);

       // amrex::MultiFab::Copy(ix,particleMeans,0,0,14,0);

        std::string asciiName = Concatenate("ascii_means",step,9);
        outputMFAscii(particleMeans, asciiName);

        std::string asciiName1 = Concatenate("ascii_charge_mean",step,9);
        outputMFAscii(chargeM, asciiName1);

        std::string asciiName2 = Concatenate("ascii_potential_mean",step,9);
        outputMFAscii(potentialM, asciiName2);

        std::string asciiPName = Concatenate("asciiParticles",step,9);
        particles.WriteParticlesAscii(asciiPName);
    }

    // particle in cplt file
    Vector<std::string> real_comp_names;
    Vector<std::string>  int_comp_names;

    real_comp_names.push_back("radius");
    real_comp_names.push_back("velx");
    real_comp_names.push_back("vely");
    real_comp_names.push_back("velz");
    real_comp_names.push_back("forcex");
    real_comp_names.push_back("forcey");
    real_comp_names.push_back("forcez");
    real_comp_names.push_back("pred_posx");
    real_comp_names.push_back("pred_posy");
    real_comp_names.push_back("pred_posz");
    real_comp_names.push_back("pred_velx");
    real_comp_names.push_back("pred_vely");
    real_comp_names.push_back("pred_velz");
    real_comp_names.push_back("pred_forcex");
    real_comp_names.push_back("pred_forcey");
    real_comp_names.push_back("pred_forcez");
    real_comp_names.push_back("vx");
    real_comp_names.push_back("vy");
    real_comp_names.push_back("vz");
    real_comp_names.push_back("fx");
    real_comp_names.push_back("fy");
    real_comp_names.push_back("fz");
    real_comp_names.push_back("ux");
    real_comp_names.push_back("uy");
    real_comp_names.push_back("uz");
    real_comp_names.push_back("mass");
    real_comp_names.push_back("R");
    real_comp_names.push_back("q");
    real_comp_names.push_back("accelFactor");
    real_comp_names.push_back("dragFactor");
    real_comp_names.push_back("ox");
    real_comp_names.push_back("oy");
    real_comp_names.push_back("oz");
    real_comp_names.push_back("ax");
    real_comp_names.push_back("ay");
    real_comp_names.push_back("az");
    real_comp_names.push_back("travelTime");
    real_comp_names.push_back("diffAv");
    real_comp_names.push_back("stepCount");
    real_comp_names.push_back("multi");
    real_comp_names.push_back("dryDiff");
    real_comp_names.push_back("wetDiff");
    real_comp_names.push_back("totalDiff");
    real_comp_names.push_back("sigma");
    real_comp_names.push_back("eepsilon");
    real_comp_names.push_back("potential");
    real_comp_names.push_back("p3m_radius");
    real_comp_names.push_back("spring");
    int_comp_names.push_back("sorted");
    int_comp_names.push_back("i");
    int_comp_names.push_back("j");
    int_comp_names.push_back("k");
    int_comp_names.push_back("species");
    int_comp_names.push_back("visible");
    int_comp_names.push_back("pinned");
//    int_comp_names.push_back("count");

    Vector<int> write_real_comp = {
        0, // radius
        1, // velx
        1, // vely
        1, // velz
        0, // forcex
        0, // forcey
        0, // forcez
        0, // pred_posx
        0, // pred_posy
        0, // pred_posz
        0, // pred_velx
        0, // pred_vely
        0, // pred_velz
        0, // pred_forcex
        0, // pred_forcey
        0, // pred_forcez
        0, // vx
        0, // vy
        0, // vz
        0, // fx
        0, // fy
        0, // fz
        0, // ux
        0, // uy
        0, // uz
        0, // mass
        0, // R
        1, // q
        0, // accelFactor
        0, // dragFactor
        0, // ox
        0, // oy
        0, // oz
        0, // ax
        0, // ay
        0, // az
        0, // travelTime
        0, // diffAv
        0, // stepCount
        0, // multi
        0, // dryDiff
        0, // wetDiff
        0, // totalDiff
        0, // sigma
        0, // eepsilon
        0, // potential
        0,  // p3m_radius
        0  // spring
    };

    Vector<int> write_int_comp = {
        0, // sorted
        0, // i
        0, // j
        0, // k
        1,  // species
        1,  // visible
        1  // pinned
    };

    t1 = ParallelDescriptor::second();
    
    particles.WritePlotFile(cplotfilename, "particles",
                            write_real_comp, write_int_comp, real_comp_names, int_comp_names);
    
    t2 = ParallelDescriptor::second() - t1;
    ParallelDescriptor::ReduceRealMax(t2);
    amrex::Print() << "Time spent writing particle plotfile " << t2 << std::endl;

}
