#include "common_functions.H"

// takes in a flattened multifab
// (i.e., a multifab with only 1 cell in 1 direction but >1 cells in the other directions)
// returns a flattened multifab that is now flattened in the AMREX_SPACEDIM-1 direction
// (z direction in 3D, y direction in 2D)
MultiFab RotateFlattenedMF(MultiFab const& mf)
{
    BoxArray const& old_ba = mf.boxArray();
    DistributionMapping const& dm = mf.DistributionMap();
    Box const& domain_box = old_ba.minimalBox();
    int short_direction;
    int short_size = domain_box.shortside(short_direction);
    if (short_size != 1) {
        Print() << "RotateFlattenedMF needs a MF with short_size==1; returning the original input MultiFab\n";
        return MultiFab(mf, amrex::make_alias, 0, mf.nComp());
    } else if (short_direction == AMREX_SPACEDIM-1) {
        return MultiFab(mf, amrex::make_alias, 0, mf.nComp());
    } else {
        IntVect old_ng = mf.nGrowVect();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(old_ng[short_direction] == 0,
                                         "Not supposed to have ghost cells in the shortest direction");
        IntVect ng;
        if (short_direction == 0) {
            ng = IntVect(AMREX_D_DECL(old_ng[1],old_ng[2],0));
        } else {
            ng = IntVect(AMREX_D_DECL(old_ng[0],old_ng[2],0));
        }
        BoxList bl = old_ba.boxList();
        for (auto& b : bl) {
            const auto lo = b.smallEnd();
            const auto hi = b.bigEnd();
            if (short_direction == 0) {
                b = Box(IntVect(AMREX_D_DECL(lo[1],lo[2],0)),
                        IntVect(AMREX_D_DECL(hi[1],hi[2],0)),
                        b.ixType());
            } else {
                b = Box(IntVect(AMREX_D_DECL(lo[0],lo[2],0)),
                        IntVect(AMREX_D_DECL(hi[0],hi[2],0)),
                        b.ixType());
            }
        }
        BoxArray new_ba(std::move(bl));
        const int ncomp = mf.nComp();
        MultiFab new_mf(new_ba, dm, ncomp, ng, MFInfo().SetAlloc(false));
        for (MFIter mfi(new_mf); mfi.isValid(); ++mfi) {
            new_mf.setFab(mfi, FArrayBox(mfi.fabbox(), ncomp, mf[mfi.index()].dataPtr()));
        }
        return new_mf;
    }
}
