#ifndef _MFUTIL_H
#define _MFUTIL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_IArrayBox.H>


using namespace amrex;


/********************************************************************************
 *                                                                              *
 * MFUTIL: the magical collection of MuliFab Utilities that AMReX doesn't       *
 *         currently have. These utils perfom common MultiFab operations (like  *
 *         duplicating and converting), and might be moved to AMReX at some     *
 *         time in the future.                                                  *
 *******************************************************************************/

namespace IBMFUtil {

    template<typename T> // T can be either MultiFab or iMultiFab
    std::unique_ptr<T> duplicate(const BoxArray & ba, const DistributionMapping & dm,
                                 const T & mf_in, const Periodicity & periodicity) {

        int ng = mf_in.nGrow();
        int nc = mf_in.nComp();

        std::unique_ptr<T> mf_out(new T(ba, dm, nc, ng));

        mf_out->copy(mf_in, 0, 0, nc, ng, ng);
        mf_out->FillBoundary(periodicity); // Might not be necessary ..................... TODO: AUDIT

        return mf_out;
    }



    template<typename T_src, typename T_dest>
    std::unique_ptr<T_dest> convert(const T_src & mf_in) {

        const BoxArray & ba            = mf_in.boxArray();
        const DistributionMapping & dm = mf_in.DistributionMap();
              int ncomp                = mf_in.nComp();
              int ngrow                = mf_in.nGrow();

        std::unique_ptr<T_dest> mf_out(new T_dest(ba, dm, ncomp, ngrow));

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(mf_in, true); mfi.isValid(); ++ mfi) {
            const IArrayBox & in_tile  =      mf_in[mfi];
                  FArrayBox & out_tile = (* mf_out)[mfi];

            for(BoxIterator bit(mfi.growntilebox()); bit.ok(); ++ bit) {

                for(int i = 0; i < ncomp; i ++)
                    out_tile(bit(), i) = in_tile(bit(), i);

            }
        }

        return mf_out;
    }
};


#endif
