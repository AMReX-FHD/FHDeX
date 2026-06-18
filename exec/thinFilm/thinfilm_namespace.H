namespace thinfilm {

    extern AMREX_GPU_MANAGED amrex::Real thinfilm_h0;
    extern AMREX_GPU_MANAGED amrex::Real thinfilm_gamma;
    extern AMREX_GPU_MANAGED amrex::Real thinfilm_pertamp;
    extern AMREX_GPU_MANAGED amrex::Real thinfilm_hamaker;

    extern AMREX_GPU_MANAGED int thinfilm_icorr;
    extern AMREX_GPU_MANAGED int thinfilm_jcorr;

    extern AMREX_GPU_MANAGED int do_fft_diag;

}
