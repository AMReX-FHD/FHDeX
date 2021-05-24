#include "DsmcParticleContainer.H"

//#include "particle_functions_K.H"
#include "paramplane_functions_K.H"
#include <math.h>

// Likely excessive
// Relates lower diagonal matrix indices to those of 1D array
int getSpeciesIndex(int species1, int species2){
	if(species1<species2){
		return species1+nspecies*(species2-1)
	} else {
		return species2+nspecies*(species1-1)
	}
}

// Use to approximate pair correlation fuunc with radial distr. func
Real g0_Ma_Ahmadi(int species1, int species2){
	const double C1 = 2.5, C2 = 4.5904
	const double C3 = 4.515439, CPOW = 0.67802
   // CHANGE: Divide by collision cell volume
   // Should be stored in the boxes of BoxArray
   // Don't think it has been set up quite yet
	Real phi1 = properties[species1].partVol
	Real phi2 = properties[speices2].partVol
	if(*species1==*species2) {
		int indx = getSpecIndx(int i, int j)
		// ADD: Solid fraction by combination of species in each collision cell
		// double phi = (PartoCellVol[i])*np[i+NSpecies*(i-1)] // auto-type casted
		double numer = 1 + C1*PoCVol + C2*pow(PoCVol,2) + C3*pow(PoCVol,3)
		numer = numer * 4.0*PoCVol
		// Declare phi
		double phiRatio = phi/phi_max
		double denom = pow(1.0 - pow(phiRatio,3),CPOW)
		return std::min(1+numer/denom,chi_max)
	} else {
	 // add Lebowitz RDF here
	}
	
}

// Same as FhdParticleContainer.H?
FhdParticleContainer::FhdParticleContainer(const Geometry & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int ncells)
    : NeighborParticleContainer<FHD_realData::count, FHD_intData::count> (geom, dmap, ba, ncells)
{
    BL_PROFILE_VAR("FhdParticleContainer()",FhdParticleContainer);

    realParticles = 0; // Do we need to keep track of this?
    simParticles = 0;
	 // Pi - replace with existing if avail
	 Real pi_usr = 4.0*atan(1.0)

    for(int i=0;i<nspecies;i++) {
        
        properties[i].mass = mass[i];
        properties[i].radius = diameter[i]/2.0;
        properties[i].partVol = pow(diameter[i],3)*pi_usr*particle_neff/6;
        // Probably useful to separate by species but not sure if Neff should be the same for each
        // needs to be checked
        properties[i].Neff = particle_neff; // From DSMC, this will be set to 1 for electolyte calcs
        // Might be unused for granular
        properties[i].R = k_B/properties[i].mass; //used a lot in kinetic stats cals, bu not otherwise necessary for electrolytes
        // Should be labeled as granular temperature for clarity
        properties[i].T = T_init[i];

        if (particle_count[i] >= 0) {

            properties[i].total = particle_count[i];
            properties[i].n0 = particle_neff*properties[i].total/domainVol;
            
            Print() << "Species " << i << " count " << properties[i].total << "\n";
            Print() << "Species " << i << " n0 " << properties[i].total << "\n";
        }
        else {
            properties[i].total = (int)amrex::Math::ceil(particle_n0[i]*domainVol/particle_neff);
            properties[i].n0 = properties[i].total/domainVol;

            Print() << "Species " << i << " count " << properties[i].total << "\n";
            Print() << "Species " << i << " n0 adjusted to " << properties[i].n0 << "\n";
        }



        Print() << "Species " << i << "\n";
        Print() << "Mass " << properties[i].mass << "\n";
        Print() << "Radius " << properties[i].radius << "\n";

        realParticles = realParticles + properties[i].total*particle_neff;
        simParticles = simParticles + properties[i].total;

    }
    
    int indx
    for(int i=0;i<nspecies;i++) {
    	for(int j=0;j<nspecies,j++){
    		indx = getSpeciesIndex(i,j);
    		propBetweenSpecies[indx] = pow(properties[i].radius + properties[j].radius,2)
    		propBetweenSpecies[indx] = propBetweenSpecies[indx]*pi_usr // replace with pi
    	}
    }

    totalCollisionCells = n_cells[0]*n_cells[1]*n_cells[2];
    domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);

    Print() << "Total real particles: " << realParticles << "\n";
    Print() << "Total sim particles: " << simParticles << "\n";

    Print() << "Collision cells: " << totalCollisionCells << "\n";
    Print() << "Sim particles per cell: " << simParticles/totalCollisionCells << "\n";


}

