#include "debug_flags.h"
#include "charm++.h"

#ifndef ORTHO_CONFIG_H
#define ORTHO_CONFIG_H

namespace cp {
    namespace ortho {

class orthoConfig
{
    public:
        //-------------------- Vars that indicate problem size / decomposition --------------------
        /// The number of states in the simulation (the dimension of the input square matrix)
        int numStates;
        /// The block size for parallelization
        int grainSize;

        //------------------- Instance, Array, Group identification etc. --------------------------
        int instanceIndex;
        /// Callback to notify bubble owner (GSpace) that a tolerance update is needed
        CkCallback uponToleranceFailure;

        /// Is an OrthoHelper chare array available to perform the step 2 computations
        bool isStep2HelperOn;
        /// Pander to the BGL NIC and split the GEMMs in ortho
        int gemmSplit;



        void pup(PUP::er &p)
        {
            p|numStates;
            p|grainSize;
            p|instanceIndex;
            p|uponToleranceFailure;
            p|isStep2HelperOn;
            p|gemmSplit;
        }
};

    } // end ortho
} // end namespace cp

#endif // ORTHO_CONFIG_H

