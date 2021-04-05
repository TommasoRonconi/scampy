A Python interface for Sub-halo Clustering and Abundance Matching
-----------------------------------------------------------------

ScamPy is a highly-optimized and flexible library for "painting" an observed population of cosmological objects on top of the DM-halo/subhalo hierarchy obtained from DM-only simulations.
The method used combines the classical Halo Occupation Distribution (HOD) with the sub-halo abundance matching (SHAM), the sinergy of the two processes is dubbed Sub-halo clustering and
abundance matching (SCAM).
The procedure itself is quite easy since it only requires to apply the two methods in sequence:

1. by applying the HOD scheme, the host sub-haloes are selected;
2. the SHAM algorithm associates to each sub-halo an observable property of choice.
