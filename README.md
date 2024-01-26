## C3-C4 FBA model
This is the code for a flux balance model of C3-C4 intermediates, which was originally based on a published genome-scale constraint-based model of C3-C4 photosynthesis (Dal'Molin et al., 2010; Mallmann et al., 2014). This model was parameterized using the physiological parameters collected in Heckmann et al. (2013). The updated model was then utilized to simulate changes in metabolic flux when the Fd-GOGAT flux is limited.

## Dependency
- R 3.0.3
- sybil library
- sybilSBML library

## How to run
The simulation script of C3-C4 intermidiate flux balance analysis with a limited Fd-GOGAT flux:

`simulation Simulations/FBA_GOGAT-limited/C2-C4_simulation.R`

## Cite
- QIMING TANG, YUHUI HUANG, XIAOXIANG NI, MING-JU AMY LYU, GENYUN CHEN, ROWAN SAGE, XIN-GUANG ZHU. Increased α-ketoglutarate links the C3-C4 intermediate state to C4 photosynthesis in the genus Flaveria. BioRxiv
- HECKMANN, D., SCHULZE, S., DENTON, A., GOWIK, U., WESTHOFF, P., WEBER, A. P. M. & LERCHER, M. J. 2013. Predicting C-4 Photosynthesis Evolution: Modular, Individually Adaptive Steps on a Mount Fuji Fitness Landscape. Cell, 153, 1579-1588.
- DAL'MOLIN, C. G. D., QUEK, L. E., PALFREYMAN, R. W., BRUMBLEY, S. M. & NIELSEN, L. K. 2010. C4GEM, a Genome-Scale Metabolic Model to Study C-4 Plant Metabolism. Plant Physiology, 154, 1871-1885.
- MALLMANN, J., HECKMANN, D., BRAUTIGAM, A., LERCHER, M. J., WEBER, A. P. M., WESTHOFF, P. & GOWIK, U. 2014. The role of photorespiration during the evolution of C-4 photosynthesis in the genus Flaveria. Elife, 3.
