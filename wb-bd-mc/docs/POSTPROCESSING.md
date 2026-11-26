# Physical post-processing contracts

[Overall scope](VALIDATION.md).

**Scientific workflow validation pending.** Unit, estimator and retained-data
reconstruction checks are development evidence for the specific operations
and cases tested. The publication analysis and transport estimates require
separate scientific validation on their intended datasets and conditions.

Use `from bdmc import postprocessing as pp`. Numerical kernels require NumPy; the 2D RDF transform additionally requires SciPy. Read the function docstrings for argument/return schemas.

| Purpose | API and physical contract |
|---|---|
| Read network frames | `iter_legacy_frames` consumes synchronized `.traj/.hash/.weight` and optional labeled `.chain` frames before start/stride selection. Ragged sparse rows are valid; partial/extra EOF, malformed rows and conflicting labels are errors. Split appended runs explicitly. Auxiliary columns are not physical velocities. |
| Geometry | `legacy_geometry` reads the actual current HEUR clock convention: simple shift = `Wi_tau_C * L_gradient * time`, with no extra C; step shift = `gamma_0 * L_gradient`. `minimum_image` finds the Euclidean nearest image of rectangular/Lees–Edwards lattices, including repeated crossings. Current interacting legacy shear is 3D xy. |
| Topology | `topology_statistics` separates loops, bridge chains, distinct edges and components including isolates. Self weight is two ends per loop; off-diagonal weights are reciprocal bridge multiplicities. End conservation is `sum(ends)=2*(loops+bridges)`. Components do not establish spatial percolation. |
| Bond tensors | `bond_observables` uses unordered bridges and explicit `chains` versus `distinct_edges` measures. `RR_sum=sum(w*q*q.T)`; per-particle and per-bridge denominators are different. Empty bridge means are undefined. |
| Chain midpoint/history | `chain_midpoints` uses PBC connectors and stable endpoint IDs. It is a geometric proxy, not resolved polymer COM or diffusion. `connector_correlation` follows the same labeled chains, not unlabeled edge survival. |
| Distributions | `histogram_density` uses selected total weight and exact bin measure, retaining outside mass. `empirical_cdf` combines duplicate samples; spline derivatives are not automatically positive PDFs. `spherical_bond_density` uses exact `(b^3-a^3)/3*dmu*dphi` volume with explicit zero-direction mass. For an apolar angular distribution include +q/-q with half weights. |
| Distance streams | `distance_stream` reads ragged selected samples and optionally the known six-decimal printer grammar for glued tokens. It infers neither absent zeros nor frame boundaries. Selected-distance PDFs are not all-pair RDFs. |
| Temporal correlations | `correlation` computes linear `<x(t)y(t+lag)>` with conjugated, zero-padded FFT or direct sums. Declare centering, actual origins and all-origin versus historical fixed-window estimation. Require a uniform increasing clock; do not silently regularize historical resets/drift. |
| Stress | `stress_components` requires an explicit energy schema. Symmetric columns are `xx yy zz xy xz yz`. `stress_correlation` preserves the correct isotropic projection `<dev(A):dev(B)>/10`, equivalent to off-diagonal weights 1/5 and normal-difference weights 1/30. Connector/repulsive decomposition needs both ordered cross terms. |
| Equilibrium response | `green_kubo` requires connected correlation, declared equilibrium and consistent volume/thermal energy. It returns finite configurational `G=V*C/kBT` and running integral, not converged total viscosity. Random-force RF remains a stochastic diagnostic; singular/solvent terms need separate theory. |
| Virial | `conservative_virial` includes every repulsive pair and correct bridge multiplicity, Gaussian/modified/FENE laws and units. Current HEUR energy unit is `C*kBT`, stress unit `C*kBT/length^d`; thermal energy in recorded units is `1/C`. Stored `RF/V=sum(q outer force_on_i)/V` has the opposite sign to pressure virial. The legacy repulsive polynomial is consistent here at cutoff 1; invalid FENE extension/coincident radial forces are rejected. |
| Transport | `unwrap_positions` uses exact counters; unsheared small-step inference requires explicit `assume_small_steps=True`. Evolving shear requires compatible image counters/per-frame shifts. Sparse wrapped samples cannot identify missed winding. `mean_square_displacement` uses explicitly unwrapped coordinates, declared time origins and optional integrated affine history. `diffusion_fit` needs a justified interval and computes slope/(2d); a positive slope alone does not establish a diffusive regime or uncertainty. |
| Structure factor | `structure_factor` includes the self term; the finite direct q=0 value is N. `rdf_to_structure_factor` uses the truncated bulk transforms `1+4*pi*rho*integral(r^2*(g-1)*sinc(qr/pi))` in 3D and `1+2*pi*rho*integral(r*(g-1)*J0(qr))` in 2D. Their q=0 integral limit differs from the finite forward peak. |
| Bell diagnostic | `bell_factor` is the selected Gaussian bridge expectation `mean(exp(d*alpha^2*r*lcap))`, not an automatically established rate or MC probability. |

For positional RDF/directional/planar pair density and coordination, use `bdmc.analysis` with its explicit dimension, shell, density and finite-N/bulk normalization contracts. A count→g→coordination round trip checks internal normalization; an independent geometry/count example is needed for an independent oracle. The structural gradient-first shear convention and Euclidean nearest image agree in the half-shortest-box radial domain and can differ farther into a skew cell.

Use the declared observables for a concrete dataset and follow the [diagnostic procedure](VALIDATION.md) when results disagree.
