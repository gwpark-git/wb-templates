# mBLA C++ template

Read `README.md`, `docs/SCIENTIFIC_GUIDE.org`, `docs/VALIDATION.org`,
`docs/FIXED_POINT_METHOD.org` and `docs/LEVEL_2.org` before scientific or numerical edits. Read
`LICENSE`, `NOTICE` and `THIRD_PARTY_NOTICES.org` before distribution.
Read `docs/PYTHON_CAKE.org` before changing the optional Python search.

Explain each affected function and non-obvious numerical step in the source.
Record the meaningful units, state, branch or update order alongside the code.

Preserve equations, SI scales, pressure convention, transport-law variants,
fixed-point/cake rules and frozen fixtures. A stopping message is not a
conservation or mesh-accuracy result. Keep historical researcher-confirmed
Python/C++ non-cake agreement separate from new comparisons. The public Python
solver has no cake extension. Never regenerate fixtures during testing.

Build/run outside this source. Use <=2 build workers and a destination-specific
output budget; no original researcher paths are needed at runtime. Default
GCC/Clang builds explicitly disable floating-point contraction; a new target
requires its own nonlinear/cake sensitivity check. Record failed comparisons.
