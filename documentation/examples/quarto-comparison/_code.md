## Source-code comparison {#source-code-comparison}

The changed operations are the viscosity and osmotic-factor arithmetic.
Compare the input, outputs and residual tolerance before changing a solver.

[Reference source](models/reference.py) / [Candidate source](models/candidate.py)

```diff
--- models/reference.py
+++ models/candidate.py
@@ -1,9 +1,9 @@
-"""Direct polynomial reference for the synthetic documentation example."""
+"""Factored polynomial implementation for the synthetic documentation example."""
 
 
 def properties(phi):
-    """Return reduced diffusivity, relative viscosity and osmotic factor."""
+    """Return the same reduced quantities with a factored evaluation order."""
     diffusion = 1.0 + 1.2 * phi
-    viscosity = 1.0 + 2.5 * phi + 2.0 * phi ** 2
-    osmotic_factor = 1.0 + 4.0 * phi + 10.0 * phi ** 2
+    viscosity = 1.0 + phi * (2.5 + 2.0 * phi)
+    osmotic_factor = 1.0 + phi * (4.0 + 10.0 * phi)
     return diffusion, viscosity, osmotic_factor
```
