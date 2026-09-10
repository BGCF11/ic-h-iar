# Numerical result tables

| File | Content |
| --- | --- |
| `Resultados_Termografia_SantaGenebra_FaseB.csv` | Pixelwise parameter estimates, train/test counts, RMSEs and numerical termination flags |
| `resultados_hiar_mc.csv` | Aggregate Monte Carlo statistics by scenario, sample size and quaternion component |

The simulation design has 4 scenarios × 3 sample sizes × 4 components, giving 48 aggregate rows when all combinations complete. The Phase B worker returns one row per retained pixel; row order follows asynchronous completion.

`Convergencia` records the optimizer success flag. The map script retains rows with `Convergencia=False` when the parameter norm is available. CSV headers and numerical precision are preserved as part of the tabular interface.

The [data dictionary](../../docs/DATA_DICTIONARY.md) defines every output column. [Local execution](../../docs/REPRODUCING.md) writes new outputs under `runs/`.
