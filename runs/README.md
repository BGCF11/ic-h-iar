# Local execution directories

The research drivers use filenames relative to the current working directory. Run each experiment from a subdirectory here, using the commands in the [reproduction guide](../docs/REPRODUCING.md).

Git ignores local run contents. The input dataset and reference results remain under `data/` and `results/`. Each new execution uses a fresh directory because the drivers write fixed output filenames.
