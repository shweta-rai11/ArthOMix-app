# 16_Nomogram

`TX_MODULES` stage `id = "nomogram"`.

- **File**: `mod_nomogram.R`.
- **Input**: `results$diagnostic`/`results$featureselection` (final gene panel).
- **Main analysis**: clinical nomogram + decision curve analysis (`rms` package).
- **Output**: `results$nomogram`.
- **UI**: Transcriptomics → Sub-modules → Nomogram.
