# MTEX Agent Guidelines

See [.agents/rules/mtex.md](file:///c:/Users/phkr/MTEX_7_twinSystem/.agents/rules/mtex.md) and [CLAUDE.md](file:///c:/Users/phkr/MTEX_7_twinSystem/CLAUDE.md) for full details.

- **Vectorization**: Everything is vectorized. Assume `this` represents $N$ objects and avoid writing loops over elements.
- **Code Organization**: Subsystems have localized guidelines in their subdirectories (e.g. `geometry/CLAUDE.md`, `EBSDAnalysis/CLAUDE.md`, `interfaces/CLAUDE.md`, `plotting/CLAUDE.md`, `SO3Fun/CLAUDE.md`, `TensorAnalysis/CLAUDE.md`, `tests/CLAUDE.md`).
- **Comments**: 1 line only, terse, lowercase, describing what the next line does.
- **Tests**: Standalone `check_*.m` under `tests/`. Never run full test suites without user confirmation.
- **Domain vocabulary**: Follow [CONTEXT.md](file:///c:/Users/phkr/MTEX_7_twinSystem/CONTEXT.md).
