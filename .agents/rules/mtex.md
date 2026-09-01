# MTEX Rules for AI Assistant

MTEX is an open-source MATLAB toolbox for crystallographic texture analysis — EBSD, pole figures, ODF reconstruction, grain boundaries, and crystal geometry. No GUI; everything is a function/class API.

## Core Architectural Principles

1. **Everything is Vectorized**:
   - One object instance holds an array of many entities — `vector3d` is a point cloud, `EBSD` a whole scan, `grain2d` every grain in a map — with its properties as arrays in lockstep.
   - Methods operate elementwise across the array via operator overloading and broadcasting.
   - Always assume `this` represents $N$ objects; avoid writing for-loops over elements when vectorized methods exist.

2. **Code Layout & Subsystems**:
   - Class-per-folder: `@ClassName/` holds a class's methods as separate `.m` files (e.g. `geometry/@vector3d/angle.m`).
   - `geometry/` — vectors (`vector3d`), orientations (`quaternion` → `rotation` → `orientation`), Miller indices (`Miller`), crystal/specimen symmetries (`crystalSymmetry`, `specimenSymmetry`), and reference frames (`@referenceFrame`).
   - `EBSDAnalysis/` — EBSD data (`EBSD`), grains (`grain2d`, `grain3d`), grain boundaries (`grainBoundary`), and reconstruction (`parentGrainReconstructor`).
   - `PoleFigureAnalysis/` — pole figure data (`PoleFigure`) and ODF-from-pole-figure solvers.
   - `SO3Fun/`, `S2Fun/`, `S1Fun/` — function spaces on SO(3), the sphere $S^2$, and the circle $S^1$. Harmonic, RBF, triangulated, and homochoric representations share common abstract interfaces.
   - `TensorAnalysis/` — elastic/plastic tensor calculations (`tensor`, `stiffnessTensor`, etc.).
   - `plotting/` — figure and plot code; all multi-axis layout, colorbars, and annotations go through `mtexFigure`.
   - `interfaces/` — file format import/export (`loadData`, `loadEBSD_*`, `exportEBSD_*`).
   - `doc/` — published example scripts and documentation tutorials.
   - `mex/` — compiled MEX binaries.
   - `obsolete/`, `old/`, `compatibility/` — deprecated wrappers that warn and forward. Do not use as references for current patterns.

3. **Subsystem Detail Files**:
   Detailed subsystem rules and traps are documented in:
   - `geometry/CLAUDE.md`
   - `EBSDAnalysis/CLAUDE.md`
   - `interfaces/CLAUDE.md`
   - `plotting/CLAUDE.md`
   - `SO3Fun/CLAUDE.md`, `S2Fun/CLAUDE.md`, `S1Fun/CLAUDE.md`
   - `TensorAnalysis/CLAUDE.md`
   - `tests/CLAUDE.md`

## MATLAB Execution Guidelines

- Run MATLAB commands from the repository root so `startup.m` picks up `startup_mtex` and the path.
- Preferred batch syntax: `matlab -batch "your_command_here"`
- MTEX requires the **Statistics and Machine Learning Toolbox** (specifically `knnsearch`).
- **Tests**: `tests/` holds standalone `check_*.m` functions. Never run full test suites (`runTests`, `check_mtex`) unprompted as they are resource-intensive. Ask before running large suites.

## Commenting & Code Style

- **One line comments**: Describe what the next lines do. Match existing style: lowercase, terse, no trailing full stop, placed above or trailing the statement.
- If code is self-explanatory, omit comments.
- Do **not** write verbose historical summaries or multi-line justification blocks in the source code. Background explanations belong in commit messages, issues, or Architecture Decision Records (`docs/adr/`).
- Line endings are **LF** (`*.m text eol=lf`).

## Key Domain Concepts & Vocabulary (from `CONTEXT.md`)

- **Property (`ebsd.prop`)**: Per-pixel data (one value per measurement point, e.g. `MAD`, `BC`, `mis2mean`). Resized in lockstep when subsetting (`ebsd(ind)`).
- **Option (`ebsd.opt`)**: Scan-level/file-level data that is not per-pixel and does not resize on subsetting.
- **Header (`ebsd.opt.header`)**: Scan-level metadata captured from vendor file preambles/headers.
- **Grain**: A phase-homogeneous, spatially connected region of EBSD pixels. A phase change between adjacent pixels is always a grain boundary.
- **notIndexed**: Degenerate phase value for unindexed diffraction patterns. Can form notIndexed grains.
- **Grain boundary**: Atomic boundary edge between neighboring pixels of different grains. Stored in walk order along chains.
- **Chain**: Maximal run of grain boundary segments from one junction to the next.
- **Junction**: Vertex where the number of meeting boundary segments $\neq 2$.
- **Triple point**: A junction where exactly 3 segments meet separating 3 distinct real grains.
- **Reference Frame**: The `@referenceFrame` carries axes and plotting conventions (ADR 0003).

## Architectural Decision Records (ADRs)

Key architectural decisions are recorded under `docs/adr/`:
- `0001-ebsd-opt-header.md`: Scan metadata convention in `ebsd.opt.header`.
- `0003-reference-frame-vs-symmetry.md`: Separation of coordinate systems and crystal symmetries.
