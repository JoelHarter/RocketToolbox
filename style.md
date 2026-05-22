# RocketToolbox Style Guide

Short rules for RocketToolbox source. Goal: expressive, general-purpose scientific tools
that are a pleasure to use interactively and in simulation code.

## Hard Rules

- **Files must end with a single trailing newline**

- **New features must be tested** — every new function or type must be validated manually
  before committing, and an automated check must be added to `test/test.jl`. The test suite
  must pass for any pull request to be merged.

- **All numeric values must be in base SI units** — m/s not km/h, N not lbf, rad not deg, etc.
  This applies everywhere: function arguments, return values, struct fields, constants.
  Unit conversions belong at the boundary where data enters or leaves the library.

## Naming

- **Use Unicode symbols — do not spell out Greek letters or other math symbols.**
  Write `Ω`, `ω`, `θ`, `φ`, `μ`, `σ` — not `Omega`, `omega`, `theta`, `phi`, `mu`, `sigma`.
  Write `α`, `β`, `γ` — not `alpha`, `beta`, `gamma`.
  This applies to variable names, function names, constants, and struct fields.

- **Prefer subscript and superscript digits over regular digits in names.**
  Write `r₁`, `x₂`, `a²`, `b³` — not `r1`, `x2`, `a2`, `b3`.
  Subscript digits (`₀₁₂₃₄₅₆₇₈₉`) and superscript digits (`⁰¹²³⁴⁵⁶⁷⁸⁹`) are preferred
  wherever the number is part of the mathematical name, not an arbitrary index.

- **Unicode in general is preferred over ASCII workarounds.**
  Use math operators (`×`, `⋅`, `⊗`, `↺`) and arrows (`→`) where they improve readability.
  Julia supports this natively — use it.

## Types and Arrays

- **Use static arrays (`SVector`, `SMatrix`, `FieldVector`) when the size is fixed by the
  mathematical meaning of the object.**
  A 3D rotation vector is always length 3 — it should be an `SVector{3}`, not a `Vector`.
  A quaternion is always length 4 — it should be a `FieldVector{4}` (or `Quat`).
  Use dynamic `Vector`/`Matrix` only when size genuinely varies at runtime or is not
  known at definition time.

## Coordinate Type Names

Position types: `ECI`, `ECEF`, `Geo`, `ENU`, `NED`, `BFrame`, `Coord3D`.
Relative vector types (velocity, acceleration, force, angular rate, etc.): `ECIΔ`, `ECEFΔ`,
`ENUΔ`, `NEDΔ`, `BFrame_rel`, `Coord3D_rel`.
6D position+velocity state types: `ECI6`, `ECEF6`.
`BFrame` is the vehicle body-fixed frame. The short acronym `BF` was considered but `BFrame`
was chosen because `BF` is not universally recognized outside aerospace.

## Scope

- **Tools are for general-purpose prototyping and data analysis.**
  RocketToolbox is a scientific utility library — small, reusable building blocks for
  math, geometry, coordinate transforms, and the like. Functions should be general
  enough to be useful across many projects. If something is specific to one vehicle,
  one mission, one application, or one algorithm, it does not belong here.
