**English** | [日本語](README.jp.md)

# RocketToolbox 🚀🛠️🧰
A Julia toolbox of standard utilities for aerospace calculation and analysis

_In development since 2026 Jan 22_

## Overview

### 📐 Constants Library
* Unit conversions
* Physical constants
* Earth and Astronomical parameters

### 🌏 Functions for transforming between coordinate frames
* Geodetic (Geo)
* Earth-Centered Earth-Fixed (ECEF)
* Earth-Centered Inertial (ECI)

### 🔄 Quaternions and Rotations
* Quaternion object definition
* Constructors for 2D, 3D principal axis, and 3D general rotation matrices

### 👽 And much more!

## Setup

### Install Packages
Julia:
```
using Pkg
Pkg.add("StaticArrays")
```

### Use in your code
Julia:
```julia
# include.jl is the single entry point — include it once to load everything.
include("[path to repo]/include.jl")
```

All functions, types, and constants are loaded into the current scope.
No module qualifiers needed: use `Earth.μ`, `Quat(...)`, `rotmat(...)`, etc. directly.