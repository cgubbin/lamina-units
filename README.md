# lamina-units

Typed physical quantities and units for the lamina crates.

`lamina-units` gives every spectral and length quantity a type that carries its
unit, so a value can't cross a crate boundary without saying what it is. It
converts between the spectral parameterisations used in optics and the canonical
spectral coordinate used by the numerical crates: vacuum angular wavenumber
`k₀ = 2π/λ₀`, in `cm⁻¹`.

```rust
use lamina_units::Wavelength;

let wavelength = Wavelength::nanometres(1550.0_f64);
let k0 = wavelength.vacuum_wavenumber().vacuum_angular_wavenumber();
let expected = 2.0 * std::f64::consts::PI / 1.55e-4; // cm⁻¹
assert!((k0.as_inverse_centimetres() - expected).abs() < 1e-6);
```

## What it provides

- **Quantities:** `Wavelength`, `Length`, `Frequency`, `AngularFrequency`,
  `VacuumWavenumber`, `VacuumAngularWavenumber`, `VacuumAngularWavenumberSquared`
  and `Angle`, each generic over its float type and tagged with its unit.
- **Units:** `LengthUnit`, `InverseLengthUnit`, `FrequencyUnit`,
  `AngularFrequencyUnit` and `AngleUnit`, with their conversion factors.
- **Spectral transforms:** `SpectralQuantity`, `SpectralCoordinate` and
  `SpectralTransform` map a caller-facing spectral quantity to `k₀` together
  with the derivatives of that map, so derivatives can be reported in the
  caller's own parameterisation.
- **Labels and axis metadata:** `UnitLabel`, `CoordinateLabel`, `AxisMetadata`
  and `AxisQuantity` name the quantity and unit of a sampled axis.
- **Constants:** the speed of light in `cm·s⁻¹`.

The optional `serde` feature derives `Serialize` and `Deserialize` for the unit
enums.

## Boundary

| lamina-units owns | lamina-units never owns |
| --- | --- |
| Typed quantities, units and their conversions | Material models or material data (`lamina-materials`) |
| The canonical spectral coordinate `k₀` in `cm⁻¹` and maps into it | Layer stacks, transfer or scattering matrices (`lamina-core`) |
| Derivatives of unit and spectral maps | In-plane coordinates that need a medium: effective index, incidence angle (`lamina-core`) |
| Labels and metadata for sampled axes | Mode searches, fields or any solver (`lamina-mode`) |

A change that adds physics (a medium, a stack, a solver) to this crate crosses
the boundary and belongs in a crate further up.

## Family

```text
lamina-units       typed quantities and units (this crate)
lamina-core        stacks, transfer and scattering matrices, fields
lamina-materials   material models implementing lamina-core's material traits
lamina-mode        modes of planar multilayers, with numerical evidence
```

The crates are versioned together until 1.0; each changelog names the
conventions a release changes.
