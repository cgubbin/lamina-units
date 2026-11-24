# lamina-units

Typed physical quantities and units for the lamina crates.

`lamina-units` gives every spectral and length quantity a type that carries its
unit, so every value says what it is. It
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

## Scope

lamina-units provides typed spectral and length quantities with their units, exact conversions between them,
and maps from any caller-facing spectral parameterisation to the canonical coordinate `k₀` in `cm⁻¹`, together
with the derivatives of those maps. It also names sampled axes with labels and metadata.

## Family

```text
lamina-units       typed quantities and units (this crate)
lamina-core        stacks, transfer and scattering matrices, fields
lamina-materials   material models implementing lamina-core's material traits
lamina-mode        modes of planar multilayers, with numerical evidence
```

The crates are versioned together until 1.0; each changelog names the
conventions a release changes.

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or
  <https://www.apache.org/licenses/LICENSE-2.0>)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or <https://opensource.org/licenses/MIT>)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
