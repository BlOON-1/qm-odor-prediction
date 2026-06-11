# QM Descriptor Definitions

## 1. Purpose of QM descriptors

QM descriptors are included to capture interaction-relevant physicochemical variation that may not be fully represented by purely topological or 2D cheminformatic features. Within this repository, they function as molecular-level descriptors for prediction and comparison tasks rather than as direct measurements of olfactory receptor mechanism.

## 2. QM calculation workflow

Planned QM descriptors are derived from a workflow based on PubChem 3D structures, manual inspection in GaussView 6.0.16, Gaussian 16W geometry optimization and frequency analysis using B3LYP/6-311+G(2d,p) with D3(BJ), and downstream parsing with Multiwfn 3.8 dev. The repository will document descriptor generation outputs and software requirements, but external software binaries are not redistributed here.

## 3. Descriptor categories

- Electronic Energetics
- Electrostatic Landscape
- Polarizability
- Aromaticity

## 4. Descriptor list

| Descriptor | Category | Brief definition |
| --- | --- | --- |
| `HOMO_energy` | Electronic Energetics | Energy of the highest occupied molecular orbital. |
| `LUMO_energy` | Electronic Energetics | Energy of the lowest unoccupied molecular orbital. |
| `HOMO-LUMO_gap` | Electronic Energetics | Difference between LUMO and HOMO energies. |
| `First_vertical_IP` | Electronic Energetics | First vertical ionization potential estimate. |
| `First_vertical_EA` | Electronic Energetics | First vertical electron affinity estimate. |
| `Mulliken_electronegativity` | Electronic Energetics | Electronegativity proxy derived from frontier-orbital quantities. |
| `Hardness` | Electronic Energetics | Chemical hardness proxy derived from frontier-orbital quantities. |
| `Electrophilicity_index` | Electronic Energetics | Electrophilicity proxy derived from frontier-orbital quantities. |
| `Nucleophilicity_index` | Electronic Energetics | Nucleophilicity proxy used for comparative molecular characterization. |
| `Maximal_value` | Electrostatic Landscape | Maximum mapped surface quantity from the electrostatic descriptor surface. |
| `Minimal_value` | Electrostatic Landscape | Minimum mapped surface quantity from the electrostatic descriptor surface. |
| `Overall_variance` | Electrostatic Landscape | Variance of mapped surface values. |
| `Overall_average_value` | Electrostatic Landscape | Average mapped surface value. |
| `Overall_surface_area` | Electrostatic Landscape | Total mapped surface area. |
| `Nonpolar_surface_area` | Electrostatic Landscape | Surface area assigned to nonpolar regions under the chosen partition rule. |
| `Polar_surface_area` | Electrostatic Landscape | Surface area assigned to polar regions under the chosen partition rule. |
| `Molecular_polarity_index` | Electrostatic Landscape | Composite polarity measure summarizing electrostatic surface variation. |
| `Magnitude_of_dipole_moment` | Electrostatic Landscape | Magnitude of the molecular dipole moment vector. |
| `HOMA` | Aromaticity | Harmonic oscillator model of aromaticity. |
| `Bird` | Aromaticity | Bird aromaticity index. |
| `AV1245` | Aromaticity | Aromaticity-related multicenter index as reported by the chosen parsing workflow. |
| `FLU` | Aromaticity | Aromatic fluctuation index. |
| `ring_atoms_1based` | Aromaticity | Ring-atom indexing field associated with aromaticity calculations. |
| `n_rings` | Aromaticity | Number of rings evaluated in the aromaticity summary. |
| `Isotropic_average_polarizability` | Polarizability | Isotropic average polarizability. |
| `Polarizability_eigenvalue_max` | Polarizability | Largest principal polarizability eigenvalue. |
| `Polarizability_eigenvalue_avg` | Polarizability | Average principal polarizability eigenvalue. |
| `Polarizability_eigenvalue_min` | Polarizability | Smallest principal polarizability eigenvalue. |
| `Polarizability_anisotropy_def1` | Polarizability | Polarizability anisotropy under the selected definition. |

## 5. Unit and scaling considerations

Descriptor units follow the conventions of the originating QM software and parsing workflow. Release tables should preserve raw descriptor values whenever possible, while downstream machine-learning pipelines should document any transformations applied. All continuous descriptors should be standardized using training-set parameters only, with no information leakage from validation or test partitions.

## 6. Known limitations

Descriptor values depend on conformational initialization, optimization settings, and parser definitions. Some descriptors are software-specific proxies rather than directly observable physical quantities. Aromaticity summaries can also depend on ring definition and implementation details. Users should therefore interpret these variables as reproducible molecular descriptors within the documented workflow, not as standalone mechanistic proof.
