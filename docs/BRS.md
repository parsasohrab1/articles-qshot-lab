# Business Requirements Specification (BRS)

**Product:** Shot-Budget Lab (articles-qshot-lab)
*Reconstructed from the README, the requirements document, the TRL assessment and the tests of this repository. Not reviewed by the owner. A patent attorney must review before publication of method details.*

## 1. Business problem
Estimating the energy of a Pauli-sum Hamiltonian on a quantum device costs measurement shots. How a fixed shot budget is split among the terms changes the error.

## 2. Business objectives
| ID | Objective |
|----|-----------|
| BO-1 | Compare uniform, weighted, oracle, adaptive and classical-shadows allocation |
| BO-2 | Report RMSE and bias against the exact value |
| BO-3 | Check that every strategy is unbiased |

## 3. Users and stakeholders (proposed; to be confirmed by the owner)
Quantum-algorithm and quantum-chemistry researchers.

## 4. Scope
In scope: simulation for 3 to 8 qubits.
Out of scope: noisy devices, commuting-group measurements, systems above 8 qubits.

## 5. Constraints and assumptions
- All data are synthetic or public benchmark data; the product makes no clinical, operational or regulatory claim.
- Maturity is TRL 4 (see docs/TRL.md).

## 6. Success measures (to be measured on real data)
RMSE per budget on a noisy simulator or device, with readout error and grouped measurements.
