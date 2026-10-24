# Article submission content

Material for the article planned for this product. The product code on the default branch is unchanged.

| Item | Value |
|------|-------|
| Target venue | QIP 2027 – Quantum Information Processing (مسیر پوستر) |
| Type | conf |
| Deadline | 2026-12-04 – فقط پوستر؛ مهلت سخنرانی گذشته است (۵ اکتبر ۲۰۲۶) |
| Official page | https://qipconference.org/2027/call/ |
| Abstract rule (Conference poster abstract) | Short plain-text abstract for the poster track; state the question, approach and what will be shown on the poster. Typical length 120–250 words; keywords 3–5. |

## Title
Shot-Budget Allocation for Pauli-Sum Observables: A Simulation Study of Weighted Allocation and Classical Shadows

## Abstract (draft)
Estimating the expectation value of a Hamiltonian on a quantum computer requires many repeated measurements, and the way a fixed budget of shots is divided among the Pauli terms determines the accuracy of the estimate. We compare uniform allocation, allocation weighted by coefficient magnitude and classical shadows for Pauli-sum observables of small molecular Hamiltonians, using simulation at equal total shot budget and reporting root-mean-square error and bias against the exact value. The Hamiltonians are computed in the STO-3G basis and not measured on a device, because no public hardware measurement data set was found. In simulation, weighted allocation and classical shadows gave about 11–12% lower RMSE than uniform allocation. Confidence intervals, further molecules (for example BeH2 with a sparse method) and the grouping of commuting terms are [to be measured]. The poster will show the protocol, the comparison and the limits of a simulation-only study, and we make no claim about performance on hardware.

**Keywords:** shot budget; Pauli observables; classical shadows; measurement allocation; quantum simulation

## Before submitting
1. Replace every `[bracketed]` placeholder with values the team measured; nothing marked `[to be measured]` has a result yet.
2. Patent rule: do not disclose methods covered by this product's filed claims; a patent attorney must review before submission.
3. Add author names and affiliations; confirm the abstract limits in the venue's author guidelines (the limits above are typical for the venue family).
4. Numbers are the results recorded for this venue in the page database (owner's own runs). Check them against the paper-drafts before use; anything marked [to be measured] has not been measured.
5. The Word file for the submission system is in the platform repository (`outputs/abstracts/`).
