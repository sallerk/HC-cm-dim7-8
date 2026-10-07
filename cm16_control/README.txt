cm16_control -- brute-force control computations for CM types of degree 16
============================================================================

Plain Python (grp16.py also uses numpy), written separately from the enumeration code of this repository.
Sections A, B and C of ../cm16_control_data.txt summarize the outputs of these scripts.

  cyc.py        characters of (Z/m)^x; the rank dim MT(A) of a CM type S of Q(zeta_m), computed as 1 + the number
                of odd characters chi with sum_{s in S} chi(s) != 0; the stabilizer of S. Prints Lenstra's two
                examples on Q(zeta32), a variant, and S = {1..8} for m = 17                    -> cyc_out.txt
  enum_cyc.py   all 256 CM types of Q(zeta_m), phi(m) = 16 (m = 17, 32, 40, 48, 60), up to (Z/m)^x: primitive or
                not (stabilizer), rank (section A)                                             -> enum_cyc_out.txt
  fermat.py     Fermat-curve factors A_[a,b,c] of degree m, classes up to scaling and permutation: stabilizer
                and rank (section B)                                                         -> fermat_out.txt
  grp16.py      Galois CM fields of degree 16 with groups C8 x| C2 (D8, SD16, M16, C8 x C2), regular
                representation: numbers of CM types by left stabilizer, right stabilizer (reflex subgroup) and
                rank of the group-ring matrix (section C)                                    -> grp16_out.txt

Run from this folder, for example:  python enum_cyc.py > enum_cyc_out.txt   (under 1 s each).
The outputs here are from a run on 2026-10-07 (Python 3.12.9, numpy 2.1.3); those of the original runs
(2026-10-05) were not kept. The counts, ranks and examples in sections A, B and C of ../cm16_control_data.txt
agree with these outputs; the 16T labels and the remarks on the literature in that file are not computed here.
Section D of ../cm16_control_data.txt (Yanai's example) is not computed by these scripts; compare controls_lit.py.
