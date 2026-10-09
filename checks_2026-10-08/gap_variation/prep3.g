# prep3.g -- fingerprints of the intermediate objects of prep.g (W, the quotient Q = W/<rho> and its permutation
# representation, the random state), to locate the step at which runs differ.
g := 7;; n := 2*g;;
fp := function(x) return HexStringInt(Sum(List(String(x), IntChar), c -> c) * 1000003 mod 2^40 + Length(String(x))); end;;
rho := PermList(List([0..n-1], i -> ((i+g) mod n) + 1));;
W := Centralizer(SymmetricGroup(n), rho);;
Print("W gens ", fp(GeneratorsOfGroup(W)), " MT ", fp(State(GlobalMersenneTwister)), "\n");
hom := NaturalHomomorphismByNormalSubgroup(W, Subgroup(W, [rho]));;
Q := ImagesSource(hom);;
Print("Q gens ", fp(GeneratorsOfGroup(Q)), " degree ", NrMovedPoints(Q), " MT ", fp(State(GlobalMersenneTwister)), "\n");
cc := ConjugacyClassesSubgroups(Q);;
Print("classes ", Length(cc), " reps ", fp(List(cc, c -> GeneratorsOfGroup(Representative(c)))), " MT ", fp(State(GlobalMersenneTwister)), "\n");
QUIT;
