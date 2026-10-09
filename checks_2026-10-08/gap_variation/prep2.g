# prep2.g -- the computation of prep.g done twice in ONE GAP process, from freshly constructed groups, with the
# global random sources reset to the same state before each run; prints fingerprints of the random states and of
# the output, so that runs in different processes and in one process can be compared.
g := 7;; n := 2*g;;
fp := function(x) return HexStringInt(Sum(List(String(x), IntChar), c -> c) * 1000003 mod 2^40 + Length(String(x))); end;;
run := function(tag)
  local rho, W, hom, Q, cc, res, c, H, gens, reps;
  Reset(GlobalMersenneTwister, 1); Reset(GlobalRandomSource, 1);
  rho := PermList(List([0..n-1], i -> ((i+g) mod n) + 1));
  W := Centralizer(SymmetricGroup(n), rho);
  hom := NaturalHomomorphismByNormalSubgroup(W, Subgroup(W, [rho]));
  Q := ImagesSource(hom);
  cc := ConjugacyClassesSubgroups(Q);
  reps := List(cc, c -> GeneratorsOfGroup(Representative(c)));
  Print(tag, " after classes: MT ", fp(State(GlobalMersenneTwister)), " RS ", fp(State(GlobalRandomSource)),
        " reps ", fp(reps), "\n");
  res := [];
  for c in cc do
    H := PreImages(hom, Representative(c));
    gens := SmallGeneratingSet(H);
    Add(res, List(gens, x -> ListPerm(x, n)));
  od;
  Print(tag, " at end: MT ", fp(State(GlobalMersenneTwister)), " RS ", fp(State(GlobalRandomSource)),
        " gens ", fp(res), "\n");
end;;
run("first");
run("second");
QUIT;
