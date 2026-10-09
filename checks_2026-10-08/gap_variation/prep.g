# prep.g -- the "stageC.py prep" computation in GAP alone (no Sage): representatives of the conjugacy classes of
# subgroups of W(B7)/<rho>, lifted, with SmallGeneratingSet; one line per class (0-based permutation lists).
g := 7;; n := 2*g;;
rho := PermList(List([0..n-1], i -> ((i+g) mod n) + 1));;
W := Centralizer(SymmetricGroup(n), rho);;
hom := NaturalHomomorphismByNormalSubgroup(W, Subgroup(W, [rho]));;
Q := ImagesSource(hom);;
cc := ConjugacyClassesSubgroups(Q);;
out := OutputTextFile("gap_out.txt", false);;
SetPrintFormattingStatus(out, false);;
for c in cc do
  H := PreImages(hom, Representative(c));
  gens := SmallGeneratingSet(H);
  if gens = [] then gens := [rho]; fi;
  WriteLine(out, String([Size(H), List(gens, x -> List(ListPerm(x, n), i -> i - 1))]));
od;
CloseStream(out);
Print("done ", Length(cc), "\n");
QUIT;
