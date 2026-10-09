\\ Independent check of LMFDB fields: irreducibility, total reality (polsturm),
\\ Galois group via polgalois (galdata present), field discriminant via nfdisc.
\\ polys.gp defines L = [[group, label, coeffs_ascending, disc_abs], ...]
default(new_galois_format, 1);
read("polys.gp");
nbad = 0;
{
for (i = 1, #L,
  my(g = L[i][1], lab = L[i][2], f = Polrev(L[i][3]), n = poldegree(f));
  my(irr = polisirreducible(f), nr = polsturm(f), G = polgalois(f), D = nfdisc(f));
  my(gt = Str(n, "T", G[3]));
  my(ok = irr && nr == n && gt == g && D == L[i][4]);
  if (!ok, nbad++);
  print(g, "\t", lab, "\tirred=", irr, "\treal_roots=", nr, "/", n,
        "\tpolgalois=", gt, " ", G[4], " order=", G[1],
        "\tnfdisc=", D, " match_disc=", D == L[i][4], "\t", if (ok, "OK", "MISMATCH"));
);
}
print("total=", #L, " mismatches=", nbad);
