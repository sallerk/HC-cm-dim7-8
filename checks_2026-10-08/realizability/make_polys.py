import json
res = json.load(open("totally_real.json"))
rows = [f'["{g}", "{v["label"]}", {v["coeffs"]}, {v["disc_abs"]}]' for g, v in res.items() if v]
open("polys.gp", "w").write("L = [" + ", ".join(rows) + "];\n")  # one line: GP ends statements at newlines
print(len(rows), "polys;", "none for:", [g for g, v in res.items() if v is None])
