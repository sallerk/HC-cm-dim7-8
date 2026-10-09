"""Fetch, from the LMFDB nf_fields API, the totally real (r2=0) field of smallest
|disc| for each transitive group nTk (n=5,6,8). Sequential, polite (sleep between
requests). Writes raw/*.json, queries.txt, totally_real.json, other_sig.json."""
import json, os, time, subprocess, datetime

BASE = "https://www.lmfdb.org/api/nf_fields/"
NT = {5: 5, 6: 16, 8: 50}
os.makedirs("raw", exist_ok=True)
new = not os.path.exists("queries.txt")
log = open("queries.txt", "a", encoding="utf-8")
if not new:
    log.write("\n--- resumed run (fetch now via curl; urllib requests got empty/non-JSON bodies) ---\n")
if new: log.write("LMFDB API queries (nf_fields). Times are UTC at request; server 'timestamp' field also recorded.\n")
if new: log.write("Note: the API returns pages of up to 100 records; _limit is ignored; _sort=disc_abs is honored (checked with _sort=-disc_abs).\n\n")


def get(url, fname):
    p = os.path.join("raw", fname)
    if os.path.exists(p):
        return json.load(open(p, encoding="utf-8"))
    for attempt in range(30):
        t = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        try:
            out = subprocess.run(["curl", "-s", "-m", "90", "-w", "\n%{http_code}", url], capture_output=True, check=True).stdout.decode("utf-8", "replace")
            body, code = out.rsplit("\n", 1)
            if not body.lstrip().startswith("{"):
                raise ValueError(f"HTTP {code}, non-JSON body: {body[:150]!r}")
            d = json.loads(body)
            open(os.path.join("raw", fname), "w", encoding="utf-8").write(body)
            log.write(f"{t}  server_ts={d.get('timestamp')}  n_records={len(d['data'])}\n  {url}\n")
            log.flush()
            return d
        except Exception as e:
            log.write(f"{t}  ERROR attempt {attempt}: {e!r}\n  {url}\n")
            log.flush()
            time.sleep(60)
    raise RuntimeError("failed " + url)


def iso(label):
    return int(label.split(".")[-1])


res, other = {}, {}
for n, kmax in NT.items():
    for k in range(1, kmax + 1):
        g = f"{n}T{k}"
        url = (f"{BASE}?degree={n}&r2=0&galois_label={g}&_format=json"
               f"&_fields=label,coeffs,disc_abs,disc_sign,r2,degree,galois_label&_sort=disc_abs")
        d = get(url, f"{g}_r2_0.json")
        recs = d["data"]
        bad = [r for r in recs if r["galois_label"] != g or r["r2"] != 0 or r["degree"] != n]
        assert not bad, (g, bad[:2])
        if recs:
            best = min(recs, key=lambda r: (int(r["disc_abs"]), iso(r["label"])))
            assert int(best["disc_abs"]) == int(recs[0]["disc_abs"])
            res[g] = {"label": best["label"], "coeffs": best["coeffs"], "disc_abs": int(best["disc_abs"]),
                      "disc_sign": best["disc_sign"], "n_records_first_page": len(recs)}
        else:
            res[g] = None
            time.sleep(25)
            url2 = (f"{BASE}?degree={n}&galois_label={g}&_format=json"
                    f"&_fields=label,disc_abs,r2,galois_label&_sort=disc_abs")
            d2 = get(url2, f"{g}_any.json")
            sig = {}
            for r in d2["data"]:
                sig.setdefault(r["r2"], r["label"])
            other[g] = {"n_records_first_page": len(d2["data"]),
                        "smallest_label_per_r2_on_first_page": {str(a): b for a, b in sorted(sig.items())}}
        print(g, res[g]["label"] if res[g] else None, flush=True)
        time.sleep(25)

json.dump(res, open("totally_real.json", "w"), indent=1)
json.dump(other, open("other_sig.json", "w"), indent=1)
log.close()
