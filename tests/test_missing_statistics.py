"""Run: python3 tests/test_missing_statistics.py [workflow_directory].

Requires Python 3 and Rscript. Tests the actual scripts, including plots and
category splitting, using temporary files. No biological files are modified.
"""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
scripts = root / "00-scripts"
rr, het, aa = "0/0:20:20,0", "0/1:20:10,10", "1/1:20:0,20"
cases = {
    "all_missing": ["./.:.:.,."] * 4,
    "mono_alt": [aa] * 4,
    "no_hets": [rr, rr, aa, aa],
    "zero_ad": [rr, "0/1:0:0,0", "0/1:0:0,0", aa],
    "missing_ad": [rr, "0/1:.:.", "0/1:.:.,.", aa],
    "partial_ad": [rr, "0/1:.:10,.", "0/1", aa],
    "mixed_ad": [rr, het, "0/1:.:.", aa],
    "no_hom_depth": ["0/0:.:.", het, het, "1/1:.:."],
    "all_hets": [het] * 4,
    "canonical": [rr, het, het, aa],
    "high_depth": ["0/0:100:100,0", "0/1:.:.", "0/1:.:.", "1/1:100:0,100"],
    "low_carriers": [rr, rr, rr, "0/1:.:."],
}

def run(*args):
    result = subprocess.run([str(x) for x in args], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    return result

def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

with tempfile.TemporaryDirectory(prefix="stacks-missing-") as temporary:
    work = Path(temporary)
    # Run individually too: an entirely unavailable or non-canonical dataset
    # must still produce its output table and all three plots.
    for label, selected in [("combined", cases), *[(k, {k: v}) for k, v in cases.items()]]:
        vcf = work / (label + ".vcf")
        header = ["##fileformat=VCFv4.2",
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
                  '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4"]
        for i, (name, calls) in enumerate(selected.items(), 1):
            header.append("\t".join(["1", str(i), name, "A", "C", ".", "PASS", ".", "GT:DP:AD", *calls]))
        vcf.write_text("\n".join(header) + "\n")
        summary = work / (label + ".tsv")
        run(sys.executable, scripts / "08_extract_snp_duplication_info.py", vcf, summary)
        stats = {r["ID"]: r for r in read_table(summary)}
        run("Rscript", scripts / "09_classify_snps.R", summary)
        category_path = Path(str(summary) + ".categorized")
        cats = {r["ID"]: r["Category"] for r in read_table(category_path)}
        for name, row in stats.items():
            if name in ("all_missing", "mono_alt", "no_hets", "zero_ad", "missing_ad", "partial_ad", "high_depth", "low_carriers"):
                assert row["MedRatio"] == row["AvgRatio"] == "NA", name
            if name in ("all_missing", "mono_alt"):
                assert row["Fis"] == "NA", name
            if name == "all_missing":
                assert all(row[k] == "NA" for k in ("PropHet", "PropHomFreq", "PropHomRare", "MedCovHet", "MedCovHom", "TotCovHet"))
            if name in ("zero_ad", "missing_ad", "partial_ad", "mixed_ad"):
                assert float(row["PropHet"]) == .5 and float(row["Fis"]) == 0
                assert int(row["NumHet"]) == 2 and int(row["NumRare"]) == 3
            if name == "zero_ad":
                assert float(row["MedCovHet"]) == float(row["TotCovHet"]) == 0
            if name in ("missing_ad", "partial_ad"):
                assert row["MedCovHet"] == row["TotCovHet"] == "NA"
            if name in ("all_hets", "no_hom_depth"):
                assert row["MedCovHom"] == "NA"
            expected = {"all_missing": "insufficient_data", "mono_alt": "insufficient_data",
                        "zero_ad": "insufficient_data", "missing_ad": "insufficient_data",
                        "partial_ad": "insufficient_data", "no_hets": "mas",
                        "mixed_ad": "canonical", "canonical": "canonical",
                        "no_hom_depth": "canonical", "all_hets": "diverged",
                        "high_depth": "highcov", "low_carriers": "mas"}
            assert cats[name] == expected[name], (name, cats[name])
        for i in range(1, 4):
            assert Path(str(summary) + f"_{i}.png").stat().st_size > 0
        run(sys.executable, scripts / "10_split_vcf_in_categories.py", vcf, category_path)
        for category in set(cats.values()):
            output = vcf.with_name(vcf.stem + "." + category + ".vcf")
            records = [line.split("\t")[2] for line in output.read_text().splitlines() if not line.startswith("#")]
            assert set(records) == {name for name, cat in cats.items() if cat == category}
print("PASS: missing-statistic extraction, independent flags, plots and category splitting")
