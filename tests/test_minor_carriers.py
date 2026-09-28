"""Run with python3 tests/test_minor_carriers.py [workflow_directory].

Requires Python 3 and Rscript. Generated files are confined to a temporary
directory. Tests use supported unphased biallelic diploid genotypes.
"""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
scripts = root / "00-scripts"
cases = [(8, 2, 0), (7, 3, 0), (8, 0, 2), (7, 0, 3), (8, 1, 1), (7, 1, 2)]
with tempfile.TemporaryDirectory(prefix="stacks-mas-") as temporary:
    work = Path(temporary)
    for swapped in (False, True):
        rows = ["##fileformat=VCFv4.2",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(f"s{i}" for i in range(10))]
        for index, (rr, het, aa) in enumerate(cases, 1):
            calls = [((0, 0), (20, 0))] * rr + [((0, 1), (10, 10))] * het + [((1, 1), (0, 20))] * aa
            fields = []
            for gt, ad in calls:
                if swapped:
                    gt, ad = tuple(1-a for a in gt), ad[::-1]
                fields.append(f"{gt[0]}/{gt[1]}:20:{ad[0]},{ad[1]}")
            alleles = ["C", "A"] if swapped else ["A", "C"]
            rows.append("\t".join(["1", str(index), f"{index}_1", *alleles,
                                    ".", "PASS", ".", "GT:DP:AD", *fields]))
        vcf = work / f"input_{swapped}.vcf"
        vcf.write_text("\n".join(rows) + "\n")
        summary = work / f"summary_{swapped}.tsv"
        subprocess.run([sys.executable, str(scripts / "08_extract_snp_duplication_info.py"),
                        str(vcf), str(summary)], check=True)
        subprocess.run(["Rscript", str(scripts / "09_classify_snps.R"), str(summary)], check=True)
        with summary.open() as handle:
            stats = list(csv.DictReader(handle, delimiter="\t"))
        with Path(str(summary) + ".categorized").open() as handle:
            categories = list(csv.DictReader(handle, delimiter="\t"))
        assert len(stats) == len(categories) == len(cases)
        for (rr, het, aa), row, category in zip(cases, stats, categories):
            expected = het + min(rr, aa)
            assert int(row["NumMinorCarriers"]) == expected
            assert int(row["NumRare"]) == het + (rr if swapped else aa)
            assert (category["Category"] == "mas") == (expected < 3)
    # Old summaries fail explicitly, rather than silently using legacy counts.
    legacy = work / "legacy.tsv"
    lines = summary.read_text().splitlines()
    legacy.write_text("\n".join(line.rsplit("\t", 1)[0] for line in lines) + "\n")
    result = subprocess.run(["Rscript", str(scripts / "09_classify_snps.R"), str(legacy)],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "Missing NumMinorCarriers" in result.stderr
print("PASS: carrier counts, allele swaps, threshold boundaries, legacy counts and old-input rejection")
