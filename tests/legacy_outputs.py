"""Produce the legacy "before" outputs with the unmodified main code (run once; the public ones are committed).

Run in the `somvar-legacy` env (py3.9, pandas 1.5.2) against a worktree of main:
  git -C <myCNV> worktree add --detach <scratch>/myCNV-main main
  conda run -n somvar-legacy python tests/legacy_outputs.py <scratch>/myCNV-main
The legacy mawk tools run per chromosome, as somVarWES ran them (pile2CNV/PON2CNV need -c 0 there).
Public outputs (tests/data/legacy/, from tests/data/synth/): pile2CNV, TN2CNV, PON2CNV, PON combine, sample combine,
rolling and ASCAT. Restricted outputs (tests/data/restricted/legacy/, git-ignored): the same steps on the
patient-derived slices, only if tests/data/restricted/ exists.
"""

import os
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

MAIN = sys.argv[1] if len(sys.argv) > 1 else "/tmp/myCNV-main"
sys.path.insert(0, os.path.join(MAIN, "code"))

from CNV_raw import PON2CNV, TN2CNV  # noqa: E402
from combineCNV import combine_sample_CNV, filter_cov, make_PON_coverage, make_PON_snp  # noqa: E402
from convert import write_ASCAT  # noqa: E402
from rollingCov import rolling_coverage  # noqa: E402
from rollingSNP import remergeCNV, rolling_snp  # noqa: E402
from yaml import safe_load  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SYNTH = os.path.join(DATA, "synth")
OUT = os.path.join(DATA, "legacy")
RES = os.path.join(DATA, "restricted")  # patient-derived inputs/outputs: git-ignored (plan D14)
ROUT = os.path.join(RES, "legacy")
STUDY_CONFIG = os.path.join(HERE, "study_config_CNV.yaml")  # copy of somVarWES configs/general/config_CNV.yaml


def bash_shell_copy(src):
    """Copy of main's shell/ with only the shebang changed to bash.

    The legacy Python wrappers execute the tools directly, so the #!/bin/sh shebang decides. On the old cluster
    /bin/sh was bash; on Debian/Ubuntu it is dash, where every option is silently ignored (AGENTS.md trap).
    """
    dst = tempfile.mkdtemp(prefix="myCNV-main-shell-")
    for name in os.listdir(src):
        with open(os.path.join(src, name)) as fh:
            text = fh.read()
        if text.startswith("#!/bin/sh"):
            text = "#!/bin/bash" + text[len("#!/bin/sh"):]
        path = os.path.join(dst, name)
        with open(path, "w") as fh:
            fh.write(text)
        os.chmod(path, 0o755)
    return dst


SHELL = bash_shell_copy(os.path.join(MAIN, "shell"))


def sh(cmd):
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, stderr=subprocess.DEVNULL)


def pon2cnv_synth_outputs():
    """Legacy filterBed | PON2CNV on the synthetic chr7 PONmatrix (covered BED)."""
    with tempfile.TemporaryDirectory() as tmp:
        snp = os.path.join(tmp, "chr7.snp")
        sh(f"gunzip -c {SYNTH}/matrix/chr7.pon.gz | bash {SHELL}/filterBed.mawk {DATA}/haev7_covered_chr7.bed -x -c chr7 | "
           f"bash {SHELL}/PON2CNV.mawk -x -o {snp} -v 0.05 -d 50 -c 0 > {OUT}/pon2cnv_synth_chr7.cov")
        shutil.copy(snp, os.path.join(OUT, "pon2cnv_synth_chr7.snp"))


SYN_CHROMS = ["chr7", "chr21", "chrX"]


def synth_chain_outputs():
    """main's PON combine, sample combine, rolling and ASCAT on the synthetic CNV PON + sample (tests/synth.py)."""
    pon_out = os.path.join(OUT, "pon_synth_CNV")
    with tempfile.TemporaryDirectory() as pon:
        shutil.copytree(os.path.join(SYNTH, "pon", "cov"), os.path.join(pon, "cov"))
        shutil.copytree(os.path.join(SYNTH, "pon", "snp"), os.path.join(pon, "snp"))
        os.makedirs(os.path.join(pon, "CNV"))
        cfg = config(PON_path=pon)
        make_PON_coverage(config=cfg, save=True)
        make_PON_snp(config=cfg, save=True)
        shutil.copytree(os.path.join(pon, "CNV"), pon_out, dirs_exist_ok=True)
    with tempfile.TemporaryDirectory() as pon:
        shutil.copytree(pon_out, os.path.join(pon, "CNV"))
        for gc in [True, False]:
            cfg = config(cov_path=os.path.join(SYNTH, "sample"), snp_path=os.path.join(SYNTH, "sample"), PON_path=pon,
                         chrom_list=SYN_CHROMS)
            cfg["coverage"]["GCnormalize"] = gc
            cov_df, snp_df = combine_sample_CNV("syn_T-N", config=cfg)
            tag = "gc" if gc else "nogc"
            cov_df.to_csv(os.path.join(OUT, f"combine_synth_{tag}.cov.gz"), sep="\t", index=False, compression="gzip")
            snp_df.to_csv(os.path.join(OUT, f"combine_synth_{tag}.snp"), sep="\t", index=False)
    raw_cov = pd.read_csv(os.path.join(OUT, "combine_synth_gc.cov.gz"), sep="\t")
    raw_snp = pd.read_csv(os.path.join(OUT, "combine_synth_gc.snp"), sep="\t")
    for tag, chroms in [("chr7", ["chr7"]), ("all", SYN_CHROMS)]:
        cfg = config()
        cov_df = filter_cov(raw_cov.loc[raw_cov["Chr"].isin(chroms)], config=cfg)
        cov_df = rolling_coverage(cov_df, config=cfg)
        snp_df = rolling_snp(raw_snp.loc[raw_snp["Chr"].isin(chroms)], cov_df, config=cfg)
        snp_df.to_csv(os.path.join(OUT, f"rolling_synth_{tag}.cnv.snp.gz"), sep="\t", index=False, compression="gzip")
        remergeCNV(snp_df, cov_df).to_csv(os.path.join(OUT, f"rolling_synth_{tag}.cnv.gz"), sep="\t", index=False,
                                          compression="gzip")
    cnv_df = pd.read_csv(os.path.join(OUT, "rolling_synth_chr7.cnv.snp.gz"), sep="\t", compression="gzip")
    os.makedirs(os.path.join(OUT, "ascat_synth_chr7"), exist_ok=True)
    write_ASCAT(cnv_df, sample="syn_T-N", outpath=os.path.join(OUT, "ascat_synth_chr7"))


def mawk_outputs():
    """Legacy filterBed | pile2CNV per chromosome on the synthetic multi-chromosome pileup."""
    bed = os.path.join(SYNTH, "multi.bed")
    covs, snps = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for chrom in ["chr7", "chr21", "chrX"]:
            snp = os.path.join(tmp, f"{chrom}.snp")
            cov = os.path.join(tmp, f"{chrom}.cov")
            sh(f"bash {SHELL}/filterBed.mawk {bed} -x -c {chrom} < {SYNTH}/multi.cpile | "
               f"bash {SHELL}/pile2CNV.mawk -x -o {snp} -v 0.2 -V 0.8 -d 15 -c 0 > {cov}")
            covs.append(pd.read_csv(cov, sep="\t", dtype=str))
            snps.append(pd.read_csv(snp, sep="\t", dtype=str))
    pd.concat(covs).to_csv(os.path.join(OUT, "pile2cnv_multi.cov"), sep="\t", index=False)
    pd.concat(snps).to_csv(os.path.join(OUT, "pile2cnv_multi.snp"), sep="\t", index=False)


def pon2cnv_restricted_outputs():
    # PON2CNV on the chr7 PONmatrix head with the padded BED (how the stored HAEv7 PON cov/ was made)
    with tempfile.TemporaryDirectory() as tmp:
        snp = os.path.join(tmp, "chr7.snp")
        sh(f"gunzip -c {RES}/chr7.pon.gz | bash {SHELL}/filterBed.mawk {DATA}/haev7_padded_chr7.bed -x -c chr7 | "
           f"bash {SHELL}/PON2CNV.mawk -x -o {snp} -v 0.05 -d 50 -c 0 > {ROUT}/pon2cnv_chr7.cov")
        shutil.copy(snp, os.path.join(ROUT, "pon2cnv_chr7.snp"))


def config(**extra):
    with open(STUDY_CONFIG) as fh:
        cfg = safe_load(fh)["CNV"]
    cfg.update(extra)
    return cfg


def tn2cnv_outputs():
    cfg = config(mawk_path=SHELL, bed_file=os.path.join(DATA, "haev7_covered_chr7.bed"),
                 gc_split_path=os.path.join(DATA, "tracks", "gc"), genmap_split_path=os.path.join(DATA, "tracks", "genmap"))
    snp = os.path.join(OUT, "tn2cnv_chr7.snp")
    cov_df, _ = TN2CNV(clean_TN_pileup_file=os.path.join(SYNTH, "chr7.cpile"), chrom="chr7", SNP_output=snp, config=cfg)
    cov_df.to_csv(os.path.join(OUT, "tn2cnv_chr7.cov.gz"), sep="\t", index=False, compression="gzip")


def pon2cnv_py_outputs():
    with tempfile.TemporaryDirectory() as pon:
        os.makedirs(os.path.join(pon, "matrix"))
        os.makedirs(os.path.join(pon, "snp"))
        shutil.copy(os.path.join(RES, "chr7.pon.gz"), os.path.join(pon, "matrix", "chr7.pon.gz"))
        cfg = config(mawk_path=SHELL, bed_file=os.path.join(DATA, "haev7_padded_chr7.bed"), PON_path=pon,
                     gc_split_path=os.path.join(DATA, "tracks", "gc"),
                     genmap_split_path=os.path.join(DATA, "tracks", "genmap"))
        cov_df, _ = PON2CNV(chrom="chr7", config=cfg)
        cov_df.to_csv(os.path.join(ROUT, "pon2cnv_py_chr7.cov.gz"), sep="\t", index=False, compression="gzip")
        shutil.copy(os.path.join(pon, "snp", "chr7.snp.gz"), os.path.join(ROUT, "pon2cnv_py_chr7.snp.gz"))


def pon_combine_outputs():
    with tempfile.TemporaryDirectory() as pon:
        shutil.copytree(os.path.join(RES, "pon_mini", "cov"), os.path.join(pon, "cov"))
        shutil.copytree(os.path.join(RES, "pon_mini", "snp"), os.path.join(pon, "snp"))
        os.makedirs(os.path.join(pon, "CNV"))
        cfg = config(PON_path=pon)
        make_PON_coverage(config=cfg, save=True)
        make_PON_snp(config=cfg, save=True)
        shutil.copytree(os.path.join(pon, "CNV"), os.path.join(ROUT, "pon_mini_CNV"), dirs_exist_ok=True)


def combine_outputs():
    for gc in [True, False]:
        cfg = config(cov_path=os.path.join(RES, "sample_mini"), snp_path=os.path.join(RES, "sample_mini"),
                     PON_path=os.path.join(ROUT, "pon_mini_legacy"), chrom_list=["chr7", "chr21", "chrX"])
        cfg["coverage"]["GCnormalize"] = gc
        cov_df, snp_df = combine_sample_CNV("mini_T-N", config=cfg)
        tag = "gc" if gc else "nogc"
        cov_df.to_csv(os.path.join(ROUT, f"combine_{tag}.cov.gz"), sep="\t", index=False, compression="gzip")
        snp_df.to_csv(os.path.join(ROUT, f"combine_{tag}.snp"), sep="\t", index=False)


def rolling_outputs():
    for tag in ["chr7", "chr7_21"]:
        cfg = config()
        cov_df = pd.read_csv(os.path.join(RES, "nhl1", f"raw_{tag}.cov.gz"), sep="\t")
        snp_df = pd.read_csv(os.path.join(RES, "nhl1", f"raw_{tag}.snp"), sep="\t")
        cov_df = filter_cov(cov_df, config=cfg)
        cov_df = rolling_coverage(cov_df, config=cfg)
        snp_df = rolling_snp(snp_df, cov_df, config=cfg)
        snp_df.to_csv(os.path.join(ROUT, f"rolling_{tag}.cnv.snp.gz"), sep="\t", index=False, compression="gzip")
        remergeCNV(snp_df, cov_df).to_csv(os.path.join(ROUT, f"rolling_{tag}.cnv.gz"), sep="\t", index=False,
                                          compression="gzip")
    cnv_df = pd.read_csv(os.path.join(ROUT, "rolling_chr7.cnv.snp.gz"), sep="\t", compression="gzip")
    os.makedirs(os.path.join(ROUT, "ascat_chr7"), exist_ok=True)
    write_ASCAT(cnv_df, sample="NHL1_A-B", outpath=os.path.join(ROUT, "ascat_chr7"))


def main():
    os.makedirs(OUT, exist_ok=True)
    mawk_outputs()
    tn2cnv_outputs()
    pon2cnv_synth_outputs()
    synth_chain_outputs()
    if not os.path.isdir(os.path.join(RES, "pon_mini")):
        print("public legacy outputs written to", OUT, "(no restricted fixtures)")
        return
    os.makedirs(ROUT, exist_ok=True)
    pon2cnv_restricted_outputs()
    pon2cnv_py_outputs()
    pon_combine_outputs()
    # the combine step needs a PON CNV/ dir: use the legacy one just made
    shutil.copytree(os.path.join(ROUT, "pon_mini_CNV"), os.path.join(ROUT, "pon_mini_legacy", "CNV"), dirs_exist_ok=True)
    combine_outputs()
    shutil.rmtree(os.path.join(ROUT, "pon_mini_legacy"))
    rolling_outputs()
    print("legacy outputs written to", OUT, "and", ROUT)


if __name__ == "__main__":
    main()
