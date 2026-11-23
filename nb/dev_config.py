"""Paths and parameters for the Gen-2 notebooks (replaces the hard-coded Dropbox paths of the legacy cells).

    from dev_config import load_dev
    dev = load_dev()            # dict: sample, pair, the input/output files, CNVconfig, ...

Reads ../config/config_dev.yaml (or $MYCNV_NB_CONFIG). The sample is `sample:` of that file or $MYCNV_NB_SAMPLE,
the data root $MYCNV_NB_DATA, the output root $MYCNV_NB_OUT. Outputs go to <out_root>/<sample key>/ in the study
layout (cov/, data/, pre/, figures/, plot/); inputs come from the cached study folder studies/<study>/CNV/<id>/.
"""

import os

import yaml

from mycnv.config import load_config

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CONFIG = os.path.join(HERE, "..", "config", "config_dev.yaml")


def _expand(paths: dict) -> dict:
    """Resolve {name} references between the path entries (in file order)."""
    out = {}
    for key, value in paths.items():
        out[key] = str(value).format(**out) if isinstance(value, str) else value
    return out


def load_dev(config_file: str | None = None, sample: str | None = None) -> dict:
    config_file = os.path.abspath(config_file or os.environ.get("MYCNV_NB_CONFIG", DEFAULT_CONFIG))
    with open(config_file) as fh:
        raw = yaml.safe_load(fh)
    paths = dict(raw["paths"])
    paths["data_root"] = os.environ.get("MYCNV_NB_DATA", paths["data_root"])
    paths["out_root"] = os.environ.get("MYCNV_NB_OUT", paths["out_root"])
    paths = _expand(paths)

    key = sample or os.environ.get("MYCNV_NB_SAMPLE", raw.get("sample"))
    s = raw["samples"][key]
    pair = s["pair"]
    study_dir = os.path.join(paths["studies"], s["study"], "CNV", str(s["id"]))
    out_dir = os.path.join(paths["out_root"], key)
    chain_dir = out_dir if raw.get("chain_from", "run") == "run" else study_dir

    dev = dict(
        config_file=config_file,
        key=key,
        sample=pair,
        study=s["study"],
        id=str(s["id"]),
        study_dir=study_dir,
        out_dir=out_dir,
        # study inputs (Snakemake outputs of the original runs)
        raw_cov=os.path.join(study_dir, "raw", f"{pair}.cov.gz"),
        raw_snp=os.path.join(study_dir, "raw", f"{pair}.snp"),  # plain text
        study_cnv_snp=os.path.join(study_dir, "data", f"{pair}.cnv.snp.gz"),
        study_cnv=os.path.join(study_dir, "data", f"{pair}.cnv.gz"),
        # notebook outputs
        roll_cov=os.path.join(out_dir, "cov", f"{pair}.roll.cov.gz"),
        cnv_snp=os.path.join(out_dir, "data", f"{pair}.cnv.snp.gz"),
        cnv=os.path.join(out_dir, "data", f"{pair}.cnv.gz"),
        pre_dir=os.path.join(out_dir, "pre"),
        fig_path=os.path.join(out_dir, "figures"),
        plot_path=os.path.join(out_dir, "plot"),
        # inputs of notebooks 5 / visualize (outputs of 3 + 4, or the study files; see chain_from)
        chain_cnv_snp=os.path.join(chain_dir, "data", f"{pair}.cnv.snp.gz"),
        chain_cnv=os.path.join(chain_dir, "data", f"{pair}.cnv.gz"),
        # PON build (0A, 0B): outputs in <out_root>/PON/
        pon_out=os.path.join(paths["out_root"], "PON"),
        pon_build=raw.get("pon_build", {}),
        CNVconfig=load_config(config_file),
        **{k: paths[k] for k in ["data_root", "out_root", "studies", "static", "pon_dir", "bed", "gc", "genmap"]},
    )
    for d in ["out_dir", "pre_dir", "fig_path", "plot_path"]:
        os.makedirs(dev[d], exist_ok=True)
    for f in ["roll_cov", "cnv_snp"]:
        os.makedirs(os.path.dirname(dev[f]), exist_ok=True)
    return dev
