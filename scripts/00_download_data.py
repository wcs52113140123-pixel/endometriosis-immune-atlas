"""
00_download_data.py — download all single-cell datasets from GEO into data/raw/.

Strategy: each GEO series ships a <GSE>_RAW.tar of supplementary files. We download
it and extract per-sample matrices. Two datasets (GSE213216, GSE179640) wrap each
sample as a nested cellranger tarball (…/outs/filtered_feature_bc_matrix.h5); for
those we pull only the filtered .h5 per sample. Others expose flat per-GSM
{barcodes,features,matrix}.mtx.gz or GSM*.h5.

Per-sample tissue labels come from GEO series metadata; a resolved map is provided
in data/interim/ (add_gsm_tissue.json / manifests) as produced during the study.
"""
import os, tarfile, urllib.request, glob, gzip, shutil

RAW = "data/raw"
DATASETS = {
    # accession : short role note
    "GSE247695": "paired eutopic vs ectopic endometriosis (benign)",
    "GSE224334": "OCCC snRNA (malignant)",
    "GSE213216": "endometriosis / ovarian endometrioma / eutopic / ovary / control (benign) [nested tarballs]",
    "GSE179640": "endometrium / peritoneum / ovary / organoids (Tan 2022, benign) [nested tarballs]",
    "GSE203191": "menstrual endometrium, endo phenotypes (benign)",
    "GSE214411": "minimal/mild endometriosis endometrium (benign)",
    "GSE224333": "OCCC (malignant, independent of GSE224334)",
    "GSE291389": "normal endometrium -> endometriosis -> EAOC axis (malignant incl. EAOC)",
}
NESTED = {"GSE213216", "GSE179640"}   # per-sample cellranger tarballs inside RAW.tar

def ftp_url(gse):
    stub = gse[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/{gse}_RAW.tar"

def fetch(gse):
    out = os.path.join(RAW, gse); os.makedirs(out, exist_ok=True)
    tar = os.path.join(RAW, f"{gse}_RAW.tar")
    if not os.path.exists(tar):
        print(f"[{gse}] downloading RAW.tar …")
        urllib.request.urlretrieve(ftp_url(gse), tar)
    with tarfile.open(tar) as t:
        t.extractall(out)
    if gse in NESTED:
        # each member is GSM..._sampleN.tar.gz containing */outs/filtered_feature_bc_matrix.h5
        for nz in glob.glob(os.path.join(out, "GSM*.tar.gz")):
            gsm = os.path.basename(nz).split("_")[0]
            with tarfile.open(nz) as tt:
                h5m = [m for m in tt.getmembers() if m.name.endswith("outs/filtered_feature_bc_matrix.h5")]
                if h5m:
                    tt.extract(h5m[0], out)
                    shutil.move(os.path.join(out, h5m[0].name), os.path.join(out, f"{gsm}.h5"))
            os.remove(nz)
        # tidy extracted sample folders
        for d in glob.glob(os.path.join(out, "*", "outs")):
            shutil.rmtree(os.path.dirname(d), ignore_errors=True)
    os.remove(tar)
    n_h5 = len(glob.glob(os.path.join(out, "*.h5")))
    n_mtx = len(glob.glob(os.path.join(out, "*matrix.mtx.gz")))
    print(f"[{gse}] {DATASETS[gse]} -> {n_h5} h5 / {n_mtx} mtx files")

if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    for gse in DATASETS:
        try: fetch(gse)
        except Exception as e: print(f"[{gse}] ERROR: {e}")
    print("done. Flat per-GSM mtx datasets keep GSM*_{barcodes,features,matrix}.mtx.gz; "
          "regroup into per-GSM folders before scanpy.read_10x_mtx if needed.")
