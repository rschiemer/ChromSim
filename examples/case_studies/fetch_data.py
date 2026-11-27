"""Download the public data sets used by the case studies into ``examples/case_studies/data/``.

    pixi run python examples/case_studies/fetch_data.py            # all data sets
    pixi run python examples/case_studies/fetch_data.py fraunhofer_lysozyme

The data belong to their authors and are **not** part of ChromSim's repository: every file is
downloaded from a pinned commit of its public GitHub repository and checked against its SHA-256,
so the case studies always see exactly the files they were written for. ``data/`` is git-ignored.
Respect each data set's license (below); several are for non-commercial use only.

Data sets
---------
unicorn5_res
    Native UNICORN 5 result files (.res): ÄKTA prime sample from PyCORN (GPL-2.0) and three
    Superdex 200 16/60 SEC runs from PyCornGUI (GPL-2.0).
unicorn7_csv
    UNICORN 7 CSV exports (UTF-16, tab): SEC run from Appia (MIT); three HisTrap (IMAC) +
    Superdex 75 runs from easy_aekta (no license stated: local testing only, do not redistribute).
unicorn5_asc
    UNICORN 5 .asc/.csv exports from chromr (MIT).
fraunhofer_lysozyme
    Measured acetone/dextran tracer pulses (with and without column, triplicates) and lysozyme
    gradient elutions (0–120 CV) on a strong cation exchanger. S. Menzel & J. Buyel, Fraunhofer
    IME, shipped with CADET-MCMC. CC BY-NC-SA 4.0 (non-commercial, attribution, share-alike).
cadet_workshop_yamamoto
    Simulated (noise-free) linear-gradient elutions of three proteins at four gradient volumes,
    from the CADET-Workshop lesson 99.1. CC BY-NC-SA 4.0.
"""

import hashlib
import sys
import urllib.parse
import urllib.request
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"

LICENSES = {
    "unicorn5_res": "GPL-2.0 (PyCORN, PyCornGUI sample files)",
    "unicorn7_csv": "MIT (Appia); no license stated (easy_aekta) - local testing only",
    "unicorn5_asc": "MIT (chromr)",
    "fraunhofer_lysozyme": "CC BY-NC-SA 4.0, (c) S. Menzel, J. Buyel (Fraunhofer IME)",
    "cadet_workshop_yamamoto": "CC BY-NC-SA 4.0 (CADET-Workshop)",
}

# (data set, GitHub repository, commit, path in the repository, SHA-256 of the file)
FILES = [
    (
        "unicorn5_res",
        "pyahmed/PyCORN",
        "d7747b6eaff818f07b5f41e950b65eb716ba4e0f",
        "samples/sample1.res",
        "15c56238a14a2bc58ee2f1707031b24661c086d87092a94cf53eb95388e27f3d",
    ),
    (
        "unicorn5_res",
        "kotarokelley/PyCornGUI",
        "98143df4810ebc911554088b96c8629e573dbc41",
        "PyCornGUI/samples/KOKO038032 1660 S200.res",
        "b9b98e2c6f9e32edc426d4f539aecc1eaf0f2226d93d6b8a4d7a20562b254c5d",
    ),
    (
        "unicorn5_res",
        "kotarokelley/PyCornGUI",
        "98143df4810ebc911554088b96c8629e573dbc41",
        "PyCornGUI/samples/KOKO135159 cut ON 1660 S200.res",
        "d777072c8605eed6984f2eebdbcbcdf2619e518c3d06433de3f30dfb938fe3be",
    ),
    (
        "unicorn5_res",
        "kotarokelley/PyCornGUI",
        "98143df4810ebc911554088b96c8629e573dbc41",
        "PyCornGUI/samples/160502 KOKO029125005099 1660 S200 2001.res",
        "6cdd3a4a21b41837bab67ef4e69e117f2cb311dd72377f3368484743fd694fc9",
    ),
    (
        "unicorn7_csv",
        "PlethoraChutney/Appia",
        "a001b2c0522ef28df480d98a514a35e2dd91191c",
        "test-files/2018_0821SEC_detergentENaC.csv",
        "f6a4dfbb9fd6ca12961c5b4ed3a19630574aa03697cc9303cbb5a80d314feda5",
    ),
    (
        "unicorn7_csv",
        "Bloeci/easy_aekta",
        "e542b2120a26ea5e970c572eab2588a5ed767d6b",
        "test_input/test_rey/20260209_RT_ebo_scaf_20260206_001-024_10mL_HisTrapExcel1mL_SD75increase 001.csv",  # noqa: E501  (upstream file name)
        "d3c696fd20418512c5e9449dae635d262f449fe9f4924d64c1d951dd4a948b3b",
    ),
    (
        "unicorn7_csv",
        "Bloeci/easy_aekta",
        "e542b2120a26ea5e970c572eab2588a5ed767d6b",
        "test_input/test_rey/20260209_RT_ebo_scaf_20260206_001-024_10mL_HisTrapExcel1mL_SD75increase 002.csv",  # noqa: E501  (upstream file name)
        "33999c6011e1ad7d5eaab4a354ad9dc5a433886730543634f6d4853f999d9842",
    ),
    (
        "unicorn7_csv",
        "Bloeci/easy_aekta",
        "e542b2120a26ea5e970c572eab2588a5ed767d6b",
        "test_input/test_rey/20260209_RT_ebo_scaf_20260206_001-024_10mL_HisTrapExcel1mL_SD75increase 003.csv",  # noqa: E501  (upstream file name)
        "c164535a69051e751880b367403c47236786891a0a17ce44a29202ec9c88dada",
    ),
    (
        "unicorn5_asc",
        "BradyAJohnston/chromr",
        "7e76afdb97054fedd738017a6b547e3e3fe0b2ca",
        "inst/extdata/example_akta_1.asc",
        "e71a3a4f3a6797c32ec23482cf9f4eb69fd2dbb636e199cbb1c5a9f75163d494",
    ),
    (
        "unicorn5_asc",
        "BradyAJohnston/chromr",
        "7e76afdb97054fedd738017a6b547e3e3fe0b2ca",
        "inst/extdata/example_akta_2.asc",
        "6736b846454f58434a5346f17bfef2d8e52403725dd993480ce9533ca447240a",
    ),
    (
        "unicorn5_asc",
        "BradyAJohnston/chromr",
        "7e76afdb97054fedd738017a6b547e3e3fe0b2ca",
        "inst/extdata/example_akta_3.csv",
        "3f773c9637f018cf659f601e346b05bc194be2bc1751309b6a9fe0837bab4cdf",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/bypassAcetone1.dat",
        "776fea78de8ce84a04167e2326e93dc1a63c545814af374db4312b10ae9c679d",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/bypassAcetone2.dat",
        "85639dd358b654f2f5a75194ae5f8fbb546ce9fa7e368ac12eb4a852bd228e72",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/bypassAcetone3.dat",
        "776fea78de8ce84a04167e2326e93dc1a63c545814af374db4312b10ae9c679d",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/bypassDextran1.dat",
        "b17547f14491564a5c7983771d023bc903fee57a53994251779ad7ad102230bf",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/bypassDextran2.dat",
        "f46f09fd9b28f381765b992baf52830335b0596af31024329061ed08632fa697",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/bypassDextran3.dat",
        "54da30ad22a1e68740a2099b73f832202f2961ce365eb703c8e14837b51b7db1",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnAcetone1.dat",
        "9a030f1c57e946e7f12a34d1ebedb18f7a2ba6ef30e074c75df9b65d9d77fdea",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnAcetone2.dat",
        "fc59704b3e392cc9849c83d59bdfc00c4abd0b3f7a1cb0c49fee446f8c92a403",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnAcetone3.dat",
        "f6451c9bf8a43ee4a73fa58174e3f22a1eb18771aad230be47097d321a64bc1f",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnDextran1.dat",
        "adf2d0cfb500d318796dccfc8e16395cd74b1b12793dc12d29e1ce920b6916b3",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnDextran2.dat",
        "6ee02607301febc0d5e23cf96647ad65ffff833af9d3529fd3955f2a9fc2548c",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnDextran3.dat",
        "1b0e177fb569cd4f45a485eb28f5427601f6301f5da5d7e9807bb744480d4b50",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnLysozyme0cv.dat",
        "0d5396aceb0e411211829ff365aa0f590486f678b7d5045fd21df1826a54088a",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnLysozyme5cv.dat",
        "2d33eb2d97feb70d30fb7f5c834b0c4be2df4f355cf7a2a7061c650e1fa27fca",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnLysozyme10cv.dat",
        "1b524b591bc08cad76bdbd1d1e887818c606ba3f669feec6f955e96ac3b4706d",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnLysozyme30cv.dat",
        "e45f1f8a35f8289ee29689f2e45419a623eb22f0ccfb87e6214a7c42db05bbaa",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnLysozyme60cv.dat",
        "d6915bad2e26a3eaec26cc6350713b863d9cefe16a27f7aba265deec6a553c47",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/columnLysozyme120cv.dat",
        "5298b12eaef6d51a2c4827377d953f6212b7f6de86602e6551ba8053c281e0c8",
    ),
    (
        "fraunhofer_lysozyme",
        "cadet/CADET-MCMC",
        "ccc8002a7a517d74103828011ee2436d44ffaa0e",
        "dataset/README.md",
        "9dd1b8db9302372d3407651023713564a72fddd108b76015c887c65bebc8233b",
    ),
    (
        "cadet_workshop_yamamoto",
        "cadet/CADET-Workshop",
        "7cf039e8e744c7b2c1144c607b959d0aa549d600",
        "99 Misc/99.1 Yamamoto's Method/experiments/three_proteins/18.8mL.csv",
        "15ce03ef4853b0dddd38378204ba885300e8c0641c9cddc2b4ef6fc3415deb50",
    ),
    (
        "cadet_workshop_yamamoto",
        "cadet/CADET-Workshop",
        "7cf039e8e744c7b2c1144c607b959d0aa549d600",
        "99 Misc/99.1 Yamamoto's Method/experiments/three_proteins/37.6mL.csv",
        "f8fc7891488f2e43f39830ec5673b7d3ccbf58d91b9d22b214ed6ec6c316ad3f",
    ),
    (
        "cadet_workshop_yamamoto",
        "cadet/CADET-Workshop",
        "7cf039e8e744c7b2c1144c607b959d0aa549d600",
        "99 Misc/99.1 Yamamoto's Method/experiments/three_proteins/56.4mL.csv",
        "c0ba2775b1bd6a9f10f502da40c4c83419c20a9a515ed1821745b26d082cfe8a",
    ),
    (
        "cadet_workshop_yamamoto",
        "cadet/CADET-Workshop",
        "7cf039e8e744c7b2c1144c607b959d0aa549d600",
        "99 Misc/99.1 Yamamoto's Method/experiments/three_proteins/75.2mL.csv",
        "e877ff7515049a7f7ff063be51e865ae455bd759d0ee012b22bcf9210336ccde",
    ),
]


def raw_url(repository: str, commit: str, path: str) -> str:
    """Download URL of a file at a fixed commit."""
    return f"https://raw.githubusercontent.com/{repository}/{commit}/{urllib.parse.quote(path)}"


def fetch(datasets: list[str] | None = None) -> None:
    """Download (or verify) the files of the given data sets; default all of them."""
    wanted = datasets or sorted(LICENSES)
    unknown = sorted(set(wanted) - set(LICENSES))
    if unknown:
        raise SystemExit(f"Unknown data set(s) {unknown}; choose from {sorted(LICENSES)}.")
    for dataset, repository, commit, path, sha256 in FILES:
        if dataset not in wanted:
            continue
        target = DATA / dataset / Path(path).name
        if target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest() == sha256:
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(raw_url(repository, commit, path), timeout=60) as response:
            content = response.read()
        if hashlib.sha256(content).hexdigest() != sha256:
            raise SystemExit(f"Checksum mismatch for {path} from {repository}; not saved.")
        target.write_bytes(content)
        print(f"{dataset}: {target.name}")
    for dataset in wanted:
        sources = sorted(
            {f"https://github.com/{r}/tree/{c}" for d, r, c, _, _ in FILES if d == dataset}
        )
        (DATA / dataset / "SOURCE.txt").write_text(
            f"License: {LICENSES[dataset]}\nSource(s):\n" + "\n".join(sources) + "\n"
        )
    print(f"Data in {DATA}")


if __name__ == "__main__":
    fetch(sys.argv[1:])
