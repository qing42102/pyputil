from pyputil.bin.modeplot import run
from snapshottest.file import FileSnapshot


def test_svgs(snapshot, datadir):
    # TODO: also check phonopy input ("qpoints.yaml")
    run({
        "input": datadir.join("POSCAR"),
        "eigs": datadir.join("gamma-dynmat.npz"),
        "supercell": [1, 1, 1],
        "config": None,
        "output": datadir,
        "gif": False,
        "parallel": None,
        "all_gifs": False,
    })

    for f in ["mode_3.svg", "mode_20.svg"]:
        output = str(datadir.join(f))
        snapshot.assert_match(FileSnapshot(output))


def test_gif(snapshot, datadir):
    run({
        "input": datadir.join("POSCAR"),
        "eigs": datadir.join("gamma-dynmat.npz"),
        "supercell": [1, 1, 1],
        "config": None,
        "output": datadir,
        "gif": 3,
        "parallel": None,
        "all_gifs": False,
    })

    output = str(datadir.join("mode_3.gif"))
    snapshot.assert_match(FileSnapshot(output))

