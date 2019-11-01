from pyputil.bin.structure_gen import run_cnt, run_flake, run_gnr
from snapshottest.file import FileSnapshot


def test_gnrs(snapshot, tmpdir):
    for args in [{
        "width": 4,
        "length": 12,
        "periodic": False,
        "zigzag": False,
        "no_hydrogen": False,
    }, {
        "width": 6,
        "length": 1,
        "periodic": True,
        "zigzag": False,
        "no_hydrogen": False,
    }, {
        "width": 3,
        "length": 1,
        "periodic": True,
        "zigzag": True,
        "no_hydrogen": False,
    }]:
        args["output"] = str(tmpdir.join("test.vasp"))
        # generate the GNR
        run_gnr(args)
        snapshot.assert_match(FileSnapshot(args["output"]))


def test_flakes(snapshot, tmpdir):
    output = str(tmpdir.join("test"))
    # generate the (three) flakes
    run_flake({
        "output": output,
        "radius": 7,
    })
    snapshot.assert_match(FileSnapshot(f"{output}-a.vasp"))
    snapshot.assert_match(FileSnapshot(f"{output}-b.vasp"))
    snapshot.assert_match(FileSnapshot(f"{output}-c.vasp"))


def test_cnts(snapshot, tmpdir):
    for args in [{
        "n": 5,
        "m": 5,
    }, {
        "n": 6,
        "m": 0,
    }, {
        "n": 2,
        "m": 5,
    }]:
        args["output"] = str(tmpdir.join("test.vasp"))
        # generate the CNT
        run_cnt(args)
        snapshot.assert_match(FileSnapshot(args["output"]))
