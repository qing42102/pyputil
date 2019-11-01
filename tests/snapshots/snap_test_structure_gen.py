# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import Snapshot
from snapshottest.file import FileSnapshot


snapshots = Snapshot()

snapshots['test_gnrs 1'] = FileSnapshot('snap_test_structure_gen/test_gnrs 1.vasp')

snapshots['test_gnrs 2'] = FileSnapshot('snap_test_structure_gen/test_gnrs 2.vasp')

snapshots['test_gnrs 3'] = FileSnapshot('snap_test_structure_gen/test_gnrs 3.vasp')

snapshots['test_flakes 1'] = FileSnapshot('snap_test_structure_gen/test_flakes 1.vasp')

snapshots['test_flakes 2'] = FileSnapshot('snap_test_structure_gen/test_flakes 2.vasp')

snapshots['test_flakes 3'] = FileSnapshot('snap_test_structure_gen/test_flakes 3.vasp')

snapshots['test_cnts 1'] = FileSnapshot('snap_test_structure_gen/test_cnts 1.vasp')

snapshots['test_cnts 2'] = FileSnapshot('snap_test_structure_gen/test_cnts 2.vasp')

snapshots['test_cnts 3'] = FileSnapshot('snap_test_structure_gen/test_cnts 3.vasp')
