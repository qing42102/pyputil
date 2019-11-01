# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import Snapshot
from snapshottest.file import FileSnapshot


snapshots = Snapshot()

snapshots['test_svgs 1'] = FileSnapshot('snap_test_modeplot/test_svgs 1.svg')

snapshots['test_svgs 2'] = FileSnapshot('snap_test_modeplot/test_svgs 2.svg')

snapshots['test_gif 1'] = FileSnapshot('snap_test_modeplot/test_gif 1.gif')
