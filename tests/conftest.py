from __future__ import unicode_literals

from distutils import dir_util
from pathlib import Path
from pytest import fixture


# ref: https://stackoverflow.com/a/1994840
@fixture
def datadir(tmpdir, request):
    """
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    """
    path = Path(request.module.__file__)
    directory = path.parent / path.stem

    if directory.is_dir():
        dir_util.copy_tree(directory, str(tmpdir))

    return tmpdir
