"""
We work with files rather than local strings because task functions
depend on actual files. But everything is under the test_data/ directory.
"""
from falcon_unzip.tasks import top as MOD
from falcon_unzip.io import deserialize
import pytest

def diff(expected_fn, got_fn):
    expected = deserialize(str(expected_fn))
    got = deserialize(str(got_fn))
    assert expected == got, '"{}" != "{}"'.format(expected_fn, got_fn)

@pytest.fixture(scope='module')
def datadir(request):
    """Return '../test_data/' directory.
    """
    return request.fspath.join('..', '..', 'test_data')


def test_p_ctg_fai2ctgs(datadir, tmpdir):
    p_ctg_fai_fn = str(datadir.join('greg200k-sv2-ccs/2-asm-falcon/p_ctg.fa.fai'))
    o_fn = str(tmpdir.join('CTGS.json'))

    MOD.p_ctg_fai2ctgs(p_ctg_fai_fn, o_fn)

    expected_fn = datadir.join('greg200k-sv2-ccs/expected/CTGS.json')
    diff(expected_fn, o_fn)
