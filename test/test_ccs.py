from falcon_unzip.mains import polish_gather as M

def test_is_haplotig():
    assert M.is_haplotig('/a/foo_bar.fa')
    assert not M.is_haplotig('/a/foo.fa')
