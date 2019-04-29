from falcon_unzip.mains import gen_rid_gathered_json as M
from falcon_kit.io import NativeIO as StringIO
import json

def test():
    sin = StringIO("""\
foo
bar
""")
    sout = StringIO()

    M.dump(sin, sout)

    got = json.loads(sout.getvalue())
    expected = [{'rid_to_phase_out': 'foo'}, {'rid_to_phase_out': 'bar'}]
    assert got == expected
