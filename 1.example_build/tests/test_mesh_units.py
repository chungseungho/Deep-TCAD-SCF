"""Check that only Gmsh ASCII 2.2 node coordinates are rescaled."""
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry.make_mesh import scale_msh2_nm_to_cm


def test_msh_scale():
    src='$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n2\n1 10 0 0\n2 20 0 0\n$EndNodes\n$Elements\n1\n1 1 2 3 4 1 2\n$EndElements\n'
    with TemporaryDirectory() as d:
        a=Path(d)/'nm.msh';b=Path(d)/'cm.msh'
        a.write_text(src)
        scale_msh2_nm_to_cm(a,b)
        raw=b.read_text()
        assert '$Elements\n1\n1 1 2 3 4 1 2\n' in raw
        line=raw.split('$Nodes')[1].splitlines()[2]
        assert abs(float(line.split()[1])-1e-6)<1e-20


if __name__=='__main__':
    test_msh_scale()
    print('PASS: msh2 coordinate scaling nm -> cm preserves physical element tags')
