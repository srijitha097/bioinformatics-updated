import unittest,tempfile,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from report import association
class ReportTests(unittest.TestCase):
    def test_association_filters_invalid_rows(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'assoc.txt'
            file.write_text('CHR SNP BP P TEST\n1 rs1 200 0.01 ADD\n1 rs2 300 NA ADD\n1 rs1 200 0.5 DOMDEV\n')
            frame=association(file)
            self.assertEqual(len(frame),1)
            self.assertEqual(frame.iloc[0]['SNP'],'rs1')
if __name__=='__main__':unittest.main()
