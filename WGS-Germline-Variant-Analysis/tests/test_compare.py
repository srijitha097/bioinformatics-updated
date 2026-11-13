import tempfile,unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_vcfs import compare
from coverage import summarize

class CompareTest(unittest.TestCase):
    def test_allele_concordance_and_bed(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            header='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            (d/'truth.vcf').write_text(header+'chr20\t102\t.\tA\tC\t50\tPASS\t.\nchr20\t103\t.\tG\tT\t50\tPASS\t.\n')
            (d/'call.vcf').write_text(header+'chr20\t102\t.\tA\tC\t50\tPASS\t.\nchr20\t104\t.\tA\tG\t50\tPASS\t.\n')
            (d/'hc.bed').write_text('chr20\t100\t105\n')
            obj=compare(d/'call.vcf',d/'truth.vcf',d/'hc.bed','chr20:100-110')
            self.assertEqual((obj['TP'],obj['FP'],obj['FN']),(1,1,1))
            self.assertEqual(obj['precision'],0.5)
    def test_depth(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'coverage.tsv';p.write_text('chr20\t1\t0\nchr20\t2\t20\n')
            self.assertEqual(summarize(p)['pct_ge_20x'],50)
if __name__=='__main__':unittest.main()
