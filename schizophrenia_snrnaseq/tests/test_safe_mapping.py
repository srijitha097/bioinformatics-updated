import importlib.util
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]/'scripts/import_helpers.py'
spec=importlib.util.spec_from_file_location('import_helpers',HERE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_diagnosis_does_not_misclassify_other_diseases():
 assert m.diagnosis_group('Schizophrenia')=='case'
 assert m.diagnosis_group('Healthy control')=='control'
 assert m.diagnosis_group('bipolar disorder')=='other'
 assert m.diagnosis_group('schizoaffective disorder')=='other'
 assert m.diagnosis_group('major depressive disorder')=='other'
