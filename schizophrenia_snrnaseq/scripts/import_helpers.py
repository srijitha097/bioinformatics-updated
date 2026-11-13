"""Pure metadata helpers that can be tested without downloading real data."""
import re
def diagnosis_group(value):
    # Fail closed: do not call bipolar, MDD, or other diagnosis a control.
    s=str(value).strip().casefold()
    if re.search(r'\bschizophrenia\b|\bschizoaffective\b|\bscz\b|\bsz\b',s):
        return 'case' if 'schizoaffective' not in s else 'other'
    if re.fullmatch(r'control|healthy|healthy control|unaffected|unaffected control|hc|con',s):
        return 'control'
    return 'other'

