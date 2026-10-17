#!/usr/bin/env python3
"""Verify frozen original input bytes; obtain inputs lawfully before running analyses."""
import hashlib,sys,pathlib
HASH={
'four_income_group_final_dataset.csv':'5dfe182725b16192ca00d93d7c34e900471e8b11855542fd6319cf08b8b958e6',
'irena_vintage_overlap.csv':'fdd711c71b851b50754d4bde3f266bd497532f60db922a122c49b43b1b40b5af'
}
if len(sys.argv)!=3:
    raise SystemExit("Usage: python verify_frozen_inputs.py /authorized/four_income_group_final_dataset.csv /authorized/irena_vintage_overlap.csv")
for s in sys.argv[1:]:
    p=pathlib.Path(s)
    name=p.name
    if name not in HASH: raise SystemExit("Unexpected input filename: "+name)
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h != HASH[name]: raise SystemExit("FROZEN HASH MISMATCH for "+name+": "+h)
    print("PASS",name,h)
