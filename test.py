import stjames
from stjames.correction import Correction
from stjames.method import Method

print(stjames.BasisSet(name="STO-3G"))
print(stjames.BasisSet(name="STO-3G", cutoff_threshold=1e-10))
print(stjames.Settings.parse_obj({"basis_set": {"name": "STO-3G"}}))

print(stjames.Settings(method="b97-3c").level_of_theory)
print(stjames.Settings(method="hf-3c").level_of_theory)
print(stjames.Settings(method=Method.HF3C, corrections=[Correction.D3BJ]).level_of_theory)
print(stjames.Settings(method=Method.B97D3, corrections=[Correction.D3BJ]).level_of_theory)
