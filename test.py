import stjames

print(stjames.BasisSet(name="STO-3G"))
print(stjames.BasisSet(name="STO-3G", cutoff_threshold=1e-10))
print(stjames.Settings.parse_obj({"basis_set": {"name": "STO-3G"}}))
