import stjames

print(stjames.BasisSet(name="STO-3G"))
print(stjames.BasisSet(name="STO-3G", cutoff_threshold=1e-10))
print(stjames.Settings.parse_obj({"basis_set": {"name": "STO-3G"}}))

print(stjames.Settings(method="hf-3c").level_of_theory)
print(stjames.Settings(method="hf").level_of_theory)
print(stjames.Settings(method="hf", corrections=["d3bj"]).level_of_theory)
print(stjames.Settings(method="b97d3", corrections=["d3bj"]).level_of_theory)
