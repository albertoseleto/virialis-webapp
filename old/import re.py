import re

txt = "Redundant internal coordinates taken from checkpoint file:\n d:\\programas ic\\gaussian\\Scratch\\gxx.chk\n Charge =  0 Multiplicity = 1\n O,0,-0.0007781877,0.,-0.0006071343\n H,0,0.0017065327,0.,0.9612549545\n H,0,0.9327992298,0.,-0.2321628667\n Recover connectivity data from disk."

values = re.findall(r"-?\d+\.\d+", txt)
print(values)