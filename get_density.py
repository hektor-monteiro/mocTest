with open("old_grid_mod.f90", "r") as f:
    code = f.read()

start = code.find("read(77")
while start != -1:
    print("-----------------------------------")
    print(code[start-100:start+100])
    start = code.find("read(77", start + 1)
