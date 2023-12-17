import psutil

def memory_usage():

    mem = psutil.virtual_memory()

    togb = 1e-9

    total = mem[0]*togb
    avail = mem[1]*togb
    used = mem[3]*togb
    free = mem[4]*togb
    active = mem[5]*togb
    inactive = mem[6]*togb

    print("Total:     ",total,"GB")
    print("Avail:     ",avail,"GB")
    print("Used:      ",used,"GB")
    print("Free:      ",free,"GB")
    print("Active:    ",active,"GB")
    print("Inactive:  ",inactive,"GB")
