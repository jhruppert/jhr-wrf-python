import psutil
import numpy as np

def memory_usage():

    mem = psutil.virtual_memory()

    togb = 1e-9

    nround=2

    total    = np.round(mem[0]*togb, nround)
    avail    = np.round(mem[1]*togb, nround)
    used     = np.round(mem[3]*togb, nround)
    free     = np.round(mem[4]*togb, nround)
    active   = np.round(mem[5]*togb, nround)
    inactive = np.round(mem[6]*togb, nround)

    print()
    print("Total:     ",total,"GB")
    print("Avail:     ",avail,"GB")
    print("Used:      ",used,"GB")
    print("Free:      ",free,"GB")
    print("Active:    ",active,"GB")
    print("Inactive:  ",inactive,"GB")
    print()
