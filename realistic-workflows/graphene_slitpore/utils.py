import numpy as np

def get_last_frame(xyz_in, last_frame):
    """Get last frame of the xyz trajectory
    """
    with open(xyz_in, 'r') as fi:
        for line_num, line_in in enumerate(fi):
            if 'MC_STEP' in line_in and str(last_frame) in line_in:
                 break
            #if line_num == 222398:
            #     import pdb; pdb.set_trace()
  

    last_lines = open(xyz_in, 'r').readlines()[line_num-1:]
    with open('last_frame.xyz', 'w') as fo:
        for line in last_lines:
            fo.write(line)
