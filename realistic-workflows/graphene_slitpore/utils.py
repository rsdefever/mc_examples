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

def ovito_wrap(xyz_in):
    from ovito.io import import_file, export_file
    from ovito.modifiers import WrapPeriodicImagesModifier #, DeleteSelectedModifier
    """Wrap xyz trajectory through Ovito

    Note : `xyz_in` must have box information
    """

    columns = ['Particle Type',
             'Position.X',
             'Position.Y',
             'Position.Z'
             ]
    
    pipeline = import_file(xyz_in, columns=columns,
            multiple_frames=True)
    
    pipeline.modifiers.append(WrapPeriodicImagesModifier())
    #pipeline.modifiers.append(DeleteSelectedModifier())
   
    # Currently last frame is hard coded
    # TODO: Change this 
    last_frame = pipeline.compute(80)
    export_file(last_frame, 'wrapped.xyz', 'xyz', columns=columns)
