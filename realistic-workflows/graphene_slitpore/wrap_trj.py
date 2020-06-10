import mbuild as mb
from ovito.io import import_file, export_file
from ovito.modifiers import WrapPeriodicImagesModifier

columns = ['Particle Type',
         'Position.X',
         'Position.Y',
         'Position.Z'
         ]

pipeline = import_file('nvt.out_wrapped.xyz', columns=columns,
        multiple_frames=True)

pipeline.modifiers.append(WrapPeriodicImagesModifier())

export_file(pipeline, 'wrapped.xyz', 'xyz', columns=columns, multiple_frames=True)
