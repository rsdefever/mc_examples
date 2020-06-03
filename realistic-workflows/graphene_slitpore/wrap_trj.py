from ovito.io import import_file, export_file
from ovito.modifiers import WrapPeriodicImagesModifier#, DeleteSelectedModifier

columns = ['Particle Type',
         'Position.X',
         'Position.Y',
         'Position.Z'
         ]

pipeline = import_file('nvt.out_wrapped.xyz', columns=columns,
        multiple_frames=True)

pipeline.modifiers.append(WrapPeriodicImagesModifier())
#pipeline.modifiers.append(DeleteSelectedModifier())

export_file(pipeline, 'wrapped.xyz', 'xyz', columns=columns, multiple_frames=True)
