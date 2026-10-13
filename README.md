# Rutilated quartz

Orientation of rutile needles in quartz, from work with Nick Tailby (AMNH) in 2019.

Needle directions were fitted as lines in a 3D image volume of one quartz grain (sample Jill m2 g5; `Jill_M2_g5/VG fit Jill m2g5.csv`), and the grain's crystal orientation was measured by EBSD (Euler angles from EDAX OIM). `EBSD_Sample_ReferenceFrames.py` rotates vectors between the sample and crystal reference frames. `Data_processing_example.ipynb` plots the needle directions on stereonets in both frames, with density contours.

```bash
uv sync
uv run --with jupyter jupyter lab Data_processing_example.ipynb
```

`Notes on EBSD Directions.txt` records how the EBSD stage axes map onto the sample axes.

## Status

This is a working research repository from 2019, kept for reference. The two `EBSD Rotation Matrix` scripts are earlier Python 2 versions of the reference-frame code.

## License

The code is under the MIT License (see LICENSE). The data, images and figures come from the collaboration and are not covered by it.
