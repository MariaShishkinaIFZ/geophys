===========================================
The 3d spatial data codec tool: 3d-codec.py
===========================================

The ``3d-codec.py`` utility is dedicated to the conversion
from and to various spatial data formats used currently
in the geophysics framework.

Currently supported formats:

* ASCII text file: one spatial point per line, spatial point
  contains spatial information ``(X, Y, Z)`` and a data value ``(D)``
  at given coordinates. ``(X,Y,Z,D)`` columns order is arbitraty
  (the order is defined by ``--iformat`` or ``--oformat`` option), other
  columns can present and will be ignored.

  .. note::
     When reading the ASCII, incoming points spatial distribution
     is verified to fit a spatial grid with orthogonal axis with uniform
     intervals between grid nodes (axes spatial scales and 0s may
     differ from one axis to another).

  .. note::
     When reading the ASCII the spatial grid parameters are extracted
     and tracked internally by the tool even if the output format
     lacks the spatial data information.

* Binary file: point data values (float) are located sequentially
  in the file in the order defined (``--iformat`` or ``--oformat``).

  .. note::
     When writing the binary file the spatial grid information is
     lost and needed to be tracked separately for now.

Usage Examples
--------------

To convert from ``./in.txt`` ascii format with column order ``zyxd`` to
the binary file ``./bin.dat`` with index change oder: ``zyx`` use following:

.. code-block:: bash

   python3 3d-codec.py --in ./in.txt --out ./bin.dat     \
      --itype ascii --iformat zyxd --otype binary --oformat "zyx"

To convert from ``./bin.dat`` binary format with index change order ``zyx``
and ``ZYX`` sizes ``(2;2;3)`` to the asci file ``./out.txt`` with column oder:
``zyx``, use following:                                       \n'

.. code-block:: bash

   python3 ../3d-codec.py --in ./bin.dat --out ./out.txt       \
         --itype binary --iformat "z(2)y(2)x(3)" --otype ascii \
         --oformat zyxd
