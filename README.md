GEOMETRY ARTIFACTS
==================

Minimal working example to reproduce visual artifact when intersection is applied on a tessellated solid. In this example, the intersection of a pyramid and a cube result on a shorter pyramid (as expected) plus some artifact.

# Basic commands

Build with `cmake`, 
```bash
source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
cmake -B build -S . -D CMAKE_INSTALL_PREFIX=install
cmake --build build -- install
export LD_LIBRARY_PATH=$PWD/install/lib:$LD_LIBRARY_PATH
```

To display the geometry using ROOT display

```bash
geoDisplay compact/arc_v0.xml 
```

To convert the geometry into GDML (output file also part of the repository)

```bash
geoConverter -compact2gdml -input compact/strange_v0.xml -output strange.gdml
```

Visualze geometry using Geant4+Qt

```bash
ddsim --compactFile compact/strange_v0.xml --runType qt --macroFile vis.mac --part.userParticleHandler=''
```

A screenshot of the G4+Qt visualization can be seen in the file `geometry_artifact.png`.

To show materials in along a given line (in this case from (0,0,-10)cm to (100,0,90)cm), only air (world material) is shown. This may point to an issue just with visualization.

```bash
materialScan compact/strange_v0.xml 0 0 -10 100 0 90
```

Raytracing option in ROOT do not remove the artifacts. I have not tried Ray-tracing display option with G4+Qt


==================

# Solution to the problem (thanks to Evgueni Tcherniaev)

The tessellated solid was constructed incorrectly. When you run Geant4 application, there is the following warning message :

```
-------- WWWW ------- G4Exception-START -------- WWWW -------

*** ExceptionHandler is not defined ***
*** G4Exception : GeomSolids1001
      issued by : G4TessellatedSolid::SetSolidClosed()
Defects in solid: kk - some facets have wrong orientation!
*** This is just a warning message. ***
-------- WWWW ------- G4Exception-END -------- WWWW -------
```

The normals of the facets in a tessellated solid should point to outside of the solid. It means that the vertices of the facets should be set in the anti-clockwise order looking from outside of the solid. In this case all lateral facets were oriented incorrectly.
