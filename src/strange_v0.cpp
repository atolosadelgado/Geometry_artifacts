//----------------------------------
//  pfRICH: Proximity Focusing RICH
//  Author: C. Dilks
//----------------------------------

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include <XML/Helper.h>

using namespace dd4hep;
using namespace dd4hep::rec;
using dd4hep::SubtractionSolid;

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>


TessellatedSolid CreatePyramid(double height_pyramid = 10 * cm, double aperture_deg = 60.)
{

  double base_side = height_pyramid * std::tan(3.14 / 180 * aperture_deg / 2);
  auto pxfc = [&](int i)
  { return base_side * std::cos(3.14 / 3. * i); };
  auto pyfc = [&](int i)
  { return base_side * std::sin(3.14 / 3. * i); };

  using Vertex = TessellatedSolid::Vertex;
  std::vector<Vertex> vertices;

  // point of pyramid, 0, 0, 0
  // base of pyramid on plane (x,y, z=height)
  // double height_pyramid = 10*cm;
  vertices.emplace_back(Vertex(0, 0, 0));
  vertices.emplace_back(Vertex(0, 0, height_pyramid));
  for (unsigned int i = 0; i < 6; i++)
  {
    vertices.emplace_back(Vertex(pxfc(i), pyfc(i), height_pyramid));
  }

  // vertex 0 is point of pyramid
  // vertex 1 is center of base of pyramid
  // vertex 2-7 are corners of base of pyramid
  TessellatedSolid shape("kk", vertices);
  for (unsigned int i = 2; i < 7; i++)
  {
    // definition of side of pyramid
    shape->AddFacet(0, i + 1, i);
    // definition of part of base of pyramid
    shape->AddFacet(1, i, i + 1);
  }
  
  //-- last facet added by hand
  // definition of side of pyramid
  shape->AddFacet(0, 2, 7);
  // definition of part of base of pyramid
  shape->AddFacet(1, 7, 2);

  shape->CloseShape(true, false, true);
  return shape;
}

// create the detector
static Ref_t createDummy3(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  OpticalSurfaceManager surfMgr = desc.surfaceManager();
  DetElement det(detName, detID);
  sens.setType("tracker");

  // Vessel
  Box vesselSolid(1 * cm, 1 * cm, 1 * cm);
  auto shape = CreatePyramid( 10*cm , 45);


  // Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
  Solid vesselSolid_final = IntersectionSolid(vesselSolid, shape);
  Volume vesselVol(detName, vesselSolid_final, desc.material("Aluminum"));

  // place mother volume (vessel)
  Volume motherVol = desc.pickMotherVolume(det);
  PlacedVolume vesselPV = motherVol.placeVolume(vesselVol, Position(0, 0, 0));
  vesselPV.addPhysVolID("system", detID);
  det.setPlacement(vesselPV);

  return det;
}
DECLARE_DETELEMENT(STRANGETYPE3, createDummy3)





