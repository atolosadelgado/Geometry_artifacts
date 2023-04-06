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

#include "DD4hep/DetFactoryHelper.h"
#include "TGeoArb8.h"
#include <iomanip>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;




UnionSolid createdummy(double height_pyramid = 10 * cm, double aperture_deg = 10.)
{


    double r_out = 210*cm; // height_pyramid + 20*cm;
    double r_in  = 190*cm; // height_pyramid;
    double phistep = 13.333 * deg;
    aperture_deg = phistep/2;
    double r_ave = 0.5*(r_out+r_in);
    double r_thickness = r_out - r_in;


    double side_out = 148.15*mm; //r_out * std::tan( M_PI / 180 * aperture_deg / 2) /2 ;
    double side_in  = side_out*r_in/r_out; //r_out * std::tan( M_PI / 180 * aperture_deg / 2) /2 ;

//     std::cout << side_out/cm << std::endl;

//     std::cin.ignore();

//     double side_in  = r_in  * std::tan( M_PI / 180 * aperture_deg / 2) / 2;
    double apothem_out = side_out*sin( 60*deg );

    auto pxfc_in = [&](int i)
    { return side_in * std::cos( M_PI / 3. * i); };

    auto pxfc = [&](int i)
    { return side_out * std::cos( M_PI / 3. * i); };
    auto pyfc = [&](int i)
    { return side_out * std::sin( M_PI / 3. * i); };

    using Vertex = TessellatedSolid::Vertex;
    std::vector<Vertex> vertices;
    vertices.reserve(16);

    for (unsigned int i = 0; i < 6; i++)
    {
        vertices.emplace_back(Vertex(pxfc(i), pyfc(i), r_in));
        vertices.emplace_back(Vertex(pxfc(i), pyfc(i), r_out));
    }
    vertices.emplace_back(Vertex(0, 0, r_in));
    vertices.emplace_back(Vertex(0, 0, r_out));


    double v[8][2];
    for( int ivx = 0; ivx<4; ++ ivx)
    {
        v[ivx][0] = pxfc(ivx);
        v[ivx][1] = pyfc(ivx);
        v[ivx+4][0] = pxfc_in(ivx);
        v[ivx+4][1] = pyfc(ivx);
    }
    EightPointSolid right_side( r_thickness ,&v[0][0]);

    double v2[8][2];
    for( int ivx = 0; ivx<4; ++ ivx)
    {
        v2[ivx][0] = pxfc(ivx + 3);
        v2[ivx][1] = pyfc(ivx + 3);
        v2[ivx+4][0] = pxfc_in(ivx + 3);
        v2[ivx+4][1] = pyfc(ivx + 3);
    }
    EightPointSolid left_side( r_thickness ,&v2[0][0]);


    return UnionSolid(left_side, right_side);


}



TessellatedSolid CreatePyramid_v2(double height_pyramid = 10 * cm, double aperture_deg = 60.)
{

    double r_out = height_pyramid + 1*cm;
    double r_in  = height_pyramid;

    double side_out = r_out * std::tan( M_PI / 180 * aperture_deg / 2);
    double apothem_out = side_out*sin( 60*deg );

    auto pxfc = [&](int i)
    { return side_out * std::cos( M_PI / 3. * i); };
    auto pyfc = [&](int i)
    { return side_out * std::sin( M_PI / 3. * i); };

    using Vertex = TessellatedSolid::Vertex;
    std::vector<Vertex> vertices;
    vertices.reserve(16);

    for (unsigned int i = 0; i < 6; i++)
    {
        vertices.emplace_back(Vertex(pxfc(i), pyfc(i), r_in));
        vertices.emplace_back(Vertex(pxfc(i), pyfc(i), r_out));
    }
    vertices.emplace_back(Vertex(0, 0, r_in));
    vertices.emplace_back(Vertex(0, 0, r_out));




    TessellatedSolid shape("kk", vertices);

    for (unsigned int i = 1; i <= 11; i+=2 )
    {
        // top base
        shape->AddFacet(13, i , ((i+2) % 12) );
        // bottom
        shape->AddFacet(12,  ((i+1) % 12) , i-1  );
    }

    // sides of the cell
    for(unsigned int side_n = 0; side_n<6; ++side_n)
    {
        int vx_bottom_left  = (2*side_n) % 12;
        int vx_top_left     = (vx_bottom_left+1 ) % 12 ;
        int vx_top_right    = (vx_bottom_left+3 ) % 12;
        int vx_bottom_right = (vx_bottom_left+2 ) % 12;
        shape->AddFacet( vx_bottom_left, vx_bottom_right, vx_top_right, vx_top_left );
    }




    shape->CloseShape(true, false, true);
    return shape;
}
// create the detector
static Ref_t tricky_pyramidal_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    xml::Component dims = detElem.dimensions();
    DetElement det(detName, detID);
    sens.setType("tracker");

    // Vessel
//     Box vesselSolid(1 * cm, 1 * cm, 1 * cm);
    auto shape = createdummy( 190*cm , 13.333);


    // Define the actual mirror as intersection of the mother volume and the hollow sphere just defined
    //   Solid vesselSolid_final = IntersectionSolid(vesselSolid, shape);

    Volume vesselVol(detName+"vol", shape, desc.material("Aluminum"));

    // place mother volume (vessel)
    Volume motherVol = desc.pickMotherVolume(det);

    Transform3D pyrTr1 =  RotationX( -90. * deg) * RotationY( -90. * deg) * Translation3D(0, 0, -190*cm);
    PlacedVolume vesselPV = motherVol.placeVolume(vesselVol, pyrTr1);
    vesselPV.addPhysVolID("system", detID);
    det.setPlacement(vesselPV);


      double hexagon_side_length = 14.815 * cm;
  /// Distance in x-direction
  double zstep = 2 * hexagon_side_length*sin( 60*deg );;
  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 13.333 * deg;

     Translation3D pyrTr2 (0, 0, zstep);
     auto pyrTr3 = RotationZ(  phistep*0.5) *Translation3D(0, 0, zstep/2);
    motherVol.placeVolume(vesselVol, pyrTr2*pyrTr1);
    motherVol.placeVolume(vesselVol, pyrTr3*pyrTr1);

    det.setPlacement(vesselPV);



    return det;
}
DECLARE_DETELEMENT(TRICKY_PYR_TYPE, tricky_pyramidal_cell)


