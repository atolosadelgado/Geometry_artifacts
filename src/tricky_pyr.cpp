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

// using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include <XML/Helper.h>

using namespace dd4hep;
using namespace dd4hep::rec;
using dd4hep::SubtractionSolid;



TessellatedSolid CreatePyramid_v2 ( double r_in = 190 * cm, double r_out = 210*cm /*, double aperture_deg = 13.333/2 * deg*/ )
{

    double r_ave = 0.5* ( r_out+r_in );
    double r_thickness = r_out - r_in;


    double side_out = 148.15*mm;
    double side_in  = side_out*r_in/r_out;

    double apothem_out = side_out*sin ( 60*deg );


    auto pxfc_in = [&] ( int i ) {
        return side_out * std::cos ( M_PI / 3. * i );
    };

    auto pxfc = [&] ( int i ) {
        return side_out * std::cos ( M_PI / 3. * i );
    };
    auto pyfc = [&] ( int i ) {
        return side_out * std::sin ( M_PI / 3. * i );
    };

    using Vertex = TessellatedSolid::Vertex;
    std::vector<Vertex> vertices;
    vertices.reserve ( 16 );

    for ( unsigned int i = 0; i < 6; i++ ) {
        vertices.emplace_back ( Vertex ( pxfc_in ( i ), pyfc ( i ), r_in ) );
        vertices.emplace_back ( Vertex ( pxfc ( i ), pyfc ( i ), r_out ) );
    }
    vertices.emplace_back ( Vertex ( 0, 0, r_in ) );
    vertices.emplace_back ( Vertex ( 0, 0, r_out ) );




    TessellatedSolid shape ( "kk", vertices );

    for ( unsigned int i = 1; i <= 11; i+=2 ) {
        // top base
        shape->AddFacet ( 13, i, ( ( i+2 ) % 12 ) );
        // bottom
        shape->AddFacet ( 12, ( ( i+1 ) % 12 ), i-1 );
    }

    // sides of the cell
    for ( unsigned int side_n = 0; side_n<6; ++side_n ) {
        int vx_bottom_left  = ( 2*side_n ) % 12;
        int vx_top_left     = ( vx_bottom_left+1 ) % 12 ;
        int vx_top_right    = ( vx_bottom_left+3 ) % 12;
        int vx_bottom_right = ( vx_bottom_left+2 ) % 12;
        shape->AddFacet ( vx_bottom_left, vx_bottom_right, vx_top_right, vx_top_left );
    }




    shape->CloseShape ( true, false, true );
    return shape;
}
// create the detector
static Ref_t tricky_pyramidal_cell ( Detector &desc, xml::Handle_t handle, SensitiveDetector sens )
{
    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();
    xml::Component dims = detElem.dimensions();
    DetElement det ( detName, detID );
    sens.setType ( "tracker" );

    // Vessel


    double r_in = 190 * cm;
    double r_out = 210*cm;
    double hex_side = 148.15*mm;

    double hex_aperture = atan( hex_side / r_out ) *rad;

    double hex_side_out = 0.99* hex_side;// r_out*sin( hex_aperture ); //hex_side;//
    double hex_side_in = hex_side_out *r_in/r_out;

    std::vector<double> zplanes = { r_in, r_out *cos(hex_aperture) };
    std::vector<double> rs = {hex_side_in*sin(60*deg), hex_side_out*sin(60*deg) };
    /// Hexagonal pyramid
    Polyhedra shape ( "mypyramid", 6, 30 * deg, 360 * deg, zplanes, rs );



    /// rotation of 90deg around Y axis, to align Z axis of pyramid with X axis of cylinder
    Transform3D pyramidTr ( RotationZYX ( 0, -90. * deg, 0. * deg ), Translation3D ( 0, 0, 0 ) );


    double phistep = hex_aperture*3;
  double zstep =  hex_side*sin(60*deg);





    // place mother volume (vessel)
    Volume motherVol = desc.pickMotherVolume ( det );

//      Volume vesselVol ( Form("%s_v%d",detName.c_str(),1), shape, desc.material ( "Aluminum" ) );
//                 PlacedVolume vesselPV = motherVol.placeVolume ( vesselVol , Translation3D(0,0,-190*cm)*RotationZ (0 ));
//             vesselPV.addPhysVolID ( "system", detID );
//             vesselPV.addPhysVolID ( "module", 0 );
//             DetElement subdet ( detName + std::to_string ( 0 ), detID + 1 + 0 );
//             det.setPlacement ( vesselPV );

//     Transform3D pyrTr1 =  RotationX( -90. * deg) * RotationY( -90. * deg) * Translation3D(0, 0, -190*cm);
    for( int ring = 0; ring<6; ++ring)
    {

        double phi_offset = 0 + 0.5 * phistep * (0 == ring%2);
        double mirror_abs_pos_z = ring * zstep;// - 4*cm * (0 == ring%2);
        Volume vesselVol ( Form("%s_v%d",detName.c_str(),ring), shape, desc.material ( "Aluminum" ) );
            vesselVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d", ring)));

        for ( int phin = 0; phin<29; ++phin ) {
            PlacedVolume vesselPV = motherVol.placeVolume ( vesselVol, Translation3D(0,0,mirror_abs_pos_z)*RotationZ ( phistep * phin + phi_offset )* RotationY ( 90. * deg ) );
            vesselPV.addPhysVolID ( "system", detID );
            vesselPV.addPhysVolID ( "module", phin );
            DetElement subdet ( detName + std::to_string ( phin ), detID + 1 + phin );
            det.setPlacement ( vesselPV );

        }

    }

//     det.setPlacement ( vesselPV );




//     det.setPlacement ( vesselPV );



    return det;
}


static Ref_t create_barrel_cell(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  xml::Component dims = detElem.dimensions();
  DetElement det(detName, detID);
  sens.setType("tracker");

  // read Martin file and store parameters by name in the map

  // mother volume corresponds to the world
  Volume motherVol = desc.pickMotherVolume(det);

  // Vessel, cylindral
  double vessel_outer_r = 210 * cm;
  double vessel_inner_r = 190 * cm;
  double vessel_length = 440 * cm;
  double vessel_wall_thickness = 1.0 * cm;
  if (vessel_outer_r <= vessel_inner_r)
    throw std::runtime_error("Ilegal parameters: vessel_outer_r <= vessel_inner_r");

  // Cooling,
  double cooling_radial_thickness = 1 * cm;

  // Cell parameters
  /// Cell is intersection of hexagonal pyramid and the cylinder
  double hexagon_side_length = 14.815 * cm;
  /// Distance in x-direction
  double zstep = 2 * hexagon_side_length;
  /// Distance in phi angle between cells
  /// since cells are regular hexagons, this distance matches
  /// the angle that one cell covers
  double phistep = 13.333 * deg;
  /// number of repetition of unique cells around the barrel
  int phinmax = 27; // 27;

  // Mirror parameters
  double thickness_sphere(10 * mm);
  double mirror_diameter_safe_shrink = 1*mm;
  double mirror_z_safe_shrink = 1*mm;
  double mirror_z_origin_Martin = vessel_outer_r - vessel_wall_thickness - 37 * cm - mirror_z_safe_shrink;

  // Light sensor parameters
  double sensor_sidex = 8 * cm;
  double sensor_sidey = 8 * cm;
  double sensor_thickness = 0.2 * cm;
  double sensor_z_origin_Martin = vessel_inner_r + vessel_wall_thickness + cooling_radial_thickness;

  // Build cylinder for gas.
  Tube gasvolSolid(vessel_inner_r + vessel_wall_thickness,
                   vessel_outer_r - vessel_wall_thickness,
                   vessel_length);

  // Use regular polyhedra for endcaps cells
  // PolyhedraRegular cellS("aa",6,0,4*cm);

  // Use pyramid for barrel cells
  std::vector<double> zplanes = { vessel_inner_r, vessel_outer_r };
  std::vector<double> rs = { 0.75*hexagon_side_length, 0.95*hexagon_side_length };
  /// Hexagonal pyramid
  Polyhedra shape("mypyramid", 6, 30 * deg, 360 * deg, zplanes, rs);
  /// rotation of 90deg around Y axis, to align Z axis of pyramid with X axis of cylinder
  Transform3D pyramidTr(RotationZYX(0, -90. * deg, 0. * deg), Translation3D(0, 0, 0));

  // Build sensor shape
  Box sensorSol(sensor_sidex / 2, sensor_sidey / 2, sensor_thickness / 2);
  Volume sensorVol(detName + "_sensor", sensorSol, desc.material("Aluminum"));

  // Build the mirror for ncell=1..18
  std::vector<int> ncell_vector = {-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18,
                                   1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
  // std::vector<int> ncell_vector = { /*-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13, -14,-15,-16,-17,-18,*/
  //                                 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
  //                                 };
  ncell_vector = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
  for (auto ncell : ncell_vector)
  {
    // The following line skips even number cells
    // if (!(ncell % 2))
    //   continue;

    // The following line skips odd number cells
    // if ((ncell % 2))
    // continue;

    /// cell shape. Coordinate system still the same as cylinder!
    Volume cellVol(detName + "_cell" + std::to_string(ncell), shape, desc.material("Aluminum"));
    cellVol.setVisAttributes(desc.visAttributes(Form("mirror_vis%d",fabs(ncell))));

    // there is no cell number 0, and cell number 1 do not need to be reflected
    if (0 == ncell || -1 == ncell)
      continue;

    // cells at z
    bool reflect_parameters = false;
    if (0 > ncell)
    {
      ncell *= -1;
      reflect_parameters = true;
    }

    // initialize parameters for creating the mirror
    double center_of_sphere_x(-999.);
    double center_of_sphere_z(-999.);
    double radius_of_sphere(-999.);

    double center_of_sensor_x(-999.);
    double angle_of_sensor(-999.);

    // convert Roger nomenclature (one cell number) to Martin nomenclature (row and col numbers)
    int name_col = ncell / 2;
    int name_row = ncell % 2 ? 1 : 2;



    // position of mirror in cylinder coordinate system
    double mirror_abs_pos_z = name_col * zstep - 0.5 * zstep * (2 == name_row);
    if (reflect_parameters)
      mirror_abs_pos_z *= -1.0;

    // row 2 is shifted half step size
    double phi_offset = 0 + 0.5 * phistep * (2 == name_row);

    for (int phin = 0; phin < phinmax; ++phin)
    {
      PlacedVolume cellPV = motherVol.placeVolume(cellVol, RotationZ(phistep * phin + phi_offset) * Translation3D(0, 0, mirror_abs_pos_z)*pyramidTr);
      cellPV.addPhysVolID("system", detID).addPhysVolID("module", 17 * phin + name_col);
      // create mirrors as separate detectors, so properties can be adjusted lated!
      det.setPlacement(cellPV);
    }
  }

  return det;
}



DECLARE_DETELEMENT ( TRICKY_PYR_TYPE, create_barrel_cell )


